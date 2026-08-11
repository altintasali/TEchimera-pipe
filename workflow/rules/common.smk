import gzip
import itertools
import os
import re
import tempfile

import pandas as pd
import yaml
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

V = config["versions"]

# -----------------------------------------------------------------------------
# Per-rule compute resources (threads/mem_mb/runtime), loaded from
# workflow/default-config/resources.yaml (built-in defaults) and optionally
# overridden by input/resources.yaml if you create one. A rule (or a
# missing key within its entry) not
# present there falls back to a small conservative default instead of
# failing -- see the "HPC / SLURM" section of the README for how these feed
# into cluster execution.
# -----------------------------------------------------------------------------
RESOURCES = config.get("resources", {})
_RESOURCE_DEFAULTS = {"threads": 1, "mem_mb": 4000, "runtime": 60}


def get_resources(rule_name):
    """Return {threads, mem_mb, runtime} for a rule name."""
    return {**_RESOURCE_DEFAULTS, **RESOURCES.get(rule_name, {})}

# -----------------------------------------------------------------------------
# Load & validate sample sheet
# columns: sample, fastq_1, fastq_2 (optional), strandedness (optional),
#          condition (optional).
#
# Repeated rows with the same `sample` name are allowed (nf-core/rnaseq
# style): they mean that sample's reads are split across lanes/runs and will
# be concatenated into one file per read before alignment (rules/fastq.smk).
# -----------------------------------------------------------------------------
if not os.path.exists(config["samples"]):
    raise WorkflowError(
        f"samples: '{config['samples']}' (from config.yaml, or --configfile/"
        "--config on the command line) does not exist.\n"
        f"Current working directory (relative paths are resolved against "
        f"this): {os.getcwd()}\n"
        f"That path resolves to: {os.path.abspath(config['samples'])}\n"
        "Snakemake resolves relative paths against the directory you *run "
        "the command from*, not the directory workflow/Snakefile lives in "
        "-- cd into the repo root (the directory that directly contains "
        "workflow/, config/, and .tests/) and try again."
    )


def _nonempty(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text or None


raw_samples = pd.read_csv(config["samples"], dtype=str, comment="#")
# jsonschema's "null" type only matches Python None, not pandas/NumPy NaN, so
# empty optional cells (fastq_2, strandedness, condition) need to be
# converted explicitly. astype(object) first is required on pandas>=2's
# str-backed columns, which otherwise silently keep NaN instead of None.
# Note: validate() below re-introduces NaN for empty cells (observed with
# snakemake 9.x), so _nonempty() also treats NaN/NA as empty defensively.
raw_samples = raw_samples.astype(object).where(pd.notnull(raw_samples), None)
validate(raw_samples, schema="../schemas/samples.schema.yaml")

# Group repeated rows for the same sample into one per-sample entry holding a
# list of fastq paths per read. Consistency is enforced eagerly at parse time
# so a typo in a lane-split sample fails fast instead of mid-run.
_rows = {}
for _, row in raw_samples.iterrows():
    sample = _nonempty(row["sample"])
    if sample is None:
        raise WorkflowError(f"{config['samples']}: a row has no sample name.")
    entry = _rows.setdefault(
        sample,
        {"fastq_1": [], "fastq_2": [], "strandedness": None, "condition": None},
    )
    fq1 = _nonempty(row["fastq_1"])
    if fq1 is None:
        raise WorkflowError(f"sample '{sample}' has a row with no fastq_1.")
    entry["fastq_1"].append(fq1)
    fq2 = _nonempty(row["fastq_2"])
    if fq2 is not None:
        entry["fastq_2"].append(fq2)
    for col in ("strandedness", "condition"):
        value = _nonempty(row.get(col))
        if value is None:
            continue
        if entry[col] is not None and entry[col] != value:
            raise WorkflowError(
                f"sample '{sample}' has inconsistent '{col}' values across "
                f"rows ({entry[col]!r} vs {value!r}) in {config['samples']}."
            )
        entry[col] = value

for sample, entry in _rows.items():
    n_r1, n_r2 = len(entry["fastq_1"]), len(entry["fastq_2"])
    if n_r2 and n_r2 != n_r1:
        raise WorkflowError(
            f"sample '{sample}' has {n_r1} fastq_1 files but {n_r2} fastq_2 "
            f"files in {config['samples']} -- every lane/run must provide "
            "both reads (or the sample must be single-end throughout)."
        )

samples = (
    pd.DataFrame([{"sample": s, **e} for s, e in _rows.items()])
    .set_index("sample", drop=False)
    .sort_index()
)

SAMPLES = list(samples["sample"])
HAS_CONDITION = "condition" in samples.columns and samples["condition"].notna().all()


def sample_fastqs(sample, read):
    """Raw fastq paths for a sample's read 1/2, across all lanes/runs."""
    files = samples.loc[sample, f"fastq_{read}"]
    if isinstance(files, list):
        return [f for f in files if f]
    return [files] if files else []

# Warn (not fail) about sample-sheet fastq paths that don't exist yet: a typo
# would otherwise only surface when that rule's job runs -- potentially hours
# into a SLURM queue. Warning keeps pre-download dry-runs and laptops running
# against cluster-only paths working (the sjdb_overhang "auto" check above
# similarly tolerates missing files if an explicit value is given).
_missing_fastqs = [
    f"{sample}: {path}"
    for sample in SAMPLES
    for read in (1, 2)
    for path in sample_fastqs(sample, read)
    if path and not os.path.exists(path)
]
if _missing_fastqs:
    logger.warning(
        f"{len(_missing_fastqs)} fastq file(s) referenced by the sample sheet "
        "do not exist yet:\n  "
        + "\n  ".join(_missing_fastqs)
        + "\nThey must be present before the workflow can run; check the "
        "paths (and that you're in the directory they're relative to)."
    )

# -----------------------------------------------------------------------------
# Reference files: transparently support gzipped fasta/gtf/te_gtf.
# STAR's --genomeFastaFiles/--sjdbGTFfile and the annotation-to-BED rules
# (chimera.smk) all expect plain-text files, so any ref.* path ending in .gz is
# decompressed once via the gunzip_reference rule (ref.smk); everything else
# in the workflow refers to the resolved (always-uncompressed) path below,
# never to config["ref"][...] directly.
#
# Decompressed files land in ref.decompressed_dir. The default is a
# directory under the results tree (results/pipeline_info/ref_decompressed):
# unlike /tmp, that lives on the same shared filesystem the snakemake
# scheduler and every cluster job can see, so gunzip_reference's output is
# always visible -- a node-local /tmp path silently breaks cluster runs with
# "output ... missing locally, parent dir not present" (the file is written
# on the compute node, but the scheduler checks for it on the submission
# node). The files are temp()-wrapped (see the gunzip_reference rule) so they
# are deleted after downstream rules consume them. Set ref.decompressed_dir
# to a scratch dir on shared storage if you would rather keep them.
#
# Because a node-local dir (or an empty value, which would resolve to
# "/<stem>" at the filesystem root) can only break runs that actually
# decompress something, we hard-fail at parse time when gzipped references
# are present -- before any job is submitted -- instead of letting the first
# gunzip_reference job die hours into a queue. Configs with plain refs never
# reach gunzip_reference, so an explicit decompressed_dir is dead config
# there and left alone.
# -----------------------------------------------------------------------------
DECOMPRESS_DIR = config["ref"].get(
    "decompressed_dir", "results/pipeline_info/ref_decompressed"
)

_HAS_GZ_REFS = any(
    str(config["ref"][k]).endswith(".gz") for k in ("fasta", "gtf", "te_gtf")
)
if "decompressed_dir" in config.get("ref", {}) and _HAS_GZ_REFS:
    _local_tmp = tempfile.gettempdir().rstrip(os.sep)
    _node_local = (
        DECOMPRESS_DIR == "/tmp"
        or DECOMPRESS_DIR.startswith("/tmp/")
        or DECOMPRESS_DIR == "/var/tmp"
        or DECOMPRESS_DIR.startswith("/var/tmp/")
        or "$TMPDIR" in DECOMPRESS_DIR
        or (_local_tmp and DECOMPRESS_DIR.startswith(_local_tmp))
    )
    if _node_local or not DECOMPRESS_DIR.strip():
        raise WorkflowError(
            f"Gzipped reference files are used, but ref.decompressed_dir is "
            f"{DECOMPRESS_DIR!r}{' (empty)' if not DECOMPRESS_DIR.strip() else ''} "
            f"-- a node-local temp directory. gunzip_reference writes its "
            f"output there, and on a multi-node cluster each job runs on a "
            f"compute node with its own /tmp that the snakemake scheduler "
            f"cannot see ('MissingOutputException ... (missing locally, parent "
            f"dir not present)'). Delete the ref.decompressed_dir line to use "
            f"the default results/pipeline_info/ref_decompressed (shared "
            f"storage, temp()-cleaned), or point it at a scratch directory on "
            f"shared storage."
        )

REFERENCE_GZ_SOURCES = {}  # decompressed stem -> original .gz path, for gunzip_reference


def _resolve_ref_path(key):
    path = config["ref"][key]
    if str(path).endswith(".gz"):
        stem = os.path.basename(str(path))[:-3]  # strip trailing ".gz"
        decompressed = f"{DECOMPRESS_DIR}/{stem}"
        if stem in REFERENCE_GZ_SOURCES and REFERENCE_GZ_SOURCES[stem] != path:
            raise ValueError(
                f"Two different ref.* .gz files decompress to the same "
                f"filename '{stem}' ({decompressed}) -- "
                "rename one of them so their basenames are unique."
            )
        REFERENCE_GZ_SOURCES[stem] = path
        return decompressed
    return path


FASTA = _resolve_ref_path("fasta")
GTF = _resolve_ref_path("gtf")
TE_GTF = _resolve_ref_path("te_gtf")

# References are required for everything downstream (index, gene model,
# chimera overlap detection), so fail fast at parse time on a missing/typo'd
# path instead of letting the first rule that touches it die hours into a
# queue. For .gz refs we check the .gz source: the decompressed target is a
# temp() rule output that legitimately doesn't exist until gunzip_reference
# runs.
_missing_refs = []
for _key, _resolved in (("fasta", FASTA), ("gtf", GTF), ("te_gtf", TE_GTF)):
    _source = config["ref"][_key]
    _check = _source if str(_source).endswith(".gz") else _resolved
    if not os.path.exists(_check):
        _missing_refs.append(
            f"ref.{_key}: {_check} (exists: False, resolves to: "
            f"{os.path.abspath(_check)})"
        )
if _missing_refs:
    raise WorkflowError(
        "Missing reference file(s) from config.yaml:\n  "
        + "\n  ".join(_missing_refs)
        + "\nFix the paths in config.yaml (they're resolved relative to the "
        "directory you run snakemake from, usually the repo root)."
    )

# -----------------------------------------------------------------------------
# STAR sjdbOverhang: auto-detected from the sample sheet's fastq files
# (max read length - 1, STAR's own recommendation) unless the user pins an
# explicit integer in config.yaml. Detection peeks at the first few reads of
# every fastq_1/fastq_2 file, gzip-aware, so it works whether reads are
# compressed or not.
# -----------------------------------------------------------------------------
def _open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _fastq_read_length(path, n_reads=20):
    """Max read length in the first n_reads of a fastq(.gz) file, or 0 if
    the file doesn't exist / can't be read (e.g. not downloaded yet)."""
    if not path or not os.path.exists(path):
        return 0
    max_len = 0
    try:
        with _open_maybe_gz(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:  # sequence line of each fastq record
                    max_len = max(max_len, len(line.strip()))
                if i >= n_reads * 4:
                    break
    except OSError:
        return 0
    return max_len


def _resolve_sjdb_overhang():
    configured = config["ref"].get("sjdb_overhang", "auto")
    if str(configured).strip().lower() != "auto":
        return int(configured), "config.yaml (explicit)"
    paths = []
    for s in SAMPLES:
        for read in (1, 2):
            paths.extend(sample_fastqs(s, read))
    max_len = max((_fastq_read_length(p) for p in paths), default=0)
    if max_len == 0:
        cwd = os.getcwd()
        checked = "\n".join(
            f"  - {p}  (exists: {os.path.exists(p)}, resolves to: "
            f"{os.path.abspath(p)})"
            for p in paths
        ) or "  (no fastq_1/fastq_2 values found in the sample sheet at all)"
        raise ValueError(
            "ref.sjdb_overhang is 'auto' but no readable fastq file was "
            f"found to auto-detect read length from. Current working "
            f"directory (paths are resolved relative to this): {cwd}\n"
            f"Paths checked from {config['samples']}:\n{checked}\n"
            "Either run Snakemake from the directory these paths are "
            "relative to (usually the repo root), fix the paths in the "
            "sample sheet, or set ref.sjdb_overhang to an explicit integer "
            "in config.yaml to skip auto-detection entirely."
        )
    return max_len - 1, f"auto-detected (max read length {max_len} - 1)"


SJDB_OVERHANG, _SJDB_OVERHANG_SOURCE = _resolve_sjdb_overhang()

# Record how key values were resolved -- a lightweight, always-current log of
# workflow parse-time decisions, separate from the per-rule logs in
# results/pipeline_info/logs/.
os.makedirs("results/pipeline_info/logs", exist_ok=True)
with open("results/pipeline_info/logs/config_resolution.log", "w") as fh:
    fh.write(f"sjdb_overhang = {SJDB_OVERHANG} ({_SJDB_OVERHANG_SOURCE})\n")
    fh.write(f"decompressed reference directory = {DECOMPRESS_DIR}\n")
    if REFERENCE_GZ_SOURCES:
        fh.write("gzipped reference files to be decompressed:\n")
        for stem, gz_path in REFERENCE_GZ_SOURCES.items():
            fh.write(f"  {gz_path} -> {DECOMPRESS_DIR}/{stem}\n")
    else:
        fh.write("no gzipped reference files (fasta/gtf/te_gtf) detected\n")

# -----------------------------------------------------------------------------
# Where the STAR index lives. Defaults to a generated directory under results/
# when star.index is not set (or empty) in config.yaml; an existing index
# directory is honored as-is (rebuilt only when missing or via
# `snakemake -R star_index`).
# -----------------------------------------------------------------------------
STAR_INDEX_DIR = (config["star"].get("index") or "").strip() or "results/star_index"

# -----------------------------------------------------------------------------
# STAR index freshness. star_index (ref.smk) marks its fasta/gtf inputs
# ancient() so a shared/prebuilt index is honored as-is regardless of file
# mtimes -- but that also means a changed reference is silently *not* rebuilt.
# So star_index records a small stamp of everything the index was built with
# (ref paths + size + mtime, sjdb_overhang, STAR version, index_extra) into
# {index}/.refs_used.txt, and we compare it at parse time: a mismatch warns
# that the current config expects a different index and -R star_index is
# needed. An index without a stamp (pre-existing/shared) is left alone.
# -----------------------------------------------------------------------------
def _stat_line(prefix, path):
    try:
        st = os.stat(path)
        return [f"{prefix}_size={st.st_size}", f"{prefix}_mtime={st.st_mtime}"]
    except OSError:
        return [f"{prefix}_size=missing", f"{prefix}_mtime=missing"]


def star_index_stamp():
    """One key=value line per input that determines the STAR index content."""
    lines = []
    for key in ("fasta", "gtf"):
        src = str(config["ref"][key])
        lines.append(f"{key}={src}")
        lines.extend(_stat_line(key, src))
    lines.append(f"sjdb_overhang={SJDB_OVERHANG}")
    lines.append(f"star={config['versions']['star']}")
    lines.append(f"index_extra={config['star'].get('index_extra', '')}")
    return "\n".join(lines) + "\n"


_STAR_INDEX_MANIFEST = os.path.join(STAR_INDEX_DIR, ".refs_used.txt")
if os.path.isdir(STAR_INDEX_DIR) and os.path.isfile(_STAR_INDEX_MANIFEST):
    try:
        with open(_STAR_INDEX_MANIFEST) as fh:
            _recorded = fh.read()
    except OSError:
        _recorded = ""
    if _recorded != star_index_stamp():
        logger.warning(
            f"The STAR index at {STAR_INDEX_DIR} was built for a different "
            "reference setup than config.yaml now specifies (its "
            ".refs_used.txt stamp no longer matches). If you changed the "
            "reference fasta/GTF, sjdb_overhang, or STAR version since, "
            "rebuild it with `snakemake -R star_index` -- otherwise the old "
            "index is silently reused. If the index is deliberately shared, "
            "run `snakemake -R star_index` once to re-record the current "
            "reference into it and silence this warning."
        )

# -----------------------------------------------------------------------------
# Per-sample strandedness resolution.
# The sample sheet's optional "strandedness" column (nf-core/rnaseq
# vocabulary: auto/forward/reverse/unstranded, so an nf-core samplesheet can
# be used as-is) takes priority per-sample; if empty/absent for a sample, it
# defaults to "auto". Resolved and validated once, eagerly, so a typo fails
# fast at parse time rather than mid-run.
# -----------------------------------------------------------------------------
STRANDEDNESS_ALIASES = {
    "auto": "auto",
    "forward": "forward",
    "reverse": "reverse",
    "unstranded": "no",  # nf-core vocabulary -> resolved no/forward/reverse
    "no": "no",  # "no" accepted directly too (the classic library value)
}


def _normalize_stranded(raw, context=""):
    key = str(raw).strip().lower()
    if key not in STRANDEDNESS_ALIASES:
        raise ValueError(
            f"Unrecognized strandedness value '{raw}'{context}. Expected one "
            "of: auto, forward, reverse, unstranded (nf-core vocabulary), or "
            "no (the classic unstranded value)."
        )
    return STRANDEDNESS_ALIASES[key]


# Default when a sample sheet has no "strandedness" column, or leaves it
# blank for a given sample: auto-detect via RSeQC.
DEFAULT_STRANDEDNESS_MODE = "auto"


def _sample_effective_mode(sample):
    if "strandedness" in samples.columns:
        val = samples.loc[sample, "strandedness"]
        if pd.notna(val) and str(val).strip() != "":
            return _normalize_stranded(val, f" for sample '{sample}' in {config['samples']}")
    return DEFAULT_STRANDEDNESS_MODE


SAMPLE_STRANDED_MODE = {s: _sample_effective_mode(s) for s in SAMPLES}
AUTO_SAMPLES = [s for s in SAMPLES if SAMPLE_STRANDED_MODE[s] == "auto"]


def _is_paired(sample):
    """A sample is paired-end if it has at least one non-empty fastq_2 value.

    This is resolved per-sample (not globally), so a single sample sheet can
    freely mix paired-end and single-end samples. Every lane/run of one
    sample must agree (enforced while loading the sample sheet above).
    """
    return len(sample_fastqs(sample, 2)) > 0


# -----------------------------------------------------------------------------
# Fastq preparation pipeline: raw lanes -> (merge) -> (trim) -> STAR.
# Samples sequenced across multiple lanes/runs are concatenated by
# rules/fastq.smk (cat_fastq) into results/fastq/. If trimming is enabled,
# rules/trimming.smk (trim_galore, nf-core/rnaseq defaults) then runs on the
# merged files and STAR consumes the trimmed output; otherwise STAR reads the
# merged (or raw, for single-lane) files directly.
# -----------------------------------------------------------------------------
TRIM_ENABLED = bool(config.get("trimming", {}).get("enabled", True))

# TrimGalore! always appends _trimmed (single-end) or _val_1/_val_2 (paired)
# to the *input* basename, and its --basename normalization only strips a
# single "_trimmed"/"_val_1" suffix. Feeding it an already-trimmed fastq
# therefore produces a no-op rename like foo_trimmed_trimmed_trimmed.fq.gz
# while the declared results/trimming/<sample>_trimmed.fq.gz output is never
# created (or a stale file is silently reused). Fail fast at parse time
# instead of surfacing a confusing MissingOutputException mid-run.
if TRIM_ENABLED:
    _TRIMMED_FASTQ_RE = re.compile(
        r"(_trimmed|_val_[12]|\.trimmed)(\.fq|\.fastq)(\.gz)?$", re.IGNORECASE
    )
    _offenders = [
        f"{sample}: {path}"
        for sample in SAMPLES
        for read in (1, 2)
        for path in sample_fastqs(sample, read)
        if _TRIMMED_FASTQ_RE.search(os.path.basename(str(path)))
    ]
    if _offenders:
        raise WorkflowError(
            "Trimming is enabled (trimming.enabled), but the sample sheet "
            "points at fastqs that already look trimmed:\n  "
            + "\n  ".join(_offenders)
            + "\nFeeding already-trimmed reads to TrimGalore! produces a "
            "no-op rename and breaks the pipeline. Point the sample sheet at "
            "untrimmed reads, or set `trimming.enabled: false` in config.yaml "
            "to use these files as-is."
        )

# Whether the merged (lane-concatenated) and trimmed fastq files are kept
# after downstream rules (STAR/alignment) are done with them. Set either to
# false in the `outputs` config section and those files are marked temp() --
# Snakemake deletes them after the last consumer finishes, saving disk on big
# runs at the cost of re-running the merge/trim steps if they're needed again.
KEEP_MERGED_FASTQ = bool(config.get("outputs", {}).get("keep_merged_fastq", True))
KEEP_TRIMMED_FASTQ = bool(config.get("outputs", {}).get("keep_trimmed_fastq", True))


def _maybe_temp(path, keep):
    """A rule output path, wrapped in temp() unless the user opted to keep
    it via the `outputs` config section. Snakemake allows temp() on wildcard
    patterns, so this is applied at parse time (the keep flag is a constant)."""
    if keep:
        return path
    return temp(path)


# The cat_fastq rule's output pattern, temp()-wrapped when the merged fastqs
# are not kept (see the `outputs` config section).
MERGED_FASTQ_OUTPUT = _maybe_temp(
    "results/fastq/{sample}_R{read}.fastq.gz", KEEP_MERGED_FASTQ
)


# Directory for STAR's per-run scratch files (the _STARtmp directory), outside
# the results tree. See the star_align rule in align.smk.
STAR_TMPDIR = config["star"].get("tmpdir", "/tmp")


def merged_fastq_path(sample, read):
    """Path of the concatenated fastq for a sample's read (results/fastq/).

    Samples with a single lane/run have no merge step: their raw path is
    returned directly so no unnecessary copy is made. Multi-lane samples are
    always merged to a .gz file (plain-text lanes are gzipped on the fly by
    the cat_fastq rule), so STAR's --readFilesCommand zcat applies uniformly.
    """
    files = sample_fastqs(sample, read)
    if len(files) <= 1:
        return files[0] if files else None
    return f"results/fastq/{sample}_R{read}.fastq.gz"


def star_fastq(sample, read):
    """The fastq path STAR aligns for a sample's read (after merge/trim)."""
    base = merged_fastq_path(sample, read)
    if base is None:
        return None
    if TRIM_ENABLED:
        if _is_paired(sample):
            return f"results/trimming/{sample}_val_{read}.fq.gz"
        return f"results/trimming/{sample}_trimmed.fq.gz"
    return base


def star_input(wildcards):
    """Return {fq1: ..., [fq2: ...]} paths for a sample."""
    d = {"fq1": star_fastq(wildcards.sample, 1)}
    if _is_paired(wildcards.sample):
        d["fq2"] = star_fastq(wildcards.sample, 2)
    return d


def star_reads_param(wildcards):
    """Space-separated --readFilesIn argument (fq1 [fq2])."""
    d = star_input(wildcards)
    reads = [d["fq1"]]
    if "fq2" in d:
        reads.append(d["fq2"])
    return " ".join(reads)


def star_read_command_param(wildcards):
    """--readFilesCommand zcat if the fastq files are gzip-compressed."""
    d = star_input(wildcards)
    if str(d["fq1"]).endswith(".gz"):
        return "--readFilesCommand zcat"
    return ""


def strandedness_input(wildcards):
    """Dependency on this sample's auto-detected strandedness call.

    Returns [] (no dependency) when this sample's effective mode is fixed
    (from the sample sheet's "strandedness" column or the config default),
    since no RSeQC/auto-detection step is needed in that case.
    """
    if SAMPLE_STRANDED_MODE[wildcards.sample] == "auto":
        return f"results/rseqc/{wildcards.sample}_strandedness.txt"
    return []


def get_strandedness_param(wildcards, input):
    """Resolve the --stranded value (no/forward/reverse) for a sample.

    If this sample's effective mode is "auto", read the value that
    workflow/scripts/determine_strandedness.py determined from the RSeQC
    infer_experiment.py output (the file is guaranteed to already exist
    because it is declared as a rule input alongside this params function).
    Otherwise, the fixed value (sample sheet override or config default) is
    used directly.
    """
    mode = SAMPLE_STRANDED_MODE[wildcards.sample]
    if mode == "auto":
        with open(input.strandedness) as fh:
            return fh.read().strip()
    return mode


# -----------------------------------------------------------------------------
# The "condition" column of the sample sheet (parsed above as HAS_CONDITION)
# drives the chimera sample-QC coloring (PCA / sample clustering,
# sample_qc.smk) and your own downstream differential analysis of the chimera
# count matrix -- the pipeline deliberately does no automatic DE on chimera
# counts (they are sparse / zero-inflated); it ships the counts matrix and
# per-junction annotations instead.
# -----------------------------------------------------------------------------


def all_trim_outputs():
    """One trimmed fastq path per sample (only when trimming is enabled),
    used by the `trimming_only` convenience target (Snakefile)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1.fq.gz"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed.fq.gz"
        )
        for s in SAMPLES
    ]


def all_fastqc_reports():
    """One FastQC .zip report path per trimmed sample (only when trimming is
    enabled), used by the multiqc rule so it scans results/trimming/ for the
    FastQC + TrimGalore! reports without depending on the trimmed fastq files
    themselves (which may be temp()-deleted after alignment)."""
    if not TRIM_ENABLED:
        return []
    return [
        (
            f"results/trimming/{s}_val_1_fastqc.zip"
            if _is_paired(s)
            else f"results/trimming/{s}_trimmed_fastqc.zip"
        )
        for s in SAMPLES
    ]


def all_raw_fastqc_reports():
    """FastQC .zip report paths for the raw (merged) input fastqs, one per
    sample/read. Run unconditionally so the MultiQC report always covers the
    untrimmed input, regardless of the optional `trimming` step."""
    return [
        f"results/fastqc/raw/{s}_R{read}_fastqc.zip"
        for s in SAMPLES
        for read in (1, 2)
        if read == 1 or _is_paired(s)
    ]


def all_benchmark_files():
    """Every benchmark file this configuration will produce, so the
    benchmark_summary rule (qc.smk) aggregates exactly the rules that ran
    into the MultiQC resource-usage section. Kept in lockstep with the
    `benchmark:` declarations in the rules -- only the conditional ones
    (merging, trimming, gunzip, RSeQC bed12 conversion, strandedness
    auto-detection) need to be gated here. The multiqc and benchmark_summary
    rules' own benchmarks are deliberately excluded to avoid a cyclic
    dependency (their resource use is negligible)."""
    files = [
        "results/pipeline_info/benchmarks/software_versions/software_versions.txt",
        "results/pipeline_info/benchmarks/star_index/star_index.txt",
    ]
    # BED12 gene-model conversion only runs when RSeQC strandedness
    # auto-detection actually needs it.
    if AUTO_SAMPLES:
        files += [
            "results/pipeline_info/benchmarks/gtf_to_genepred/gtf_to_genepred.txt",
            "results/pipeline_info/benchmarks/genepred_to_bed12/genepred_to_bed12.txt",
        ]
    files += [
        "results/pipeline_info/benchmarks/annotation_to_bed/annotation_to_bed.txt",
        "results/pipeline_info/benchmarks/chimera_counts/chimera_counts.txt",
    ]
    for stem in REFERENCE_GZ_SOURCES:
        files.append(f"results/pipeline_info/benchmarks/gunzip_reference/{stem}.txt")
    for s in SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/star_align/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_sort/{s}.txt",
            f"results/pipeline_info/benchmarks/samtools_index/{s}.txt",
            f"results/pipeline_info/benchmarks/parse_chimeric_junctions/{s}.txt",
            f"results/pipeline_info/benchmarks/junction_qc/{s}.txt",
            f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R1.txt",
        ]
        if WRITE_IGV_BED:
            files.append(f"results/pipeline_info/benchmarks/chimera_igv_bed/{s}.txt")
        if _is_paired(s):
            files.append(f"results/pipeline_info/benchmarks/fastqc_raw/{s}_R2.txt")
        # Lane-merging (cat_fastq) only for samples with multiple lanes/runs.
        for read in (1, 2):
            if (read == 1 or _is_paired(s)) and len(sample_fastqs(s, read)) > 1:
                files.append(
                    f"results/pipeline_info/benchmarks/cat_fastq/{s}_R{read}.txt"
                )
        if TRIM_ENABLED:
            files.append(
                f"results/pipeline_info/benchmarks/"
                f"{'trim_galore_pe' if _is_paired(s) else 'trim_galore_se'}/{s}.txt"
            )
    for s in AUTO_SAMPLES:
        files += [
            f"results/pipeline_info/benchmarks/rseqc_infer_experiment/{s}.txt",
            f"results/pipeline_info/benchmarks/determine_strandedness/{s}.txt",
        ]
    if WRITE_COUNTS:
        transform = CHIMERA_QC["pca_transform"]
        files += [
            f"results/pipeline_info/benchmarks/sample_qc_transform/{transform}.txt",
            f"results/pipeline_info/benchmarks/sample_qc/{transform}.txt",
        ]
    return sorted(set(files))


def allocated_resources_by_rule():
    """{rule: {"threads", "mem_mb"}} for every rule that has benchmark files
    (the benchmark_summary rule's input), read from resources.yaml -- the
    per-job allocation against which the benchmark_summary script computes
    CPU/RAM efficiency. Iterated in Snakemake's own rule order (workflow.rules,
    an OrderedDict) so the resource-usage table lists rules in the order they
    appear in the workflow, not alphabetically."""
    benchmark_rules = {
        os.path.basename(os.path.dirname(path)) for path in all_benchmark_files()
    }
    out = {}
    for r in workflow.rules:
        if r.name in benchmark_rules:
            out[r.name] = get_resources(r.name)
    return out


# -----------------------------------------------------------------------------
# Conda environments for every tool in the workflow, generated from
# config["versions"]. Editing a version string there and re-running is all
# that's needed -- Snakemake/conda resolves and downloads the matching build
# automatically (via `--sdm conda` / `--use-conda`), no manual install
# required. This intentionally does NOT rely on snakemake-wrapper-bundled
# environments, because a wrapper's pinned tag moves all of its tools'
# versions together -- it can't give you e.g. STAR 2.7.11b + samtools 1.21
# independently if that particular wrapper tag happens to pin a different
# samtools. Generating one small env per tool from explicit version strings
# is the only way to match an arbitrary external reference (e.g. a specific
# nf-core/rnaseq release's tool versions) exactly.
# -----------------------------------------------------------------------------
# Absolute path: rules in this repo are included from workflow/Snakefile, and
# Snakemake resolves relative paths declared inside an included file against
# that file's directory (workflow/rules/), not the run directory -- so a
# relative "workflow/envs/generated" here would make rule env: references point
# at workflow/rules/workflow/envs/generated/ (which never exists). The env
# files are written at parse time relative to the process CWD, so an absolute
# path is the one form both sides agree on.
GENERATED_ENV_DIR = os.path.abspath("workflow/envs/generated")
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)

# Absolute path to workflow/scripts: shell directives that run the workflow's
# own python scripts need a path that resolves the same way from the run
# directory (the repo root) as the included rules file is parsed from.
SCRIPTS_DIR = os.path.abspath("workflow/scripts")


def _write_env(name, dependencies, pip_dependencies=None):
    dependencies = list(dependencies)
    if pip_dependencies:
        dependencies.append({"pip": list(pip_dependencies)})
    path = f"{GENERATED_ENV_DIR}/{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(
            {"channels": ["bioconda", "conda-forge"], "dependencies": dependencies},
            fh,
            sort_keys=False,
        )
    return path


STAR_ENV = _write_env("star", [f"star={V['star']}"])
SAMTOOLS_ENV = _write_env("samtools", [f"samtools={V['samtools']}"])
# trim-galore brings cutadapt (its core dependency) along automatically;
# fastqc is added explicitly because trim_galore's --fastqc_args (nf-core/
# rnaseq default) shells out to it, and its reports feed the MultiQC report.
TRIM_GALORE_ENV = _write_env(
    "trim_galore", [f"trim-galore={V['trim_galore']}", f"fastqc={V['fastqc']}"]
)
# Standalone FastQC env for the always-on raw FastQC (qc.smk), independent of
# the optional trimming step.
FASTQC_ENV = _write_env("fastqc", [f"fastqc={V['fastqc']}"])
# python>=3.9 floor: without it, conda's solver can backtrack all the way to
# an ancient RSeQC/MultiQC build (seen in practice: Python 3.6, from ~2021)
# to find something that resolves at all -- and those old builds pull in a
# pysam/htslib linked against OpenSSL 1.0, which doesn't exist on modern
# systems (`ImportError: libcrypto.so.1.0.0: cannot open shared object
# file`). Pinning a modern floor forces the solver toward current,
# self-consistent builds instead.
RSEQC_ENV = _write_env("rseqc", [f"rseqc={V['rseqc']}", "python>=3.9"])
MULTIQC_ENV = _write_env("multiqc", [f"multiqc={V['multiqc']}", "python>=3.9"])

# Chimera sample-QC (PCA / sample clustering) runs in R with DESeq2
# (nf-core/rnaseq style); deseq2 + r-base come from conda. A python>=3.9
# floor guards the solver against ancient, mutually-incompatible builds.
CHIMERA_QC_ENV = _write_env(
    "chimera_qc",
    [
        f"bioconductor-deseq2={V['deseq2']}",
        f"r-base={V['r_base']}",
        "r-pheatmap",
        "python>=3.9",
    ],
)

UCSC_TOOLS_ENV = _write_env(
    "ucsc_tools",
    [
        f"ucsc-gtftogenepred={V['ucsc_gtftogenepred']}",
        f"ucsc-genepredtobed={V['ucsc_genepredtobed']}",
    ],
)
