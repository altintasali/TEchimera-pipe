# -----------------------------------------------------------------------------
# Optional TrimGalore! adapter/quality trimming before STAR, mirroring the
# nf-core/rnaseq defaults (Trim Galore! is nf-core's default trimmer, run
# with --paired/--gzip and FastQC executed inside via --fastqc_args).
# Controlled by the `trimming` section of config.yaml; set enabled: false to
# skip this step entirely (STAR then reads the merged/raw fastqs directly).
#
# Two rules because Snakemake requires a rule's output patterns to be static
# and TrimGalore! names paired vs single-end outputs differently. Only one of
# them ever triggers per sample (paired vs single-end is a per-sample
# property), so they can share the shell and helper functions below.
# -----------------------------------------------------------------------------


def _trim_inputs(wildcards):
    d = {"r1": merged_fastq_path(wildcards.sample, 1)}
    if _is_paired(wildcards.sample):
        d["r2"] = merged_fastq_path(wildcards.sample, 2)
    return d


def _trim_reads_param(wildcards):
    d = _trim_inputs(wildcards)
    return d["r1"] + (" " + d["r2"] if "r2" in d else "")


def _trim_paired_flag(wildcards):
    return "--paired" if _is_paired(wildcards.sample) else ""


def _nextseq_param(wildcards):
    n = int(config["trimming"].get("trim_nextseq", 0) or 0)
    return f"--nextseq {n}" if n > 0 else ""


def _trim_outdir(wildcards):
    return "results/trimming"


_TRIM_SHELL = (
    "mkdir -p {params.outdir} && "
    "trim_galore {params.paired} --gzip --cores {threads} "
    "{params.nextseq} {params.extra} "
    "--fastqc_args '-t {threads}' "
    "--basename {wildcards.sample} --output_dir {params.outdir} "
    "{params.reads} "
    "> {log} 2>&1"
)


rule trim_galore_pe:
    # Paired-end: --basename fixes the output names ({sample}_val_1/2.fq.gz)
    # so they don't depend on the input file's basename -- the merged path
    # (results/fastq/) and a raw single-lane path differ, but both must map
    # to the same outputs. The FastQC reports (run inside TrimGalore! via
    # --fastqc_args) are declared as outputs too: the multiqc rule depends on
    # them rather than on the (possibly temp()-deleted) trimmed fastqs.
    input:
        unpack(_trim_inputs),
    output:
        fq1=_maybe_temp("results/trimming/{sample}_val_1.fq.gz", KEEP_TRIMMED_FASTQ),
        fq2=_maybe_temp("results/trimming/{sample}_val_2.fq.gz", KEEP_TRIMMED_FASTQ),
        fq1_fastqc="results/trimming/{sample}_val_1_fastqc.zip",
        fq2_fastqc="results/trimming/{sample}_val_2_fastqc.zip",
    params:
        paired=_trim_paired_flag,
        reads=_trim_reads_param,
        nextseq=_nextseq_param,
        extra=config["trimming"].get("extra", ""),
        outdir=_trim_outdir,
    threads: get_resources("trim_galore_pe")["threads"]
    resources:
        mem_mb=get_resources("trim_galore_pe")["mem_mb"],
        runtime=get_resources("trim_galore_pe")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/trim_galore_pe/{sample}.txt",
    log:
        "results/pipeline_info/logs/trimming/{sample}.log",
    conda:
        TRIM_GALORE_ENV
    shell:
        _TRIM_SHELL


rule trim_galore_se:
    # Single-end: TrimGalore! names the output {sample}_trimmed.fq.gz.
    input:
        unpack(_trim_inputs),
    output:
        fq=_maybe_temp("results/trimming/{sample}_trimmed.fq.gz", KEEP_TRIMMED_FASTQ),
        fastqc="results/trimming/{sample}_trimmed_fastqc.zip",
    params:
        paired=_trim_paired_flag,
        reads=_trim_reads_param,
        nextseq=_nextseq_param,
        extra=config["trimming"].get("extra", ""),
        outdir=_trim_outdir,
    threads: get_resources("trim_galore_se")["threads"]
    resources:
        mem_mb=get_resources("trim_galore_se")["mem_mb"],
        runtime=get_resources("trim_galore_se")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/trim_galore_se/{sample}.txt",
    log:
        "results/pipeline_info/logs/trimming/{sample}.se.log",
    conda:
        TRIM_GALORE_ENV
    shell:
        _TRIM_SHELL
