#!/usr/bin/env python3
"""Deterministically generates the tiny synthetic dataset under .tests/,
used by config/test.yaml for a fast end-to-end smoke test of the workflow
(mirrors nf-core's `-profile test` pattern).

Not part of the workflow DAG -- run manually if you ever need to
regenerate/resize the test data:

    python3 .tests/generate_test_data.py

Layout produced:
    .tests/resources/genome.fa           one 50kb contig "chrT"
    .tests/resources/genome.gtf          100 genes, "chrT" 1001-31000 (+)
    .tests/resources/te_annotation.gtf   one TE,   "chrT" 40001-46000 (+)
    .tests/reads/{sample}_R{1,2}.fastq.gz            paired-end samples
    .tests/reads/{sample}_R1.fastq.gz                single-end samples
    .tests/reads/treatment_rep1_L{001,002}_R{1,2}.fastq.gz   (paired-end, lane-split)
    .tests/reads/control_rep3_L{001,002}_R1.fastq.gz         (single-end, lane-split)
    .tests/samples.csv

Single-end samples (control_rep3, treatment_rep3) deliberately sit in the
same two conditions as the paired-end samples so the test run exercises the
workflow's per-sample single/paired-end branching -- trimming (trim_galore
vs trim_galore_se), STAR --readFilesIn, RSeQC strandedness auto-detection
(which supports single-end BAMs), and the chimera screen's per-sample
junction parsing. Every sample also carries a small number of SPLIT READS
(gene segment + TE segment in one read, see CHIMERA_READS) so STAR emits
real gene-TE chimeric junctions and the chimera annotation / counts /
sample-QC stages run on non-empty data instead of silently passing on
zero-event inputs.

Genome size note: earlier versions of this script used a 10kb genome, which
turned out to sit in STAR's known small-genome crash zone (custom
references under ~10-20kb can trigger "double free or corruption" while
writing the SAindex, even with --genomeSAindexNbases correctly downscaled --
see e.g. alexdobin/STAR issues #2158, #965, #1472, #350). 50kb plus an
explicit --genomeChrBinNbits override (config/test.yaml) avoids that.

DESeq2 sizing note: the count design matters. A contrast with a handful of
genes whose per-gene dispersions all sit near the minimum makes DESeq2's
mean-dispersion curve fit throw ("all gene-wise dispersion estimates are
within 2 orders of magnitude from the minimum value"). So the synthetic data
uses (a) 100 genes, (b) base expression spread log-uniformly over 2..500
counts, and (c) within-group (between-replicate) multipliers of 0.5x/1.5x
(and 1.0x/3.0x for treatment, which is 2x control on average). That spread
keeps the gene-wise dispersion estimates spanning well over 2 orders of
magnitude, which is what DESeq2's fit requires -- verified empirically.
"""
import gzip
import io
import math
import random
from pathlib import Path

HERE = Path(__file__).parent
GENOME_LEN = 50_000
CONTIG = "chrT"
N_GENES = 100
GENE_LEN = 300  # per-gene length; must stay > FRAGMENT_LEN for sampling
GENE_START = 1_001  # 1-based, inclusive start of the gene block
GENE_END = GENE_START + N_GENES * GENE_LEN - 1  # 1-based, inclusive
TE_START, TE_END = 40_001, 46_001
READ_LEN = 50
FRAGMENT_LEN = 200  # insert size for paired-end reads

# Log-uniform base read-pair counts across the 100 genes (2 .. 500).
BASE_MIN, BASE_MAX = 2, 500

# Per-sample multiplier applied to every gene's base count. Treatment is 2x
# control on average (1.0/2.0/3.0 vs 0.5/1.0/1.5); the 0.5x/1.5x within-group
# spread is what keeps DESeq2's dispersion fit from degenerating (see
# docstring). The *_rep3 samples are single-end (see SINGLE_END_SAMPLES).
SAMPLE_MULT = {
    "treatment_rep1": 1.0,
    "treatment_rep2": 3.0,
    "control_rep1": 0.5,
    "control_rep2": 1.5,
    "treatment_rep3": 2.0,
    "control_rep3": 1.0,
}

# The TE goes the other direction (up in control) so both features aren't
# perfectly confounded; magnitudes just need to be non-trivial.
SAMPLE_TE_PAIRS = {
    "treatment_rep1": 15,
    "treatment_rep2": 25,
    "treatment_rep3": 20,
    "control_rep1": 60,
    "control_rep2": 90,
    "control_rep3": 75,
}

# Samples sequenced single-end: one read per fragment (the 5' read) instead
# of a read pair, and only fastq_1 entries in samples.csv. Mixed into the
# same conditions as the paired-end samples so one test run covers both.
SINGLE_END_SAMPLES = {"control_rep3", "treatment_rep3"}

# Samples whose reads are written as two lane files (nf-core/rnaseq style)
# instead of one, exercising the workflow's lane/run merging (cat_fastq):
# treatment_rep1 is written as *_L001_* and *_L002_* and listed as two rows
# (same sample name) in .tests/samples.csv -- and control_rep3 exercises the
# same merge for a single-end sample (R1 only).
LANE_SPLIT_SAMPLES = {"treatment_rep1", "control_rep3"}
LANE_NAMES = ("L001", "L002")

# Chimeric reads. The samples above only ever sample reads from within a
# single locus, so STAR finds no gene-TE chimeric junctions and the chimera
# screen/QC stages run on empty inputs. To exercise them, every sample gets a
# handful of SPLIT READS (one read carrying a gene segment + a TE segment),
# which STAR's chimeric-alignment detection reports as junctions. Two
# junctions, both present in all six samples (with per-sample count spread so
# the DESeq2 QC view has something to show):
#   gene_to_te  -- read = gene 25bp + TE 25bp
#   te_to_gene  -- read = TE 25bp + gene 25bp
# Values are (gene_to_te reads, te_to_gene reads) per sample. Both segments
# are 25bp >= the test config's --chimSegmentMin 12 / overhang 12, and both
# breakpoints sit well inside the gene/TE features (breakpoint_tolerance=0).
CHIMERA_GENE_INDEX = 49  # GENE50, mid-block
CHIMERA_READS = {
    "treatment_rep1": (6, 4),
    "treatment_rep2": (12, 8),
    "treatment_rep3": (9, 6),
    "control_rep1": (3, 5),
    "control_rep2": (8, 7),
    "control_rep3": (5, 9),
}
CHIMERA_SEGMENT = 25

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def make_genome(seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))


def write_fasta(path, name, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def write_gtf(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")


def sample_read_pairs(genome, region_start, region_end, n_pairs, seed):
    """Simulate n_pairs FR paired-end reads (exact substrings, no errors)
    from within [region_start, region_end) of the + strand of `genome`
    (0-based half-open)."""
    rng = random.Random(seed)
    pairs = []
    max_start = region_end - FRAGMENT_LEN
    if max_start <= region_start:
        raise ValueError("region too small for the requested fragment length")
    for _ in range(n_pairs):
        frag_start = rng.randint(region_start, max_start)
        frag = genome[frag_start : frag_start + FRAGMENT_LEN]
        r1 = frag[:READ_LEN]
        r2 = revcomp(frag[-READ_LEN:])
        pairs.append((r1, r2))
    return pairs


def sample_single_reads(genome, region_start, region_end, n_reads, seed):
    """Simulate n_reads single-end reads (the 5' read of a fragment, exact
    substrings, no errors) from within [region_start, region_end) of the +
    strand of `genome` (0-based half-open). Same fragment sampling as
    sample_read_pairs so single-end counts stay comparable to paired-end
    (one read per fragment in both cases)."""
    rng = random.Random(seed)
    reads = []
    max_start = region_end - FRAGMENT_LEN
    if max_start <= region_start:
        raise ValueError("region too small for the requested fragment length")
    for _ in range(n_reads):
        frag_start = rng.randint(region_start, max_start)
        reads.append(genome[frag_start : frag_start + READ_LEN])
    return reads


def write_fastq_gz(path, reads, name_prefix, mate=1):
    # mtime=0 keeps the gzip header stable so regeneration is byte-identical.
    with gzip.GzipFile(path, "wb", mtime=0) as gz:
        with io.TextIOWrapper(gz, encoding="ascii") as fh:
            for i, seq in enumerate(reads):
                qual = "I" * len(seq)
                name = f"{name_prefix}_{i}"
                if mate is not None:
                    name += f"/{mate}"
                fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def base_counts():
    """n_genes log-uniform integer base counts in [BASE_MIN, BASE_MAX]."""
    lo, hi = math.log(BASE_MIN), math.log(BASE_MAX)
    return [
        round(math.exp(lo + (hi - lo) * i / (N_GENES - 1))) for i in range(N_GENES)
    ]


def gene_region(i):
    """0-based half-open sampling region for gene i (see gene_region0 docstring)."""
    start0 = GENE_START - 1 + i * GENE_LEN
    return start0, start0 + GENE_LEN


def chimera_reads(genome, n_gt, n_tg, seed):
    """Split reads for the chimera fixture: n_gt gene+TE and n_tg TE+gene
    50bp reads. The breakpoints land in the inner half of gene CHIMERA_GENE
    and the inner half of the TE (250bp clear of either end), so they sit
    comfortably inside the GTF features (breakpoint_tolerance=0 in the test
    config)."""
    g0 = GENE_START - 1 + CHIMERA_GENE_INDEX * GENE_LEN
    gene_window = (g0 + GENE_LEN // 4, g0 + 3 * GENE_LEN // 4)
    te_window = (TE_START - 1 + 250, TE_END - 251)
    rng = random.Random(seed)
    reads = []
    for _ in range(n_gt):
        g = rng.randint(*gene_window)
        t = rng.randint(*te_window)
        reads.append(
            genome[g : g + CHIMERA_SEGMENT] + genome[t : t + CHIMERA_SEGMENT]
        )
    for _ in range(n_tg):
        t = rng.randint(*te_window)
        g = rng.randint(*gene_window)
        reads.append(
            genome[t : t + CHIMERA_SEGMENT] + genome[g : g + CHIMERA_SEGMENT]
        )
    return reads


def main():
    resources_dir = HERE / "resources"
    reads_dir = HERE / "reads"
    resources_dir.mkdir(parents=True, exist_ok=True)
    reads_dir.mkdir(parents=True, exist_ok=True)

    genome = make_genome(seed=1)
    write_fasta(resources_dir / "genome.fa", CONTIG, genome)

    gene_rows = []
    for i in range(N_GENES):
        gstart = GENE_START + i * GENE_LEN
        gend = gstart + GENE_LEN - 1
        name = f"GENE{i + 1}"
        gene_rows.append(
            (
                CONTIG,
                "test",
                "exon",
                str(gstart),
                str(gend),
                ".",
                "+",
                ".",
                f'gene_id "{name}"; transcript_id "{name}"; gene_name "{name}";',
            )
        )
    write_gtf(resources_dir / "genome.gtf", gene_rows)

    # TEtranscripts' curated TE GTF format: "exon" feature, gene_id/
    # transcript_id/family_id/class_id attributes (matches the format of
    # the authors' official prebuilt TE GTFs, e.g. mm10/hg38 rmsk.gtf).
    write_gtf(
        resources_dir / "te_annotation.gtf",
        [
            (
                CONTIG,
                "test",
                "exon",
                str(TE_START),
                str(TE_END),
                ".",
                "+",
                ".",
                'gene_id "TE1"; transcript_id "TE1_dup1"; '
                'family_id "TE1fam"; class_id "TE1class";',
            )
        ],
    )

    # 6 samples, 2 conditions x 3 reps: control_rep1/2 + treatment_rep1/2 are
    # paired-end, control_rep3/treatment_rep3 are single-end. Read depths
    # follow the design in the module docstring so DESeq2's dispersion fit
    # works (100 genes, log-uniform counts, within-group spread).
    base = base_counts()
    te_region = (TE_START - 1, TE_END)  # 1-based -> 0-based half-open

    rows = ["sample,fastq_1,fastq_2,strandedness,condition"]
    for sidx, (sample, mult) in enumerate(SAMPLE_MULT.items()):
        is_se = sample in SINGLE_END_SAMPLES
        r1_reads, r2_reads = [], []
        for i in range(N_GENES):
            n_reads = max(1, round(base[i] * mult))
            if is_se:
                r1_reads.extend(
                    sample_single_reads(
                        genome, *gene_region(i), n_reads, seed=sidx * 10_000 + i
                    )
                )
            else:
                for r1, r2 in sample_read_pairs(
                    genome, *gene_region(i), n_reads, seed=sidx * 10_000 + i
                ):
                    r1_reads.append(r1)
                    r2_reads.append(r2)
        if is_se:
            r1_reads.extend(
                sample_single_reads(
                    genome, *te_region, SAMPLE_TE_PAIRS[sample],
                    seed=sidx * 10_000 + 9_999,
                )
            )
        else:
            for r1, r2 in sample_read_pairs(
                genome, *te_region, SAMPLE_TE_PAIRS[sample],
                seed=sidx * 10_000 + 9_999,
            ):
                r1_reads.append(r1)
                r2_reads.append(r2)

        # Chimeric split reads (see the CHIMERA_READS comment). For paired-end
        # samples the split reads go in R1 and a 50bp mate is sampled from the
        # same gene so R1/R2 stay paired; STAR still splits R1 across loci.
        n_gt, n_tg = CHIMERA_READS[sample]
        chim = chimera_reads(
            genome, n_gt, n_tg, seed=sidx * 10_000 + 5_000
        )
        if is_se:
            r1_reads.extend(chim)
        else:
            chim_rng = random.Random(sidx * 10_000 + 5_001)
            g0 = GENE_START - 1 + CHIMERA_GENE_INDEX * GENE_LEN
            for split in chim:
                r1_reads.append(split)
                g = chim_rng.randint(g0 + GENE_LEN // 4, g0 + 3 * GENE_LEN // 4)
                r2_reads.append(
                    revcomp(genome[g : g + READ_LEN])
                )

        if sample in LANE_SPLIT_SAMPLES:
            # Write the sample's reads as two lane files (paired-end reads
            # must stay paired, so split the aligned lists by index; the
            # single-end control_rep3 writes R1-only lanes). The deterministic
            # test sample has far more than two reads, so both lanes end up
            # non-empty.
            half = len(r1_reads) // 2
            for lane, lane_r1 in (
                (LANE_NAMES[0], r1_reads[:half]),
                (LANE_NAMES[1], r1_reads[half:]),
            ):
                fq1 = reads_dir / f"{sample}_{lane}_R1.fastq.gz"
                write_fastq_gz(fq1, lane_r1, f"{sample}_{lane}", mate=None if is_se else 1)
                if is_se:
                    rows.append(
                        f"{sample},.tests/reads/{fq1.name},,auto,"
                        f"{'treatment' if 'treatment' in sample else 'control'}"
                    )
                else:
                    fq2 = reads_dir / f"{sample}_{lane}_R2.fastq.gz"
                    write_fastq_gz(fq2, r2_reads[:half] if lane == LANE_NAMES[0] else r2_reads[half:], f"{sample}_{lane}")
                    rows.append(
                        f"{sample},.tests/reads/{fq1.name},.tests/reads/{fq2.name},auto,"
                        f"{'treatment' if 'treatment' in sample else 'control'}"
                    )
        else:
            fq1 = reads_dir / f"{sample}_R1.fastq.gz"
            write_fastq_gz(fq1, r1_reads, sample, mate=None if is_se else 1)
            if is_se:
                rows.append(
                    f"{sample},.tests/reads/{sample}_R1.fastq.gz,,auto,"
                    f"{'treatment' if 'treatment' in sample else 'control'}"
                )
            else:
                fq2 = reads_dir / f"{sample}_R2.fastq.gz"
                write_fastq_gz(fq2, r2_reads, sample)
                rows.append(
                    f"{sample},.tests/reads/{sample}_R1.fastq.gz,"
                    f".tests/reads/{sample}_R2.fastq.gz,auto,"
                    f"{'treatment' if 'treatment' in sample else 'control'}"
                )

    (HERE / "samples.csv").write_text("\n".join(rows) + "\n")
    print("Wrote test data to", HERE)

    # STAR's own formulas (see the STAR manual, section 2.2.5) -- used to
    # compute config/test.yaml's star.index_extra. Re-run this script and
    # copy the printed values there if you change GENOME_LEN.
    sa_index_nbases = min(14, int(math.log2(GENOME_LEN) / 2 - 1))
    chr_bin_nbits = min(18, int(math.log2(GENOME_LEN)))
    print(
        f"Recommended star.index_extra for GENOME_LEN={GENOME_LEN}: "
        f'"--genomeSAindexNbases {sa_index_nbases} '
        f'--genomeChrBinNbits {chr_bin_nbits}"'
    )


if __name__ == "__main__":
    main()
