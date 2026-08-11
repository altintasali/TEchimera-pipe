# -----------------------------------------------------------------------------
# Chimera screen: gene-TE chimeric junction detection and annotation.
#
# The STAR alignment run already emits {sample}_Chimeric.out.junction (see
# align.smk --chimOutType). This stage annotates those junctions against the
# gene/TE tracks and aggregates them into the all-events catalog and the
# counts matrix the sample-QC stage (sample_qc.smk) consumes.
#
# Rules:
#   annotation_to_bed        GTF + TE GTF -> BED tracks (once)
#   parse_chimeric_junctions per-sample junction annotation
#   chimera_counts           merge per-sample tables -> all_events + counts
#   junction_qc              per-sample QC summary (MultiQC custom content)
#   chimera_igv_bed          per-sample IGV track (config-gated)
# -----------------------------------------------------------------------------
import os

WRITE_COUNTS = bool(config["chimera"]["outputs"]["write_counts_matrix"])
WRITE_IGV_BED = bool(config["chimera"]["outputs"]["write_igv_bed"])
CHIMERA_QC = config["chimera"]["qc"]


def chimera_counts_input():
    return [
        f"results/chimera/{s}_junctions.tsv"
        for s in SAMPLES
    ]


def all_chimera_outputs():
    """Chimera artifacts for the `all` target (Snakefile)."""
    files = [
        f"results/chimera/{s}_junctions.tsv"
        for s in SAMPLES
    ]
    files += [
        "results/chimera/all_events.tsv",
        "results/chimera/counts_matrix.tsv",
    ]
    files += [
        f"results/chimera/qc/{s}_junction_qc.tsv"
        for s in SAMPLES
    ]
    if WRITE_IGV_BED:
        files += [
            f"results/chimera/igv/{s}_junctions.bed"
            for s in SAMPLES
        ]
    return files


def all_sample_qc_outputs():
    """Sample-QC artifacts for the `all` target (Snakefile). Only produced
    when a counts matrix is written (it is the QC view's input)."""
    if not WRITE_COUNTS:
        return []
    transform = CHIMERA_QC["pca_transform"]
    return [
        f"results/chimera/qc/{transform}_counts.tsv",
        f"results/chimera/qc/pca_{transform}_mqc.json",
        f"results/chimera/qc/heatmap_{transform}_mqc.json",
    ]


rule annotation_to_bed:
    # Converts the gene GTF + the curated TE GTF into the BED tracks the
    # breakpoint-overlap test runs against (genes.bed, exons.bed, te.bed).
    # Pure-python (annotation_to_bed.py), so it runs in the base environment.
    input:
        gtf=GTF,
        te_gtf=TE_GTF,
    output:
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
    params:
        outdir="results/reference",
    threads: get_resources("annotation_to_bed")["threads"]
    resources:
        mem_mb=get_resources("annotation_to_bed")["mem_mb"],
        runtime=get_resources("annotation_to_bed")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/annotation_to_bed/annotation_to_bed.txt",
    log:
        "results/pipeline_info/logs/chimera/annotation_to_bed.log",
    shell:
        "python {SCRIPTS_DIR}/annotation_to_bed.py "
        "--gtf {input.gtf} --te-gtf {input.te_gtf} "
        "--outdir {params.outdir} > {log} 2>&1"


rule parse_chimeric_junctions:
    # Annotates one sample's STAR chimeric junctions against the gene/TE
    # tracks. See parse_chimeric_junctions.py for the full column spec.
    # Pure-python, so it runs in the base environment.
    input:
        junctions="results/star/{sample}_Chimeric.out.junction",
        genes="results/reference/genes.bed",
        exons="results/reference/exons.bed",
        te="results/reference/te.bed",
        strandedness=strandedness_input,
    output:
        "results/chimera/{sample}_junctions.tsv",
    params:
        tolerance=config["chimera"]["breakpoint_tolerance"],
        canonical_flag=lambda wc: (
            "--require-canonical"
            if config["chimera"]["require_canonical_junction"]
            else ""
        ),
        library=get_strandedness_param,
    threads: get_resources("parse_chimeric_junctions")["threads"]
    resources:
        mem_mb=get_resources("parse_chimeric_junctions")["mem_mb"],
        runtime=get_resources("parse_chimeric_junctions")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/parse_chimeric_junctions/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera/parse/{sample}.log",
    shell:
        "python {SCRIPTS_DIR}/parse_chimeric_junctions.py "
        "--junctions {input.junctions} "
        "--genes {input.genes} --exons {input.exons} --te {input.te} "
        "--sample {wildcards.sample} "
        "--breakpoint-tolerance {params.tolerance} "
        "{params.canonical_flag} "
        "--library-strandedness {params.library} "
        "--out {output} > {log} 2>&1"


rule chimera_counts:
    # Merges every sample's junction table into the all-events catalog and
    # the event x sample counts matrix (chimera_counts.py). The counts matrix
    # feeds the sample-QC PCA/clustering stage.
    input:
        tables=chimera_counts_input(),
    output:
        events="results/chimera/all_events.tsv",
        counts="results/chimera/counts_matrix.tsv",
    params:
        sample_names=lambda wc, input: " ".join(SAMPLES),
    threads: get_resources("chimera_counts")["threads"]
    resources:
        mem_mb=get_resources("chimera_counts")["mem_mb"],
        runtime=get_resources("chimera_counts")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_counts/chimera_counts.txt",
    log:
        "results/pipeline_info/logs/chimera/counts.log",
    shell:
        "python {SCRIPTS_DIR}/chimera_counts.py "
        "--tables {input.tables} "
        "--sample-names {params.sample_names} "
        "--out-events {output.events} "
        "--out-counts {output.counts} > {log} 2>&1"


rule junction_qc:
    # Per-sample junction QC summary (junction_qc.py) for the MultiQC custom
    # content table.
    input:
        "results/chimera/{sample}_junctions.tsv",
    output:
        "results/chimera/qc/{sample}_junction_qc.tsv",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("junction_qc")["threads"]
    resources:
        mem_mb=get_resources("junction_qc")["mem_mb"],
        runtime=get_resources("junction_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/junction_qc/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera/junction_qc/{sample}.log",
    shell:
        "python {SCRIPTS_DIR}/junction_qc.py "
        "--table {input} --sample {params.sample} --out {output} > {log} 2>&1"


rule chimera_igv_bed:
    # Optional IGV track per sample (one BED12-style row per junction), so
    # candidates can be inspected visually. Gated by
    # config["chimera"]["outputs"]["write_igv_bed"].
    input:
        "results/chimera/{sample}_junctions.tsv",
    output:
        "results/chimera/igv/{sample}_junctions.bed",
    params:
        sample=lambda wc: wc.sample,
    threads: get_resources("chimera_igv_bed")["threads"]
    resources:
        mem_mb=get_resources("chimera_igv_bed")["mem_mb"],
        runtime=get_resources("chimera_igv_bed")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/chimera_igv_bed/{sample}.txt",
    log:
        "results/pipeline_info/logs/chimera/igv/{sample}.log",
    script:
        "../scripts/junctions_to_igv_bed.py"
