# -----------------------------------------------------------------------------
# Chimera sample-QC: PCA / sample-clustering view of the chimera counts
# matrix (results/chimera/counts_matrix.tsv), DESeq2-normalized
# (vst/rlog) or log2, colored by the sample sheet's "condition" column.
#
# Rules:
#   sample_qc_transform  normalize the counts matrix for the QC view
#   sample_qc            PCA + clustering plots from the transformed matrix
#
# QC filters (chimera.qc: min_samples_present / min_total_counts / min_events
# / pca_transform) apply ONLY to this view -- the all-events catalog and the
# counts matrix are never reduced. Both rules run in the dedicated R env
# (CHIMERA_QC_ENV, common.smk): deseq2 + r-pheatmap.
# -----------------------------------------------------------------------------
rule sample_qc_transform:
    input:
        counts="results/chimera/counts_matrix.tsv",
    output:
        "results/chimera/qc/{transform}_counts.tsv",
    params:
        samples=config["samples"],
        min_samples_present=CHIMERA_QC["min_samples_present"],
        min_total_counts=CHIMERA_QC["min_total_counts"],
    threads: get_resources("sample_qc_transform")["threads"]
    resources:
        mem_mb=get_resources("sample_qc_transform")["mem_mb"],
        runtime=get_resources("sample_qc_transform")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/sample_qc_transform/{transform}.txt",
    log:
        "results/pipeline_info/logs/chimera/sample_qc/transform_{transform}.log",
    conda:
        CHIMERA_QC_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/chimera_sample_qc.R "
        "--transform {input.counts} {params.samples} {wildcards.transform} "
        "{params.min_samples_present} {params.min_total_counts} "
        "{output} > {log} 2>&1"


rule sample_qc:
    # PCA + clustering heatmap of the transformed counts, colored by
    # condition (sample sheet's "condition" column; absent -> one "all" group).
    input:
        transformed="results/chimera/qc/{transform}_counts.tsv",
    output:
        pca="results/chimera/qc/pca_{transform}.svg",
        heatmap="results/chimera/qc/heatmap_{transform}.svg",
    params:
        samples=config["samples"],
        min_events=CHIMERA_QC["min_events"],
    threads: get_resources("sample_qc")["threads"]
    resources:
        mem_mb=get_resources("sample_qc")["mem_mb"],
        runtime=get_resources("sample_qc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/sample_qc/{transform}.txt",
    log:
        "results/pipeline_info/logs/chimera/sample_qc/plots_{transform}.log",
    conda:
        CHIMERA_QC_ENV
    shell:
        "Rscript {SCRIPTS_DIR}/chimera_sample_qc.R "
        "--plots {input.transformed} {params.samples} {params.min_events} "
        "{output.pca} {output.heatmap} > {log} 2>&1"
