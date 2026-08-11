def _fastqc_base_name(fname):
    """The name FastQC derives from an input fastq for its output .zip:
    the basename with recognized file extensions stripped (a leading .gz
    and then .fastq/.fq/.fasta/.fa/.txt). Used to find the zip FastQC
    produced so it can be renamed to the canonical {sample}_R{read} name."""
    base = os.path.basename(fname)
    for ext in (".gz", ".bz2", ".fastq", ".fq", ".fasta", ".fa", ".txt"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


rule fastqc_raw:
    # Always-on FastQC of the raw (merged) input fastqs, independent of the
    # optional `trimming` step, so the MultiQC report covers the untrimmed
    # input even when trimming is disabled. FastQC names its .zip after the
    # input file's basename, which only matches the canonical
    # {sample}_R{read} name for the merged results/fastq/ files -- a raw
    # single-lane path would produce a differently-named zip, so it is
    # renamed to the declared output here.
    input:
        fq=lambda wc: merged_fastq_path(wc.sample, int(wc.read)),
    output:
        zip="results/fastqc/raw/{sample}_R{read}_fastqc.zip",
    params:
        outdir="results/fastqc/raw",
    threads: get_resources("fastqc_raw")["threads"]
    resources:
        mem_mb=get_resources("fastqc_raw")["mem_mb"],
        runtime=get_resources("fastqc_raw")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/fastqc_raw/{sample}_R{read}.txt",
    log:
        "results/pipeline_info/logs/fastqc/raw/{sample}_R{read}.log",
    conda:
        FASTQC_ENV
    run:
        os.makedirs(params.outdir, exist_ok=True)
        shell(
            "fastqc -t {threads} -o {params.outdir} {input.fq} "
            "> {log} 2>&1"
        )
        produced = os.path.join(
            params.outdir, _fastqc_base_name(input.fq) + "_fastqc.zip"
        )
        if os.path.abspath(produced) != os.path.abspath(output.zip):
            os.replace(produced, output.zip)


rule benchmark_summary:
    # Aggregates every benchmark file this configuration produces into a
    # self-describing MultiQC custom content table ("Resource usage",
    # id resource_usage), rendered after the RSeQC section inside the
    # multiqc_report.html via module_order in multiqc_config.yaml.
    input:
        all_benchmark_files(),
    output:
        "results/pipeline_info/benchmark_summary_mqc.json",
    params:
        allocated=allocated_resources_by_rule(),
    threads: get_resources("benchmark_summary")["threads"]
    resources:
        mem_mb=get_resources("benchmark_summary")["mem_mb"],
        runtime=get_resources("benchmark_summary")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/benchmark_summary/benchmark_summary.txt",
    log:
        "results/pipeline_info/logs/multiqc/benchmark_summary.log",
    conda:
        MULTIQC_ENV
    script:
        "../scripts/benchmark_summary.py"


rule multiqc:
    # Aggregates STAR alignment logs, (if trimming is enabled) TrimGalore!
    # + FastQC reports, the always-on raw FastQC reports, (if strandedness
    # auto-detection was used) RSeQC infer_experiment.py reports, the
    # per-rule benchmark/resource summary, and the pinned tool versions into
    # one HTML report. Runs the MultiQC version pinned in
    # config["versions"]["multiqc"]. The custom config (multiqc_config.yaml)
    # trims "_val_1"/"_val_2"/"_trimmed" off the FastQC sample names so every
    # module's rows merge into one clean row per sample in General Stats.
    input:
        expand("results/star/{sample}_Log.final.out", sample=SAMPLES),
        expand("results/rseqc/{sample}_infer_experiment.txt", sample=AUTO_SAMPLES),
        all_fastqc_reports(),
        all_raw_fastqc_reports(),
        "results/pipeline_info/benchmark_summary_mqc.json",
        "results/versions/techimera_mqc_versions.yml",
    output:
        html="results/qc/multiqc_report.html",
        data=directory("results/qc/multiqc_report_data"),
    params:
        extra="",
        indirs=lambda wc, input: sorted({os.path.dirname(f) for f in input}),
    threads: get_resources("multiqc")["threads"]
    resources:
        mem_mb=get_resources("multiqc")["mem_mb"],
        runtime=get_resources("multiqc")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/multiqc/multiqc.txt",
    log:
        "results/pipeline_info/logs/multiqc/multiqc.log",
    conda:
        MULTIQC_ENV
    shell:
        "multiqc {params.extra} --force "
        "-c workflow/default-config/multiqc_config.yaml "
        "-o results/qc -n multiqc_report "
        "{params.indirs} "
        "> {log} 2>&1\n"
