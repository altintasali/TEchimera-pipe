rule cat_fastq:
    # Concatenates the raw fastq files of a sample sequenced across multiple
    # lanes/runs (multiple rows for the same sample in the sample sheet, the
    # nf-core/rnaseq convention) into one file per read, before trimming/
    # alignment. Only ever requested for samples with >1 file per read (see
    # merged_fastq_path in common.smk: a single-lane sample reads its raw
    # fastq directly), so no unnecessary copy is made. For gzip inputs plain
    # `cat` is used (a concatenation of .gz members is itself a valid gzip
    # file); plain-text lanes are piped through gzip so the merged output is
    # always .gz and STAR's --readFilesCommand zcat applies uniformly.
    input:
        lambda wc: sample_fastqs(wc.sample, wc.read),
    output:
        MERGED_FASTQ_OUTPUT
    params:
        gz=lambda wc: all(
            str(f).endswith(".gz") for f in sample_fastqs(wc.sample, wc.read)
        ),
    threads: get_resources("cat_fastq")["threads"]
    resources:
        mem_mb=get_resources("cat_fastq")["mem_mb"],
        runtime=get_resources("cat_fastq")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/cat_fastq/{sample}_R{read}.txt",
    log:
        "results/pipeline_info/logs/cat_fastq/{sample}_R{read}.log",
    shell:
        "(if [ {params.gz} = True ]; then cat {input} > {output}; "
        "else cat {input} | gzip -c > {output}; fi) 2> {log}"
