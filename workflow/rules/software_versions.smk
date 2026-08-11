# -----------------------------------------------------------------------------
# Software versions section for MultiQC.
#
# Writes the pinned tool versions (config["versions"]) in the format
# MultiQC's "Software Versions" section expects (a `*_mqc_versions.yaml` /
# `.yml` file, per MultiQC 1.33's search pattern -- it scans for that
# filename, not `multiqc_software_versions.yml`). MultiQC then shows the
# versions in the Software Versions report section and tags the STAR, FastQC,
# samtools etc. rows in the General Stats table with them, so a report always
# records exactly which releases each step ran with.
# -----------------------------------------------------------------------------
rule software_versions:
    output:
        "results/versions/techimera_mqc_versions.yml",
    threads: get_resources("software_versions")["threads"]
    resources:
        mem_mb=get_resources("software_versions")["mem_mb"],
        runtime=get_resources("software_versions")["runtime"],
    benchmark:
        "results/pipeline_info/benchmarks/software_versions/software_versions.txt",
    log:
        "results/pipeline_info/logs/software_versions.log",
    run:
        import yaml

        versions = {
            "STAR": V["star"],
            "SAMtools": V["samtools"],
            "Trim Galore!": V["trim_galore"],
            "FastQC": V["fastqc"],
            "RSeQC": V["rseqc"],
            "MultiQC": V["multiqc"],
            "DESeq2": V["deseq2"],
        }
        with open("results/versions/techimera_mqc_versions.yml", "w") as fh:
            yaml.safe_dump(
                {"software_versions": versions}, fh, default_flow_style=False
            )
