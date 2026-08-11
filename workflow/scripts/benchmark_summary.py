#!/usr/bin/env python3
"""Aggregate Snakemake benchmark files into a self-describing MultiQC custom
content file rendered as the "Resource usage" section of the workflow's
MultiQC report.

Each input is a Snakemake benchmark .txt written by Snakemake >=8
(tab-separated rows with the header
`s  h:m:s  max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time`,
where max_rss is peak resident memory in MB and mean_load is the average CPU
load as a percentage of one core, e.g. 200 == two cores; `cpu_time` is the
total CPU time in seconds. For robustness the parser also tolerates the legacy
7-column format (`s  seconds  threads  cpu_percent  max_rss  bytes_read
bytes_written`, max_rss in KB) where present.

Files are grouped by rule (the benchmark subdirectory name). Per-rule summary
stats -- job count, mean/max wall time, and for CPU and RAM the allocated
amount (from the workflow's resources.yaml, passed in via `params.allocated`),
the mean/max amount actually used, and the mean used / allocated efficiency --
are written as a table under the fixed id "resource_usage" so the module_order
in multiqc_config.yaml places it after the RSeQC section.
"""
import json
import os
import statistics
from collections import defaultdict

try:
    # Only defined when run by Snakemake's `script:` directive.
    snakemake
except NameError:  # pragma: no cover
    snakemake = None


def _as_float(val, default=0.0):
    """Coerce a benchmark cell to float; benchmark files can hold 'NA' for
    memory/CPU when a job was too fast to sample, and missing columns return
    None."""
    if val is None:
        return default
    if isinstance(val, (int, float)):
        return float(val)
    try:
        return float(str(val).strip())
    except ValueError:
        return default


def _parse_benchmark(path):
    """Yield one dict per run row in a Snakemake benchmark file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if header is None:
                header = line.split("\t")
                continue
            yield dict(zip(header, line.split("\t")))


def main():
    allocated = {}
    if snakemake is not None:
        params = getattr(snakemake, "params", None)
        if params is not None:
            allocated = getattr(params, "allocated", {}) or {}

    rows_by_rule = defaultdict(list)
    for path in snakemake.input:
        rule = os.path.basename(os.path.dirname(path))
        rows_by_rule[rule].extend(_parse_benchmark(path))

    data = {}
    # Snakemake's own rule order (params.allocated is built from the
    # workflow.rules OrderedDict); any rule missing there falls back to the
    # input file order.
    order = {name: i for i, name in enumerate(allocated)}
    for rule in sorted(rows_by_rule, key=lambda r: order.get(r, len(order))):
        rows = rows_by_rule[rule]
        walltimes = [_as_float(r.get("s")) for r in rows]
        # snakemake >=8 writes mean_load (%; 100 == one core); the legacy
        # format had cpu_percent instead.
        loads = []
        for r in rows:
            pct = _as_float(r.get("mean_load")) if "mean_load" in r else _as_float(r.get("cpu_percent"))
            loads.append(pct / 100.0)
        # max_rss is in MB in the snakemake >=8 format (KB in the legacy one,
        # indistinguishable without extra config -- assume MB, i.e. divide by
        # 1024 to get GB).
        rss_gb = [_as_float(r.get("max_rss")) / 1024.0 for r in rows]

        alloc = allocated.get(rule, {})
        cpu_alloc_cores = int(alloc.get("threads") or 0)
        ram_alloc_gb = float(alloc.get("mem_mb") or 0) / 1024.0

        row = {
            "n": len(rows),
            "walltime_mean_h": round(statistics.mean(walltimes) / 3600.0, 3),
            "walltime_max_h": round(max(walltimes) / 3600.0, 3),
            "cpu_alloc_cores": cpu_alloc_cores,
            "cpu_used_mean_cores": round(statistics.mean(loads), 3),
            "cpu_used_max_cores": round(max(loads), 3),
            "ram_alloc_gb": round(ram_alloc_gb, 3),
            "ram_used_mean_gb": round(statistics.mean(rss_gb), 3),
            "ram_used_max_gb": round(max(rss_gb), 3),
        }
        if cpu_alloc_cores > 0:
            row["cpu_eff"] = round(row["cpu_used_mean_cores"] / cpu_alloc_cores, 3)
        if ram_alloc_gb > 0:
            row["ram_eff"] = round(row["ram_used_mean_gb"] / ram_alloc_gb, 3)
        data[rule] = row

    summary = {
        "id": "resource_usage",
        "section_name": "Resource usage",
        "description": (
            "Per-rule job count, wall time, and resource efficiency -- for "
            "CPU and RAM: the allocated amount (resources.yaml), the mean/max "
            "amount actually used (Snakemake benchmark files in "
            "results/pipeline_info/benchmarks/), and the mean used/allocated "
            "ratio. Useful for sizing resources on your cluster before a full "
            "run."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "resource_usage_table",
            "title": "Resource usage",
            "col1_header": "Rule",
            "sort_rows": False,
        },
        "headers": {
            "n": {
                "title": "N",
                "description": "Number of jobs (samples or files)",
                "format": "{:,d}",
                "min": 0,
            },
            "walltime_mean_h": {
                "title": "Wall time mean (h)",
                "format": "{:.3f}",
                "min": 0,
            },
            "walltime_max_h": {
                "title": "Wall time max (h)",
                "format": "{:.3f}",
                "min": 0,
            },
            "cpu_alloc_cores": {
                "title": "CPU allocated (cores)",
                "description": "Threads allocated per job (resources.yaml)",
                "format": "{:,d}",
                "min": 0,
            },
            "cpu_used_mean_cores": {
                "title": "CPU used mean (cores)",
                "description": "Mean average CPU cores used per job (mean_load / 100)",
                "format": "{:.3f}",
                "min": 0,
            },
            "cpu_used_max_cores": {
                "title": "CPU used max (cores)",
                "description": "Max average CPU cores used across jobs (mean_load / 100)",
                "format": "{:.3f}",
                "min": 0,
            },
            "cpu_eff": {
                "title": "CPU efficiency",
                "description": "Mean CPU cores used / allocated threads",
                "format": "{:.0%}",
                "min": 0,
                "max": 1,
            },
            "ram_alloc_gb": {
                "title": "RAM allocated (GB)",
                "description": "Memory allocated per job (resources.yaml mem_mb)",
                "format": "{:.3f}",
                "min": 0,
            },
            "ram_used_mean_gb": {
                "title": "RAM used mean (GB)",
                "description": "Mean peak resident memory used per job (max_rss)",
                "format": "{:.3f}",
                "min": 0,
            },
            "ram_used_max_gb": {
                "title": "RAM used max (GB)",
                "description": "Max peak resident memory used across jobs (max_rss)",
                "format": "{:.3f}",
                "min": 0,
            },
            "ram_eff": {
                "title": "RAM efficiency",
                "description": "Mean peak RAM used / allocated memory",
                "format": "{:.0%}",
                "min": 0,
                "max": 1,
            },
        },
        "data": data,
    }

    with open(snakemake.output[0], "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
