#!/usr/bin/env bash
# =============================================================================
# Generates Snakemake's execution report (--report) for a finished run: the
# DAG, every rule's benchmark (runtime / peak RSS) charts, and env details,
# all bundled into one shareable HTML file.
#
# Usage (run from the repo root):
#   ./workflow/scripts/benchmark-report.sh                        # -> report.html
#   ./workflow/scripts/benchmark-report.sh -o out/report.html     # custom output path
#   ./workflow/scripts/benchmark-report.sh --configfile config/test.yaml
#
# All extra arguments are passed straight through to snakemake, so
# --configfile / --cores / etc. all work. `-o/--output` is consumed here.
#
# Note: Snakemake's --report re-evaluates the whole DAG, so it only has
# benchmark data to aggregate once the workflow has actually run (see the
# "Resource usage & reports" README section). It must be run from the repo
# root, like snakemake itself.
# =============================================================================
set -euo pipefail

out="report.html"
args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)
            [[ $# -ge 2 ]] || { echo "usage: $0 [-o/--output FILE] [snakemake args...]" >&2; exit 2; }
            out="$2"
            shift 2
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

exec snakemake --report "${out}" "${args[@]}"
