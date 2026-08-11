#!/usr/bin/env python3
"""Per-sample junction QC summary from a parse_chimeric_junctions.tsv table.

Writes a tiny two-column TSV (metric, value) with the event counts a user
wants at a glance before any deeper QC: total events, how many are gene<->TE
candidates, the direction/chimera-type composition, canonical vs
non-canonical, antisense, and strand-match. Used by the MultiQC custom-content
module (parse_junction_qc) so the report can show per-sample chimera QC.
"""
import argparse
import collections
import os


def load(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = list(load(args.table))
    n = len(rows)
    gene_te = [r for r in rows if r.get("direction") in ("gene_to_te", "te_to_gene")]
    direction = collections.Counter(r.get("direction", ".") for r in rows)
    chimera_type = collections.Counter(r.get("chimera_type", ".") for r in gene_te)
    canonical = collections.Counter(r.get("canonical", "no") for r in rows)
    antisense = collections.Counter(r.get("antisense_flag", ".") for r in gene_te)
    strand_match = collections.Counter(r.get("gene_strand_match", "NA") for r in gene_te)

    lines = [
        ("sample", args.sample),
        ("events_total", str(n)),
        ("events_gene_te", str(len(gene_te))),
    ]
    for key in ("gene_to_te", "te_to_gene", "gene_to_gene", "te_to_te",
                "gene_to_other", "other_to_gene", "te_to_other", "other"):
        lines.append((f"direction_{key}", str(direction.get(key, 0))))
    for key in ("te_initiated", "te_terminated", "te_exonized", "trans", "."):
        lines.append((f"chimera_type_{key}", str(chimera_type.get(key, 0))))
    lines.append(("canonical_yes", str(canonical.get("yes", 0))))
    lines.append(("canonical_no", str(canonical.get("no", 0))))
    lines.append(("antisense_yes", str(antisense.get("yes", 0))))
    for key in ("yes", "no", "NA"):
        lines.append((f"strand_match_{key}", str(strand_match.get(key, 0))))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("metric\tvalue\n")
        for metric, value in lines:
            fh.write(f"{metric}\t{value}\n")
    print(f"{args.sample}: junction QC ({n} events) -> {args.out}")


if __name__ == "__main__":
    main()
