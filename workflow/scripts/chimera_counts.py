#!/usr/bin/env python3
"""Merge the per-sample junction tables (parse_chimeric_junctions.py) into the
chimera all-events catalog and the event x sample counts matrix.

Outputs:

  all_events.tsv   one row per unique event_id across all samples, with the
                   annotation columns taken from the first sample that saw it
                   (annotations are breakpoint-deterministic, so they agree
                   across samples) plus n_samples / total_reads.
  counts_matrix.tsv  event_id x sample read-count matrix (0 where a sample has
                   no reads supporting the event). Written only when the
                   chimera.outputs.write_counts_matrix config is true.

Nothing is filtered here: the full union of annotated events is shipped (the
QC filters in the chimera.qc config section apply only to the PCA/clustering
view in sample_qc.smk).
"""
import argparse
import os
import sys

ANNOTATION_COLUMNS = [
    "event_id", "chrom", "donor_breakpoint", "donor_strand",
    "acceptor_breakpoint", "acceptor_strand", "junction_type", "canonical",
    "repeat_flag", "direction", "gene_id", "gene_strand", "te_id",
    "te_family", "te_class", "chimera_type", "antisense_flag",
    "library_strand", "transcript_strand", "gene_strand_match",
]


def load(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True, nargs="+")
    ap.add_argument("--sample-names", required=True, nargs="+")
    ap.add_argument("--out-events", required=True)
    ap.add_argument("--out-counts", required=False)
    args = ap.parse_args()

    if len(args.tables) != len(args.sample_names):
        sys.exit("error: --tables and --sample-names must have equal length")

    events = {}
    sample_index = {s: i for i, s in enumerate(args.sample_names)}
    for path, sample in zip(args.tables, args.sample_names):
        for row in load(path):
            eid = row["event_id"]
            ev = events.setdefault(eid, {"sample": sample, "counts": {}})
            ev["counts"][sample] = int(row["reads"])
            if ev["sample"] == sample:
                # first-seen annotations (event_id is breakpoint-deterministic)
                for col in ANNOTATION_COLUMNS:
                    ev[col] = row[col]
    # keep the first sample that saw each event (stable order)
    order = sorted(events, key=lambda e: (events[e]["sample"], e))

    with open(args.out_events, "w") as fh:
        header = ANNOTATION_COLUMNS + ["n_samples", "total_reads"]
        fh.write("\t".join(header) + "\n")
        for eid in order:
            ev = events[eid]
            row = [ev.get(c, ".") for c in ANNOTATION_COLUMNS]
            row += [len(ev["counts"]), sum(ev["counts"].values())]
            fh.write("\t".join(str(x) for x in row) + "\n")

    if args.out_counts:
        with open(args.out_counts, "w") as fh:
            fh.write("event_id\t" + "\t".join(args.sample_names) + "\n")
            for eid in order:
                ev = events[eid]
                counts = [ev["counts"].get(s, 0) for s in args.sample_names]
                fh.write(eid + "\t" + "\t".join(str(c) for c in counts) + "\n")

    print(f"{len(order)} unique events across {len(args.sample_names)} samples "
          f"-> {args.out_events}")


if __name__ == "__main__":
    main()
