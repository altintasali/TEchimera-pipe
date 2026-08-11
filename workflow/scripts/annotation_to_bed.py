#!/usr/bin/env python3
"""Convert the reference GTF files into the BED tracks the chimera junction
screen overlaps breakpoints against.

Three outputs (all BED, sorted, in --outdir):

  genes.bed   BED6 per gene:      chrom, start, end, gene_id, score, strand
  exons.bed   BED6 per exon:      chrom, start, end, gene_id, score, strand
  te.bed      BED8 per TE locus:  chrom, start, end, te_id, score, strand,
                                  family, class

Gene exons come from the annotation GTF (ref.gtf, the same file STAR's index
uses). TE insertions come from ref.te_gtf -- the curated TE GTF (TEtranscripts
rmsk-based format, `exon` features carrying gene_id/family_id/class_id
attributes), i.e. the same input the user already maintains. Nothing here
depends on the library being stranded or paired-end: the tracks carry strand
annotations, but the breakpoint-overlap test itself is strand-agnostic (strand
is recorded per junction by the parser for downstream use).

Usage:
    annotation_to_bed.py --gtf GENE.GTF --te-gtf TE.GTF --outdir OUTDIR
"""
import argparse
import os
import sys


def parse_attrs(attr_text):
    """GTF attribute field -> {key: value} (values unquoted)."""
    attrs = {}
    for chunk in attr_text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(None, 1)
        if len(parts) == 2:
            key, value = parts
            attrs[key] = value.strip('"')
        elif len(parts) == 1:
            attrs[parts[0]] = ""
    return attrs


def emit(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--te-gtf", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if not os.path.isfile(args.gtf) or not os.path.isfile(args.te_gtf):
        sys.exit(f"error: --gtf / --te-gtf must point at existing files "
                 f"(got {args.gtf!r}, {args.te_gtf!r})")
    os.makedirs(args.outdir, exist_ok=True)

    # gene_id -> {strand, chrom, start, end, family, class}
    genes = {}
    exons = []
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = (
                cols[0], cols[1], cols[2], cols[3], cols[4], cols[5],
                cols[6], cols[7], cols[8],
            )
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            gid = a.get("gene_id")
            if not gid:
                continue
            g = genes.setdefault(
                gid, {"chrom": chrom, "strand": strand, "start": s, "end": e}
            )
            if chrom != g["chrom"]:
                g["chrom"] = chrom  # keep first; multi-contig genes are pathological
            g["start"] = min(g["start"], s)
            g["end"] = max(g["end"], e)
            if feat == "exon":
                exons.append((chrom, s, e, gid, strand))

    # TE GTF: features are the insertion loci themselves (exon feature in the
    # TEtranscripts rmsk-derived GTFs); keep family/class for the report.
    te_loci = {}
    with open(args.te_gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = cols[:9]
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            a = parse_attrs(attrs)
            gid = a.get("gene_id")
            if not gid:
                continue
            g = te_loci.setdefault(
                gid,
                {
                    "chrom": chrom, "strand": strand, "start": s, "end": e,
                    "family": a.get("family_id", "."), "class": a.get("class_id", "."),
                },
            )
            g["start"] = min(g["start"], s)
            g["end"] = max(g["end"], e)

    gene_rows = sorted(
        (g["chrom"], g["start"] - 1, g["end"], gid, ".", g["strand"])
        for gid, g in genes.items()
    )
    exon_rows = sorted(
        (chrom, s - 1, e, gid, ".", strand) for chrom, s, e, gid, strand in exons
    )
    te_rows = sorted(
        (g["chrom"], g["start"] - 1, g["end"], gid, ".", g["strand"],
         g["family"], g["class"])
        for gid, g in te_loci.items()
    )

    emit(os.path.join(args.outdir, "genes.bed"), gene_rows)
    emit(os.path.join(args.outdir, "exons.bed"), exon_rows)
    emit(os.path.join(args.outdir, "te.bed"), te_rows)
    print(
        f"wrote {len(gene_rows)} genes, {len(exon_rows)} exons, "
        f"{len(te_rows)} TE loci to {args.outdir}"
    )


if __name__ == "__main__":
    main()
