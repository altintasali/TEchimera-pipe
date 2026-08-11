#!/usr/bin/env python3
"""Parse STAR's Chimeric.out.junction into a per-sample, per-event table of
gene-TE chimera candidates.

STAR's chimeric detection (--chimOutType Junctions WithinBAM SoftClip, see
rules/align.smk) reports, per chimeric read, one line whose 5' segment is the
"donor" and 3' segment the "acceptor". Each segment's breakpoint is tested
against the gene-exon and TE-insertion tracks (annotation_to_bed.py); a read
whose donor hits a gene exon and acceptor a TE (or vice versa) is a
gene<->TE chimera candidate.

This step is deliberately ANNOTATE-ONLY: every junction is written with its
annotations and flags, and nothing is filtered away here (except, optionally,
events whose STAR junction type is non-canonical when require_canonical_
junction is enabled -- see the config schema; the default keeps everything).
Filtering / evidence decisions are left to the user downstream -- the pipeline
ships the full event table and the counts matrix instead.

Output columns (results/chimera/{sample}.junctions.tsv):
    event_id, sample, chrom, donor_breakpoint, donor_strand,
    acceptor_breakpoint, acceptor_strand, junction_type, canonical,
    repeat_flag, reads, donor_hits, acceptor_hits, direction, gene_id,
    gene_strand, te_id, te_family, te_class, chimera_type, antisense_flag,
    library_strand, transcript_strand, gene_strand_match

junction_type/canonical: STAR's column-6 value (0 non-canonical .. 6) and a
derived GT/AG-ish yes/no. TE-involved splicing is often non-canonical, so
this is reported but never filtered (see require_canonical_junction config
for the opt-in).

direction: gene_to_te / te_to_gene (the two primary classes), or one of
gene_to_other / te_to_other / te_to_te / gene_to_gene / other.

chimera_type (only for gene<->TE events): te_initiated (TE upstream of the
gene's TSS on the gene strand), te_terminated (TE downstream of the gene),
te_exonized (TE within the gene body), trans (different chromosomes).

antisense_flag: "yes" when an annotated gene overlaps the TE insertion on the
opposite strand of the assigned gene -- the embedded-TE / sense-antisense
ambiguity class (run twice + IGV curation is the standard treatment for these).

Strand evidence (library_strand / transcript_strand / gene_strand_match):
for stranded libraries the read's aligned strand plus the library type
(forward = read same strand as transcript; reverse = read opposite strand)
gives the transcription strand, which is compared to the annotated gene
strand. For unstranded libraries these columns are NA.
"""
import argparse
import bisect
import os


def load_bed(path, n_extra=0):
    """Return {chrom: sorted list of (start, end, extras_tuple)}."""
    tracks = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            extras = tuple(cols[3 : 6 + n_extra])
            tracks.setdefault(chrom, []).append((start, end, extras))
    for lst in tracks.values():
        lst.sort()
    return tracks


def overlapping(track, chrom, start0, end0):
    """Features in `track` whose [start, end) overlaps [start0, end0).
    `track` must be sorted by start. Returns list of (start, end, extras)."""
    if chrom not in track:
        return []
    starts = [x[0] for x in track[chrom]]
    lo = bisect.bisect_right(starts, start0) - 1
    hits = []
    for i in range(max(lo, 0), len(track[chrom])):
        s, e, extras = track[chrom][i]
        if s >= end0:
            break
        if e > start0:
            hits.append((s, e, extras))
    return hits


CANONICAL_TYPES = {1, 2, 3, 4, 5, 6}  # anything but 0 (non-canonical)


def opp(strand):
    return {"+": "-", "-": "+", ".": "."}.get(strand, ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True, help="STAR {sample}_Chimeric.out.junction")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--te", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--breakpoint-tolerance", type=int, default=0)
    ap.add_argument(
        "--require-canonical", action="store_true",
        help="Only keep junctions whose STAR junction type is canonical "
        "(GT/AG-ish, type != 0). Default: keep everything and record the type.",
    )
    ap.add_argument(
        "--library-strandedness",
        choices=["no", "forward", "reverse", "auto"],
        default="no",
        help="Per-sample strandedness value as resolved by the workflow "
        "(no = unstranded, forward = read on transcript strand, "
        "reverse = read opposite transcript strand). 'auto' is treated as "
        "unstranded here.",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genes = load_bed(args.genes)
    exons = load_bed(args.exons)
    te = load_bed(args.te, n_extra=2)

    tol = max(args.breakpoint_tolerance, 0)
    lib = args.library_strandedness
    if lib == "auto":
        lib = "no"

    gene_meta = {}
    for chrom, feats in genes.items():
        for s, e, ex in feats:
            gene_meta.setdefault(ex[0], (chrom, s, e, ex[2]))
    te_meta = {}
    for chrom, feats in te.items():
        for s, e, ex in feats:
            te_meta.setdefault(ex[0], (chrom, s, e, ex[2], ex[3], ex[4]))

    def locus(bp):
        # STAR breakpoint columns are 1-based and inclusive: the donor's last
        # base is at donor_bp, the acceptor's first base at acceptor_bp.
        # Half-open interval [bp-1-tol, bp+1+tol) covers that base.
        return (bp - 1 - tol, bp + 1 + tol)

    events = {}

    with open(args.junctions) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0]
            try:
                donor_bp, acceptor_bp = int(cols[1]), int(cols[3])
            except ValueError:
                continue
            donor_strand, acceptor_strand = cols[2], cols[4]
            jtype = cols[5]
            repeat_flag = cols[6]

            d0, d1 = locus(donor_bp)
            a0, a1 = locus(acceptor_bp)

            donor_exons = overlapping(exons, chrom, d0, d1)
            donor_tes = overlapping(te, chrom, d0, d1)
            acceptor_exons = overlapping(exons, chrom, a0, a1)
            acceptor_tes = overlapping(te, chrom, a0, a1)

            donor_genes = sorted({ex[2][0] for ex in donor_exons})
            acceptor_genes = sorted({ex[2][0] for ex in acceptor_exons})
            donor_te_ids = sorted({t[2][0] for t in donor_tes})
            acceptor_te_ids = sorted({t[2][0] for t in acceptor_tes})

            donor_gene_hit = bool(donor_genes)
            donor_te_hit = bool(donor_te_ids)
            acceptor_gene_hit = bool(acceptor_genes)
            acceptor_te_hit = bool(acceptor_te_ids)

            if donor_gene_hit and acceptor_te_hit:
                direction = "gene_to_te"
                gene_id = donor_genes[0]
                te_id = acceptor_te_ids[0]
                gene_side_strand = donor_strand
            elif donor_te_hit and acceptor_gene_hit:
                direction = "te_to_gene"
                gene_id = acceptor_genes[0]
                te_id = donor_te_ids[0]
                gene_side_strand = acceptor_strand
            elif donor_gene_hit and acceptor_gene_hit:
                direction = "gene_to_gene"
                gene_id = donor_genes[0]
                te_id = None
                gene_side_strand = donor_strand
            elif donor_te_hit and acceptor_te_hit:
                direction = "te_to_te"
                gene_id = None
                te_id = donor_te_ids[0]
                gene_side_strand = donor_strand
            elif donor_gene_hit:
                direction = "gene_to_other"
                gene_id = donor_genes[0]
                te_id = None
                gene_side_strand = donor_strand
            elif acceptor_gene_hit:
                direction = "other_to_gene"
                gene_id = acceptor_genes[0]
                te_id = None
                gene_side_strand = acceptor_strand
            elif donor_te_hit:
                direction = "te_to_other"
                gene_id = None
                te_id = donor_te_ids[0]
                gene_side_strand = donor_strand
            else:
                direction = "other"
                gene_id = None
                te_id = None
                gene_side_strand = "."

            key = (chrom, donor_bp, donor_strand, acceptor_bp, acceptor_strand,
                   direction)
            ev = events.setdefault(key, {"reads": 0, "gene_id": gene_id,
                                          "te_id": te_id})
            ev["reads"] += 1
            if ev["reads"] == 1:
                ev.update(
                    {
                        "donor_hits": ",".join(donor_genes) or ".",
                        "acceptor_hits": ",".join(acceptor_genes) or ".",
                        "donor_te_hits": ",".join(donor_te_ids) or ".",
                        "acceptor_te_hits": ",".join(acceptor_te_ids) or ".",
                        "junction_type": jtype,
                        "repeat_flag": repeat_flag,
                        "gene_side_strand": gene_side_strand,
                    }
                )
            else:
                # first-seen gene/TE win for stability; keep the union of hits
                ev["donor_hits"] = _union(ev["donor_hits"], donor_genes)
                ev["acceptor_hits"] = _union(ev["acceptor_hits"], acceptor_genes)
                ev["donor_te_hits"] = _union(ev["donor_te_hits"], donor_te_ids)
                ev["acceptor_te_hits"] = _union(ev["acceptor_te_hits"], acceptor_te_ids)

    # Resolve per-event annotations (gene span/strand, TE family/class, type,
    # antisense, strand evidence).
    rows = []
    for (chrom, donor_bp, donor_strand, acceptor_bp, acceptor_strand,
         direction), ev in sorted(events.items()):
        if args.require_canonical:
            try:
                if int(ev["junction_type"]) not in CANONICAL_TYPES:
                    continue
            except ValueError:
                continue
        gene_id = ev["gene_id"]
        te_id = ev["te_id"]

        gene_strand = "."
        gene_span = None
        if gene_id is not None and gene_id in gene_meta:
            _, gs, ge, gene_strand = gene_meta[gene_id]
            gene_span = (gs, ge)

        te_family = te_class = "."
        te_span = None
        if te_id is not None and te_id in te_meta:
            _, ts, tee, _, te_family, te_class = te_meta[te_id]
            te_span = (ts, tee)

        chimera_type = "."
        antisense = "."
        if direction in ("gene_to_te", "te_to_gene") and gene_span and te_span:
            gs, ge, gst = gene_span[0], gene_span[1], gene_strand
            ts, te = te_span[0], te_span[1]
            if gst == "+":
                if te < gs:
                    chimera_type = "te_initiated"
                elif ts > ge:
                    chimera_type = "te_terminated"
                else:
                    chimera_type = "te_exonized"
            elif gst == "-":
                if ts > ge:
                    chimera_type = "te_initiated"
                elif te < gs:
                    chimera_type = "te_terminated"
                else:
                    chimera_type = "te_exonized"
            # antisense: annotated gene overlapping the TE insertion on the
            # strand opposite the assigned gene
            for s, e, ex in overlapping(genes, chrom, te_span[0], te_span[1]):
                if ex[0] != gene_id and ex[2] == opp(gene_strand):
                    antisense = "yes"
                    break

        # strand evidence
        transcript_strand = "NA"
        match = "NA"
        if lib == "forward":
            transcript_strand = ev["gene_side_strand"]
        elif lib == "reverse":
            transcript_strand = opp(ev["gene_side_strand"])
        if transcript_strand != "NA" and gene_strand != ".":
            match = "yes" if transcript_strand == gene_strand else "no"

        try:
            canonical = "yes" if int(ev["junction_type"]) in CANONICAL_TYPES else "no"
        except ValueError:
            canonical = "no"
        event_id = f"{chrom}:{donor_bp}:{donor_strand}:{acceptor_bp}:{acceptor_strand}:{direction}"

        rows.append(
            [
                event_id, args.sample, chrom, donor_bp, donor_strand,
                acceptor_bp, acceptor_strand, ev["junction_type"], canonical,
                ev["repeat_flag"], ev["reads"],
                f"gene:{ev['donor_hits']}|te:{ev['donor_te_hits']}",
                f"gene:{ev['acceptor_hits']}|te:{ev['acceptor_te_hits']}",
                direction,
                gene_id if gene_id is not None else ".",
                gene_strand,
                te_id if te_id is not None else ".",
                te_family, te_class, chimera_type, antisense,
                lib, transcript_strand, match,
            ]
        )

    header = [
        "event_id", "sample", "chrom", "donor_breakpoint", "donor_strand",
        "acceptor_breakpoint", "acceptor_strand", "junction_type", "canonical",
        "repeat_flag", "reads", "donor_hits", "acceptor_hits", "direction",
        "gene_id", "gene_strand", "te_id", "te_family", "te_class",
        "chimera_type", "antisense_flag", "library_strand", "transcript_strand",
        "gene_strand_match",
    ]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"{args.sample}: {len(rows)} events from {args.junctions}")


def _union(csv, ids):
    have = set(x for x in csv.split(",") if x and x != ".")
    have.update(ids)
    return ",".join(sorted(have))


if __name__ == "__main__":
    main()
