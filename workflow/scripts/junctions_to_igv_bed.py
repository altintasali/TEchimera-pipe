"""Turn a per-sample chimera junction table (parse_chimeric_junctions.py)
into a BED track for IGV: one BED6-ish row per gene-TE event, spanning the
donor->acceptor breakpoint.

The script runs under Snakemake's `script:` directive, so it reads the
snakemake.input / snakemake.output / snakemake.params globals.

BED columns written (BED6, score = supporting reads, strand = donor strand):
    chrom  donor_bp  acceptor_bp  event_id  reads  donor_strand
"""
import os

with open(snakemake.input[0]) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    rows = [
        dict(zip(header, line.rstrip("\n").split("\t")))
        for line in fh
        if line.strip()
    ]

os.makedirs(os.path.dirname(str(snakemake.output[0])), exist_ok=True)
with open(snakemake.output[0], "w") as fh:
    fh.write('track name="chimera_junctions" description="gene-TE chimera '
             'junctions ({sample})" itemRgb="On"\n')
    for r in rows:
        if r.get("direction") not in ("gene_to_te", "te_to_gene"):
            continue
        try:
            donor, acceptor = int(r["donor_breakpoint"]), int(r["acceptor_breakpoint"])
        except (KeyError, ValueError):
            continue
        # BED is 0-based half-open; STAR breakpoints are 1-based inclusive,
        # so span [donor-1, acceptor).
        start, end = min(donor - 1, acceptor), max(donor - 1, acceptor)
        fh.write("\t".join([
            r.get("chrom", "."), str(start), str(end),
            r.get("event_id", "."), r.get("reads", "1"),
            r.get("donor_strand", "."),
        ]) + "\n")
