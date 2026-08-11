"""Turn RSeQC infer_experiment.py output into the resolved strandedness value.

RSeQC infer_experiment.py reports, for single-end data:
    Fraction of reads explained by "++,--":   <forward-type fraction>
    Fraction of reads explained by "+-,-+":   <reverse-type fraction>
and for paired-end data:
    Fraction of reads explained by "1++,1--,2+-,2-+": <forward-type fraction>
    Fraction of reads explained by "1+-,1-+,2++,2--": <reverse-type fraction>

Mapping to the workflow's resolved strandedness vocabulary {no,forward,reverse}:
  forward-type fraction dominates -> "forward" (e.g. QIAseq stranded / "second-strand")
  reverse-type fraction dominates -> "reverse" (e.g. Illumina TruSeq stranded / "first-strand")
  neither dominates                -> "no" (unstranded)
"""

import re
import sys

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input.txt) as fh:
    text = fh.read()

fractions = {}
for line in text.splitlines():
    m = re.search(r'explained by "([^"]+)":\s*([\d.]+)', line)
    if m:
        pattern, value = m.group(1), float(m.group(2))
        fractions[pattern] = value

if not fractions:
    raise ValueError(
        f"Could not parse any 'explained by' fractions from {snakemake.input.txt}. "
        "Check the RSeQC infer_experiment.py output format."
    )

forward_value = 0.0
reverse_value = 0.0
for pattern, value in fractions.items():
    first_token = pattern.split(",")[0]  # "++" / "1++" / "+-" / "1+-"
    if first_token.endswith("++"):
        forward_value = value
    elif first_token.endswith("+-"):
        reverse_value = value

min_fraction = float(snakemake.params.min_fraction)

print(f"Parsed fractions: {fractions}", file=sys.stderr)
print(
    f"forward-type fraction: {forward_value:.4f}, "
    f"reverse-type fraction: {reverse_value:.4f}, "
    f"min_fraction threshold: {min_fraction}",
    file=sys.stderr,
)

if forward_value >= min_fraction:
    stranded = "forward"
elif reverse_value >= min_fraction:
    stranded = "reverse"
else:
    stranded = "no"

print(f"--> calling library strandedness: {stranded}", file=sys.stderr)

with open(snakemake.output.txt, "w") as fh:
    fh.write(stranded + "\n")
