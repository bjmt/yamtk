#!/bin/bash

# Output format: BED6+4
# Score col: -10*log10(P-value)
# Custom cols:
#   7. Log-odds score
#   8. Score as a percent of the max
#   9. P-value
#   10. Q-value (if absent, then simply ".")

awk '
  BEGIN {
    IFS = "\t"
    OFS = "\t"
  }

  /^[^##]/ {
    score = int(-10 * (log($6) / log(10)))
    score = score > 1000 ? 1000 : score
    qvalue = NF > 9 ? $10 : "."
    print $1,$2-1,$3,$5,score,$4,$7,$8,$6,qvalue
  }
  ' < /dev/stdin

