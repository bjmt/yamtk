#!/bin/bash

# Output format: BED6+4
# Score col: -10*log10(P-value)
# Custom cols:
#   7. Log-odds score
#   8. Score percent
#   9. Q-value (if absent, then simply ".")
#   10. Sequence match

awk '
  BEGIN {
    IFS = "\t"
    OFS = "\t"
  }

  /^[^##]/ {
    score = log($6) / log(10)
    score = int(-10 * score)
    score = score > 1000 ? 1000 : score
    qvalue = "."
    if (NF > 9) {
      qvalue = $10
    }
    print $1,$2-1,$3,$5,score,$4,$7,$8,qvalue,$9
  }
  ' < /dev/stdin

