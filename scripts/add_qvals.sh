#!/bin/bash

# 1. Sort by P-value and descending score.
# 2. Check for and retrieve the MaxPossibleHits field in the header.
# 3. Calculate the adjusted P-values and print.
# 4. Sort by P-value in inverse order.
# 5. Adjust Q-values so that they are monotonically decreasing.

awk '/^##/ { print ; next }
  {
    if (NF != 9 && NF != 10) {
      print "Error: Input is malformed; expected 9-10 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }
  { print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k6,6n -k7,7nr" }' \
    < /dev/stdin \
  | awk '

    BEGIN {
      IFS = "\t"
      OFS = "\t"
      maxN = -1
      nHeader = 0
    }

    /^##/ {
      nHeader += 1
      if ($0 ~ /MaxPossibleHits/) {
        if (maxN != -1) {
          print "Error: Found multiple MaxPossibleHits fields in header" > "/dev/stderr"
          exit 1
        }
        if ($0 ~ /Dedupped[=]true/) {
          print "Error: Q-values will be incorrect if hits were dedupped" > "/dev/stderr"
          exit 1
        }
        for (i = 1; i <= NF; i++) {
          if ($i ~ /MaxPossibleHits/) {
            split($i, tmpVar, /[=]/)
            maxN = int(tmpVar[2])
          }
        }
        if (maxN < 1) {
          print "Error: Missing valid MaxPossibleHits header info" > "/dev/stderr"
          exit 1
        }
        print
      } else if ($0 ~ /##seqname/ && NF == 9) {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,"qvalue"
      } else {
        print
      }
      next
    }

    {
      qval = ($6 * maxN) / (NR - nHeader)
      if (qval > 1) {
        qval = 1
      }
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,qval | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k6,6nr -k7,7n"
    }

  ' \
  | awk '
    BEGIN {
      IFS = "\t"
      OFS = "\t"
      minQval = 1
    }
    /^##/ {
      print
      next
    }
    {
      qval = $10
      if (qval > minQval) {
        qval = minQval
      } else {
        minQval = qval
      }
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,qval
    }
  '

