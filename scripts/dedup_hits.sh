#!/bin/bash

# Warning: While rather memory efficient, this is still quite a slow solution
# to the problem. One way to (somewhat forcefully) speed it up would be to
# index the positions of motifs during sorting and dedup them in separate
# threads.

# 1. Sort by motif, seqname, strand, and descending score.
# 2. For each motif/seqname/strand combination, check if the hit overlaps
#    a previous (higher scoring) hit; if so, skip.

awk '/^##/ { print ; next }
  {
    if (NF != 9 && NF != 10) {
      print "Error: Input is malformed; expected 9-10 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }
  { print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k5,5 -k1,1 -k4,4 -k7,7nr" }' \
    < /dev/stdin \
  | awk '
    BEGIN {
      IFS = "\t"
      OFS = "\t"
      chr1 = ""
      strand4 = ""
      motif5 = ""
      n = 0
      skip = 0
    }
    /^##/ {
      if ($0 ~ /^##MotifCount/ && $0 !~ /Dedupped[=]true/) {
        print $0 " Dedupped=true"
      } else {
        print
      }
      next
    }
    {
      if (NR % 100000 == 0) printf("\rFinished %'"'"'d lines ...", NR) > "/dev/stderr"
      if ($1 == chr1 && $4 == strand4 && $5 == motif5) {
        for (i = 1; i <= n; i++) {
          if (($2 >= start2[i] && $2 <= stop3[i]) || ($3 >= start2[i] && $3 <= stop3[i])) {
            skip = 1
            break
          }
        }
        if (skip == 1) {
          skip = 0
          next
        }
        n += 1
        start2[n] = $2
        stop3[n] = $3
        print
      } else {
        chr1 = $1
        strand4 = $4
        motif5 = $5
        n = 1
        start2[n] = $2
        stop3[n] = $3
        print
      }
    }
    END { if (NR > 99999) printf(" done.\n") > "/dev/stderr" }
  '

