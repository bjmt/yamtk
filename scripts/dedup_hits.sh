#!/bin/bash

# For each motif/chr/strand, save hits as they come and for each new
# hit go back through saved results and see if it overlaps anything.

awk '/^##/ { print ; next }
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
      if ($1 == chr1 && $4 == strand4 && $5 == motif5) {
        for (i = 1; i <= n; i++) {
          if (($2 >= start[i] && $2 <= stop[i]) || ($3 >= start[i] && $3 <= stop[i])) {
            skip = 1
            break
          }
        }
        if (skip == 1) {
          skip = 0
          next
        }
        n += 1
        start[n] = $2
        stop[n] = $3
        print
      } else {
        chr1 = $1
        strand4 = $4
        motif5 = $5
        n = 1
        start[n] = $2
        stop[n] = $3
        print
      }
    }
  '

