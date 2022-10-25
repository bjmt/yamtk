#!/bin/bash

awk '/^##/ { print ; next }
  { 
    if (NF == 9 || NF == 10) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k5,5 -k1,1 -k4,4 -k2,2n"
    } else if (NF == 11 || NF == 12) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k7,7 -k3,3 -k6,6 -k4,4n"
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }' < /dev/stdin

