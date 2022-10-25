#!/bin/bash

awk '/^##/ { print ; next }
  { 
    if (NF == 9 || NF == 10) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k1,1 -k2,2n -k3,3n -k4,4 -k5,5"
    } else if (NF == 11 || NF == 12) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k3,3 -k4,4n -k5,5n -k6,6 -k7,7"
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }' < /dev/stdin

