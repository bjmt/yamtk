#!/bin/bash

awk '/^##/ { print ; next }
  { 
    if (NF == 9 || NF == 10) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k6,6n -k7,7nr"
    } else if (NF == 11 || NF == 12) {
      print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k8,8n -k9,9nr"
    } else {
      print "Error: Input is malformed; expected 9-12 fields, found " NF > "/dev/stderr"
      exit 1
    }
  }' < /dev/stdin

