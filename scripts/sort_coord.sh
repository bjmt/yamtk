#!/bin/bash

awk '/^##/ { print ; next }
  { print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k1,1 -k2,2n -k3,3n -k4,4 -k5,5" }' \
    < /dev/stdin

