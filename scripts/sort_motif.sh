#!/bin/bash

awk '/^##/ { print ; next }
  { print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k5,5 -k1,1 -k4,4 -k2,2n" }' \
    < /dev/stdin

