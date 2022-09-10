#!/bin/bash

awk '/^##/ { print ; next }
  { print | "sort '"${SORT_ARGS}"' -t$'\''\t'\'' -k6,6n -k7,7nr" }' \
    < /dev/stdin

