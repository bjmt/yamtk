#!/usr/bin/env bash
# Bug #13: shuffle_euler decremented an edge for the terminal vertex,
# potentially corrupting the A-out-edge count of that vertex.
# Test: k=2 shuffle preserves exact overlapping dinucleotide composition.
# Under make debug (ASan) any use of uninit euler_path[i_final] would show.
source "$TESTDIR/lib.sh"

long_fa=$(mktemp /tmp/yamtk_euler_XXXXXX)
trap 'rm -f "$long_fa"' EXIT

# 200-bp sequence to stress the Eulerian path code
python3 -c "
import random; random.seed(0)
seq=''.join(random.choice('ACGT') for _ in range(200))
print('>seq200')
print(seq)
" > "$long_fa"

out=$("$YAMTK" shuf -i "$long_fa" -k 2 -s 99 2>/dev/null)

# Count overlapping dinucleotides; concatenate multi-line FASTA first
count_di() {
    awk '/^>/{next} {seq=seq $0} END{for(i=1;i<length(seq);i++){d=substr(seq,i,2); c[d]++} for(d in c) print c[d], d}' "$@" | sort -k2
}
out_fa=$(mktemp /tmp/yamtk_euler_XXXXXX)
trap 'rm -f "$long_fa" "$out_fa"' EXIT
echo "$out" > "$out_fa"

in_di=$(count_di "$long_fa")
out_di=$(count_di "$out_fa")

if [ "$in_di" != "$out_di" ]; then
    echo "not ok - Eulerian k=2 shuffle on 200bp changed overlapping dinucleotide composition" >&2
    diff <(echo "$in_di") <(echo "$out_di") >&2; exit 1
fi
echo "ok - Eulerian k=2 shuffle on 200bp preserves dinucleotide composition (terminal vertex fix)"
