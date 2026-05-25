#!/usr/bin/env bash
# Target sequences much longer than any pool seq → pool mode must warn on
# stderr and skip the target (no record emitted), but still succeed for any
# matchable target. Here every target is too long, so we should also see a
# non-zero exit code because zero outputs were produced.
source "$TESTDIR/lib.sh"

tmpdir=$(mktemp -d -t yamtkbkgnomatch.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

cat >"$tmpdir/targets.fa" <<'EOF'
>tgt1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>tgt2
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
EOF

cat >"$tmpdir/pool.fa" <<'EOF'
>p1
ACGT
>p2
GGCC
EOF

assert_stderr_contains "bkg pool no-match warns on unfillable target" \
    "no pool match" \
    "$YAMTK" bkg -i "$tmpdir/targets.fa" -p "$tmpdir/pool.fa" -s 1 -T 5

assert_exit_nonzero "bkg pool no-match exits non-zero when nothing emitted" \
    "$YAMTK" bkg -i "$tmpdir/targets.fa" -p "$tmpdir/pool.fa" -s 1 -T 5
