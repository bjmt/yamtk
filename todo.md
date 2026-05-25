## Todo

- wrappers: yamtkr, pyamtk

- yamtk bkg: add a flag to print gc histograms

- yamme -x, yamref -x
- yamshuf: DNA/RNA shuffler
  + streme paper says it preserves positions of separator characters; I should do this too
- yamme: motif elicitation
  + (future) Seed-width trimming during refinement: after each refinement pass,
    check IC of the two outermost columns and drop them if below threshold; continue
    with the narrowed PPM. Would improve p-value sensitivity (fewer uninformative
    PWM degrees of freedom) in addition to the output-time IC trimming already
    implemented (MIN_IC_BITS).
- yamscan: Make only using -x and -s output the subset sequences instead of just
  info about the ranges? --> actually make it a separate program to manipulate seqs
- get rid of infinite for loops, always use the `#define`'d bounds

- multithread low-mem mode?

- possible performance boosts:
  + 5-way PWM stride → 4-way layout. The PWM stores 5 scores per position (4 bases + ambiguity slot) at pwm[base +
    pos*5]. Strict 4-way layout would reduce per-position memory from 20 B to 16 B and might enable SIMD. But this is a
    deep refactor (touches all 4 scoring tools + parsing), and the inner loop only touches one of the 5 slots per
    position — the wasted bytes only matter at cache-line boundaries. Estimated 5–10% in practice, not the 50%+ Agent A
    claimed.
  + SIMD inner loop. ARM NEON on Apple Silicon, AVX on x86. Real wins are plausible (2–4×) but only after the layout
    change. Hard.

