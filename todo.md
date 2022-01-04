## Todo

- No overlaps (as a separate program)

- multithreading

- Qvalues (as a separate program)
  + BH pvalues: `min(pvals / ((rank(pvals) / nPossibleHits) * 100), 1)`
  + Bonferroni pvalues: `min(pvals * nPossibleHits, 1)`
  + Do a first pass of file to create a Qval table, then read through
    a second time to re-generate results with added Q-values

```
for (int i = PvalLen - 1; i > -1; i--) {
  Qval[i] = max(min(Pval[i] * nPossibleHits / rank, Qval[i+1]), 0.0)
  // rank += PvalCount[i]; // --> Maybe do this before hand? Just a simple cumsum.
                           // Though maybe this is necessary to keep tables separate
                           // per motif.
}
```

Should there be one Pval table per motif? Or a single one?

Creating a Pval table:

  `PvalCount[round(abs(log10(pvalue)) * 1000)]++;`

Table size: ~300 * 1000 = 300,000

            300,000 x 8 bytes = 2,400,000 => 2.4 megabytes

            
           ~300 * 100 = 30,000

            
           ~300 * 10 = 3,000
           3,000 x 8 bytes = 24,000 => 0.23 megabytes
           3,000 x 8 bytes = 24,000 => 0.23 megabytes
