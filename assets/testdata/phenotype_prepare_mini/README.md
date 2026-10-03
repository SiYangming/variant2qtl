Mini molecular phenotype matrix for `phenotype_prepare` stub/smoke tests.

- `expr.tsv`: genes × samples (`s1`–`s8` plus `s_drop`); `GENE_DROP` is all-missing and should be filtered.
- `genes.bed`: chr/start/end/`gene_id` for FastQTL BED coordinates.
- `samples.txt`: keep list (`s1`–`s8`); `s_drop` is excluded by intersection.
