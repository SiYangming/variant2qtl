# SuSiE fine-mapping mini testdata

Synthetic single-locus summary statistics for stub / smoke tests.

| File | Role |
| --- | --- |
| `sumstats.tsv` | Columns `variant_id`, `beta`, `se`, `z`, `n` |
| `ld.txt` | Optional 10×10 identity LD (tab-delimited, no header) |

Regenerate LD: `python3 -c "import numpy as np; np.savetxt('ld.txt', np.eye(10), fmt='%.6f', delimiter='\\t')"`
