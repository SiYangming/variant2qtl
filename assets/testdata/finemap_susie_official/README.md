# SuSiE fine-mapping — official-style testdata

Derived from the **susieR** vignette / package recipe (not a hand-made toy locus):

- Prefer regenerate with R (uses `susieR::N3finemapping` when available):

  ```bash
  conda create -y -n susier -c conda-forge r-susier=0.14.2
  conda run -n susier Rscript scripts/make_finemap_susie_official_testdata.R
  ```

- Python fallback (same vignette simulation recipe; committed by default):

  ```bash
  python3 scripts/make_finemap_susie_official_testdata.py
  ```

| File | Role |
| --- | --- |
| `sumstats.tsv` | `variant_id`, `beta`, `se`, `z`, `n` from univariate OLS |
| `ld.txt` | In-sample correlation matrix (no header) |
| `SOURCE.txt` | Provenance stamp from the generator |

Used by tagged real nf-test: `--tag susie_real` (not run in default CI).
