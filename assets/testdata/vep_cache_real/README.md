# VEP real-smoke cache

Not committed (large). Fetch with:

```bash
bash scripts/fetch_vep_cache_testdata.sh
```

Then run:

```bash
nf-test test modules/nf-core/ensemblvep/vep/tests/main.vep_real.nf.test --tag vep_real --profile docker
```
