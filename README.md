# numt_detection
A pipeline for NUMT detection using a pileup strategy.

## One-command end-to-end run

You can now run discovery + downstream best-hit analysis in one step:

```bash
bash run_numt_end2end.sh --config numt_pipeline.config
```

## Config-driven parameters

All input paths and parameters are centralized in a config file.

1. Copy the example config:

```bash
cp numt_pipeline.config.example numt_pipeline.config
```

2. Edit `numt_pipeline.config` with your sample paths and thresholds.

The same config is consumed by:
- `run_numt_discovery.sh`
- `process_numt_candidates_one_ready.sh`
- `run_numt_end2end.sh`

If CRAM/CRAI may exist in two locations, set both primary and fallback paths in config:
- `INPUT_BAM_CRAM` / `INPUT_INDEX`
- `INPUT_BAM_CRAM_ALT` / `INPUT_INDEX_ALT`

When primary files are missing, the pipeline automatically falls back to the alternate paths.

## Run all samples with one Slurm array submission

Sample list must be a 3-column TSV:
1. sample ID
2. real species name
3. reference species name

Then submit once:

```bash
bash submit_numt_end2end_array.sh --config numt_pipeline.config --concurrent 50
```

This will automatically create a Slurm job array and process all rows in `SAMPLES_FILE`.
For each sample, CRAM is searched in `CRAM_ROOT_1` first, then `CRAM_ROOT_2`.
