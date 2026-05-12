# NUMT Detection Pipeline

This repository contains a configurable NUMT (nuclear mitochondrial DNA segment) detection workflow based on mt-read extraction, remapping, and best-hit filtering.

## Repository Layout

- `scripts/`: executable pipeline scripts.
  - `run_numt_discovery.sh`: runs single-sample NUMT discovery.
  - `discover_numt_sinks.py`: detects candidate loci from remapped BAM pileup.
  - `process_numt_candidates_one_ready.sh`: best-hit filtering for one sample BED.
  - `run_numt_end2end.sh`: one-command discovery + best-hit for one sample config.
  - `submit_numt_end2end_array.sh`: Slurm array launcher for multi-sample end-to-end runs.
  - `load_numt_modules.sh`: centralized module and conda environment loading.
  - `run_numt_one_from_3col.sh`, `submit_numt_array.sh`, `submit_numt_besthit_array_ready.sh`, `build_numt_besthit_table.sh`: legacy/auxiliary scripts.
- `config/`: configuration templates.
  - `numt_pipeline.config.example`: example config for array and single-sample modes.
- `data/`: input examples.
  - `primate_numt_test_list.txt`: example 3-column sample table.

## Quick Start

1. Copy the config template:

```bash
cp config/numt_pipeline.config.example numt_pipeline.config
```

2. Edit `numt_pipeline.config` with your paths and thresholds.

3. Run one sample end-to-end:

```bash
bash scripts/run_numt_end2end.sh --config numt_pipeline.config
```

## Pipeline Steps and Script Mapping

### Step 1: Discovery (single sample)
Use `scripts/run_numt_discovery.sh`.

Main operations inside this script:
1. Extract mt-overlapping reads and mates from BAM/CRAM (`samtools view -P`).
2. Convert BAM to FASTQ (`samtools fastq`).
3. Remap reads to nuclear-only reference (`bwa mem` + `samtools sort/index`).
4. Detect NUMT candidate loci by pileup analysis (`scripts/discover_numt_sinks.py`).
5. Write outputs:
   - `SAMPLE.numt_candidates.bed`
   - `SAMPLE.numt_candidates.tsv`

### Step 2: Best-hit filtering (single sample)
Use `scripts/process_numt_candidates_one_ready.sh`.

Main operations:
1. Merge nearby candidate regions.
2. Extract merged nuclear loci sequences.
3. Build chrM BLAST DB.
4. Run BLAST and keep best-hit loci.
5. Write outputs including:
   - `SAMPLE.numt_vs_chrM.besthit.tsv`
   - `SAMPLE.highconf_numt.bed`

### Step 3: End-to-end wrapper (single sample)
Use `scripts/run_numt_end2end.sh`.

This wrapper runs in order:
1. `scripts/run_numt_discovery.sh`
2. `scripts/process_numt_candidates_one_ready.sh`

### Step 4: Slurm array (multi-sample)
Use `scripts/submit_numt_end2end_array.sh`.

```bash
bash scripts/submit_numt_end2end_array.sh --config numt_pipeline.config --concurrent 50
```

Requirements for array mode (config keys):
- `SAMPLES_TSV`
- `CRAM_ROOT_1`, `CRAM_ROOT_2`
- `WHOLE_REF_DIR`, `NUCLEAR_ONLY_REF_DIR`, `CHRM_REF_DIR`
- `DISCOVERY_OUTROOT`, `BESTHIT_OUTDIR`

The script resolves per-sample CRAM/CRAI, generates a temporary per-sample config, and calls `scripts/run_numt_end2end.sh`.

## Environment Setup

Environment loading is centralized in `scripts/load_numt_modules.sh`.

Default modules/env:
- `SAMtools/1.21-GCC-13.3.0`
- `BWA/0.7.18-GCCcore-13.3.0`
- `BEDTools/2.31.1-GCC-13.3.0`
- `miniconda/24.11.3`
- conda env `blast_env`

You can override before running:

```bash
export NUMT_MODULE_SAMTOOLS="SAMtools/1.22-GCC-14.2.0"
export NUMT_CONDA_ENV="blast_env"
```
