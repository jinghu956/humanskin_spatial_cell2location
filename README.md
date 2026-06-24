# humanskin_spatial_cell2location

This repository contains scripts, logs, summary tables, model outputs, and figures for a human skin spatial transcriptomics cell2location analysis.

## Contents

- `cell2location_run/scripts/`: analysis scripts for spatial RNA preparation, reference processing, cell2location model fitting, abundance summaries, QC, and plotting.
- `cell2location_run/logs/`: run logs and analysis notes, including a Chinese summary log.
- `cell2location_run/results/`: lightweight CSV result tables and generated figures.
- `cell2location_run/models/`: trained model directories that are small enough to version in git.
- `protein_spatial_cluster.R`: related spatial/protein analysis script.

## Data Availability

Original input data and large reference files are intentionally excluded from this repository:

- `input/`
- `reference_atlas/`
- large derived `.h5ad` files

These files are required to fully rerun the workflow from raw data, but are not uploaded to avoid storing raw data and large binary files in GitHub.

## Main Outputs

Key derived result tables include:

- `cell2location_run/results/cell_abundance_broad_adult.csv`
- `cell2location_run/results/cell_abundance_broad_adult_with_coords.csv`
- `cell2location_run/results/cell_abundance_broad_adult_relative_with_coords.csv`
- `cell2location_run/results/spot_dominant_celltype.csv`
- `cell2location_run/results/marker_abundance_correlation_summary.csv`
- `cell2location_run/results/marker_abundance_correlation_heatmap_values.csv`

Main figures are under:

- `cell2location_run/results/figures/`

## Environment

The analysis was run in the `cell2loc_env` conda environment. See `cell2location_run/logs/00_environment_check.log` for package versions and GPU information.
