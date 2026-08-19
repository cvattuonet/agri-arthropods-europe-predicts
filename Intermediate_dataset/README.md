# Intermediate_dataset

This folder is populated by running scripts 1-5 of the pipeline. Its contents
are not committed to the repository (they're regenerable and, in most cases,
too large for git) except for the three summary tables below, which are kept
so the final site-level data can be explored without re-running the full
pipeline:

- `sites_for_abundance_models.csv` - produced by script 5, used as input to
  scripts 7, 8 and 9 (abundance models).
- `sites_for_richness_models.csv` - produced by script 5, used as input to
  scripts 7, 8 and 9 (richness models).
- `abundance_sites_for_vizualization_nonscale.csv` - produced by script 5,
  unscaled site-level table kept for visualization/data characterisation.

Everything else in this folder (raw extraction checkpoints, per-predictor
`.rds` files, full observation-level tables, etc.) is regenerated automatically
when you run scripts 1 through 5 in order against the downloaded databases
(see `Databases/*/README.md`).
