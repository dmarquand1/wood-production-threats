# Code for "Global threats of wood production to vertebrate biodiversity"

Scripts are in three folders, run in numbered order:
- scripts/01_extent: threat summaries, Fig. 1, Figs S6–S7
- scripts/02_hotspots: maps and null models, Figs 2–3, Figs S1–S5
- scripts/03_traits: trait imputation and brms models (one folder per taxon), Fig. 4
  Model fits (fit_*_hpc.R) run one chain per call: Rscript fit_logging_hpc.R <1-4>

Data are not included. IUCN Red List assessments and range maps (v2025-2) are
available from the IUCN, bird ranges from BirdLife International, and trait data
and phylogenies from the sources cited in the paper.