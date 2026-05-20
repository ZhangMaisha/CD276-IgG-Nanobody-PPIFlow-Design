# CD276 IgG and Nanobody Design Using PPIFlow

This repository contains the design configurations, selected structures, helper scripts, and summary outputs used to evaluate CD276/B7-H3-targeting IgG antibody and nanobody candidates. The project workflow moves from epitope support to large-scale design generation, Prodigy ranking, AF3 recheck, and final 0-50 ns MD evaluation. It depends on the upstream PPIFlow repository: https://github.com/Mingchenchen/PPIFlow

## Workflow

1. Define the CD276 epitope
2. Generate antibody and nanobody backbones with PPIFlow
3. Run AbMPNN sequence optimisation
4. Run Flowpacker side-chain packing
5. Score and rank candidates with Prodigy
6. Validate selected candidates with AF3
7. Run MD for selected structures and retain the continuous 0-50 ns summary as the primary MD result

## Repository Layout

- `configs/`: example YAML files, generated batch YAML snapshots, and GROMACS MDP files
- `input/`: selected structures and AF3 batch JSON inputs used for downstream validation
- `scripts/design/`: design-generation and candidate-export helpers
- `scripts/prodigy/`: Prodigy scoring, recheck parsing, and plotting helpers
- `scripts/md/`: MD setup, extension, and post-processing helpers
- `scripts/archive/`: archived helper scripts retained for reference only
- `results/epitope_analysis/`: epitope-support outputs
- `results/prodigy/initial_screening/`: full-candidate Prodigy ranking summaries
- `results/prodigy/af3_recheck/`: AF3 recheck tables, source workbook, and derived figures
- `results/md_analysis/primary_0_50ns/`: primary reported MD summary and comparison figures
- `ppiflow_af3_merged.yaml`: merged Conda environment for the PPIFlow and AF3-related workflow

## Main Outputs

- `results/epitope_analysis/discotope/file_A_discotope3.csv`: residue-level DiscoTope output supporting epitope selection
- `results/prodigy/initial_screening/antibody_ranked_all.csv`: full initial Prodigy ranking for IgG candidates
- `results/prodigy/initial_screening/nanobody_ranked_all.csv`: full initial Prodigy ranking for nanobody candidates
- `results/prodigy/af3_recheck/final_AF3_prodigy_summary.csv`: consolidated Prodigy comparison before and after AF3 re-evaluation
- `results/prodigy/af3_recheck/figures/Figure2A_binding_energy.png`: summary figure for the selected AF3-rechecked candidates
- `results/md_analysis/primary_0_50ns/`: primary MD result directory containing both summary tables and comparison figures

## Notes

This repository does not vendor the PPIFlow source code. Please install or clone PPIFlow separately from `https://github.com/Mingchenchen/PPIFlow` when reproducing the workflow.

Large intermediate files, raw model outputs, and MD trajectory/state files are intentionally excluded. This repository keeps selected inputs, summary tables, and derived figures rather than full AF3 download bundles or full MD trajectories.

Several execution scripts still assume the original external workspace layout under `/teams/.../PPIFlow-main`. They document how the runs were carried out, but they are not standalone scripts for this export repository.

`scripts/prodigy/plot_prodigy_recheck.R` expects the manually consolidated Prodigy workbook to be available in the working directory when the script is run. A copy is retained under `results/prodigy/af3_recheck/source_tables/`.

`scripts/archive/10ns_md_analysis.R` is kept only as an archived reference script. The 10 ns-only MD result folders are stored outside this repository because the continuous 0-50 ns MD summary is treated as the primary reported MD result.
