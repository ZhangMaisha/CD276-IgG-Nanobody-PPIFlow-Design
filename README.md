# CD276 IgG and Nanobody Design Using PPIFlow

This repository contains configuration files, helper scripts, and summary outputs for computational design and evaluation of CD276/B7-H3-targeting IgG antibody and nanobody candidates.

## Workflow

1. Define the CD276 epitope
2. Generate antibody and nanobody backbones with PPIFlow
3. Run AbMPNN sequence optimisation
4. Run Flowpacker side-chain packing
5. Score and rank candidates with Prodigy
6. Validate selected candidates with AF3
7. Run 10 ns GROMACS MD for selected structures

## Repository Layout

- `configs/demo_cd276/`: PPIFlow, batch, and Flowpacker-related task and step YAML files
- `configs/mdp/`: GROMACS MDP parameter files
- `scripts/`: pipeline helper scripts for batch runs, ranking export, and MD setup
- `results/prodigy/`: Prodigy scoring and ranked CSV summaries
- `results/af3_selected/`: selected AF3 input structures
- `results/md_analysis/`: MD analysis outputs, logs, and plots for selected candidates
- `ppiflow_af3_merged.yaml`: merged Conda environment for the PPIFlow and AF3-related workflow

## Key Files

- `results/prodigy/antibody_ranked_all.csv`
- `results/prodigy/nanobody_ranked_all.csv`
- `scripts/run_all_designs.sh`
- `scripts/run_prodigy_flowpacked_rank_all.py`
- `scripts/run_both_selected_md.sh`

## Notes

Large intermediate files, raw model outputs, and MD trajectory files should not be committed. The current `.gitignore` excludes `.DS_Store` files and `*.xtc` trajectories.
