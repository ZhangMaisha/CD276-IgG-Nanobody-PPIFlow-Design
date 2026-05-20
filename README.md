# CD276 IgG and Nanobody Design Using PPIFlow

This repository contains configuration files, selected input structures, helper scripts, and summary outputs for computational design and evaluation of CD276/B7-H3-targeting IgG antibody and nanobody candidates. It depends on the upstream PPIFlow repository: https://github.com/Mingchenchen/PPIFlow

## Workflow

1. Define the CD276 epitope
2. Generate antibody and nanobody backbones with PPIFlow
3. Run AbMPNN sequence optimisation
4. Run Flowpacker side-chain packing
5. Score and rank candidates with Prodigy
6. Validate selected candidates with AF3
7. Run MD for selected structures and retain the continuous 0-50 ns summary as the primary MD result

## Repository Layout

- `configs/examples/`: minimal demo task and step YAML files kept as readable reference examples
- `configs/generated_batches/ppiflow/`: generated batch YAML snapshots for stage-1 PPIFlow design runs
- `configs/generated_batches/abmpnn_flowpacker/`: generated batch YAML snapshots for AbMPNN + Flowpacker continuation runs
- `configs/mdp/`: GROMACS MDP parameter files
- `input/epitope_inputs/`: input structures used for epitope-support analyses such as DiscoTope
- `input/af3_batches/`: AF3 batch JSON inputs for top-ranked candidate validation
- `input/af3_selected/`: selected structures used as AF3 inputs
- `input/md_inputs/`: selected structures prepared as MD starting inputs
- `scripts/design/`: design-generation and candidate-export helpers for PPIFlow, AbMPNN, and Flowpacker stages
- `scripts/prodigy/`: Prodigy scoring, recheck parsing, and plotting helpers
- `scripts/md/`: MD setup, extension, and post-processing helpers
- `scripts/archive/`: older helper scripts kept for historical reference rather than the primary repository workflow
- `results/epitope_analysis/discotope/`: DiscoTope output used to support epitope definition
- `results/prodigy/initial_screening/`: Prodigy ranking summaries for the full flowpacked candidate sets
- `results/prodigy/af3_recheck/`: AF3 recheck summaries, source tables, and derived figures for selected candidates
- `results/md_analysis/primary_0_50ns/`: primary reported MD result, including summary tables and 0-50 ns comparison figures
- `ppiflow_af3_merged.yaml`: merged Conda environment for the PPIFlow and AF3-related workflow

## Key Files

- `input/epitope_inputs/fold_cd276_monomer_model_0.pdb`
- `input/af3_batches/antibody_lt_minus9_rank002_030.json`
- `input/af3_batches/nanobody_lt_minus9_rank001_030.json`
- `results/prodigy/initial_screening/antibody_ranked_all.csv`
- `results/prodigy/initial_screening/nanobody_ranked_all.csv`
- `results/epitope_analysis/discotope/file_A_discotope3.csv`
- `results/prodigy/af3_recheck/final_AF3_prodigy_summary.csv`
- `results/prodigy/af3_recheck/figures/Figure2A_binding_energy.png`
- `results/prodigy/af3_recheck/source_tables/PRODIGY_manual_summary.xlsx`
- `results/md_analysis/primary_0_50ns/summary/MD_metric_summary.csv`
- `scripts/design/run_all_designs.sh`
- `scripts/prodigy/run_prodigy_flowpacked_rank_all.py`
- `scripts/prodigy/plot_prodigy_recheck.R`
- `scripts/prodigy/parse_final_prodigy.py`
- `scripts/md/analyze_md_extra_metrics.sh`
- `scripts/md/analyze_interface_by_index.sh`
- `scripts/md/analyze_0_50ns_advanced.sh`
- `scripts/md/extend_both_selected_10ns.sh`
- `scripts/md/resume_both_to_50ns_gpu_safe.sh`
- `scripts/md/run_both_selected_md.sh`

## Notes

This repository does not vendor the PPIFlow source code. Please install or clone PPIFlow separately from `https://github.com/Mingchenchen/PPIFlow` when reproducing the workflow.

Large intermediate files, raw model outputs, and MD trajectory/state files should not be committed. The repository keeps lightweight summary tables, selected inputs, XVG analysis outputs, and combined plots, while excluding full AF3 download bundles, full MD trajectories, raw GROMACS runtime energy/log files, and archived 10 ns-only MD result folders.

Several execution scripts are preserved as historical workflow helpers and still assume the original external workspace layout under `/teams/.../PPIFlow-main`. They document how the runs were carried out, but they are not standalone scripts for this export repository.

Additional helper scripts are included for MD post-processing, interface-specific analysis, Prodigy summary parsing, and GPU-based trajectory extension to 50 ns.

`scripts/prodigy/plot_prodigy_recheck.R` expects the manually consolidated Prodigy workbook to be available in the working directory when the script is run. A copy is retained under `results/prodigy/af3_recheck/source_tables/`.

`scripts/archive/10ns_md_analysis.R` is retained as an older plotting script from the broader MD analysis workspace, but the corresponding 10 ns result folders are archived outside this repository because the continuous 0-50 ns MD summary is treated as the primary reported MD result.
