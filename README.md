# CD276 IgG and Nanobody Design Using PPIFlow

This repository contains scripts, configurations, and summary results for computational design and evaluation of CD276/B7-H3-targeting IgG antibody and nanobody candidates.

## Workflow

1. Epitope definition: CD276 residues C242-C248 (PGFSLAQ) and C282-C285 (FPDL)
2. PPIFlow generation: 300 IgG/scFv + 300 nanobody candidates
3. AbMPNN sequence optimisation
4. Flowpacker side-chain packing
5. Prodigy scoring and ranking
6. AF3 Server validation for candidates with predicted ΔG < -9 kcal/mol
7. GROMACS 10 ns MD for selected IgG and nanobody candidates

## Key result files

- `results/prodigy/antibody_ranked_all.csv`
- `results/prodigy/nanobody_ranked_all.csv`
- `results/af3_candidates/antibody_lt_minus9.csv`
- `results/af3_candidates/nanobody_lt_minus9.csv`
- `results/af3_json_batches/`
- `results/md_summary/`

## Notes

Large outputs, raw PDB batches, MD trajectories, model weights, checkpoints, Flowpacker `cluster.pth`, and AF3 model parameters are not included.
