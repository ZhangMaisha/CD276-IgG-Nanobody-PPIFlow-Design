#!/usr/bin/env bash
set -euo pipefail

cd /teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main

mkdir -p md_runs/antibody_rank5_10ns
mkdir -p md_runs/nanobody_rank7_10ns

cp md_inputs/af3_selected/antibody_rank5/input.pdb md_runs/antibody_rank5_10ns/input.pdb
cp md_inputs/af3_selected/nanobody_rank7/input.pdb md_runs/nanobody_rank7_10ns/input.pdb

echo "First run nanobody:"
./run_gromacs_10ns.sh md_runs/nanobody_rank7_10ns 2>&1 | tee md_runs/nanobody_rank7_10ns/run.log

echo "Then run antibody:"
./run_gromacs_10ns.sh md_runs/antibody_rank5_10ns 2>&1 | tee md_runs/antibody_rank5_10ns/run.log
