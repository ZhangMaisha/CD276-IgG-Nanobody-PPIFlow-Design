#!/usr/bin/env bash
set -euo pipefail

echo "===================================="
echo "START REMAINING AbMPNN + Flowpacker"
date
echo "===================================="

# antibody batch 2-12，因为 batch 1 已经完成
for i in $(seq 2 12)
do
    echo
    echo "####################################"
    echo "RUN ANTIBODY BATCH $i"
    echo "####################################"
    ./run_abmpnn_flowpacker_batch.sh antibody "$i"
done

# nanobody batch 1-12
for i in $(seq 1 12)
do
    echo
    echo "####################################"
    echo "RUN NANOBODY BATCH $i"
    echo "####################################"
    ./run_abmpnn_flowpacker_batch.sh nanobody "$i"
done

echo
echo "===================================="
echo "ALL REMAINING AbMPNN + Flowpacker FINISHED"
date
echo "===================================="

echo "Antibody flowpacked count:"
find output/abmpnn_flowpacker_designs/IgG_batch_* -path "*stage1/flowpacker_output/run_1/*.pdb" | wc -l

echo "Nanobody flowpacked count:"
find output/abmpnn_flowpacker_designs/NB_batch_* -path "*stage1/flowpacker_output/run_1/*.pdb" | wc -l
