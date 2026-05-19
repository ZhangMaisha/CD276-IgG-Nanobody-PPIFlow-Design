#!/usr/bin/env bash
set -euo pipefail

echo "=================================="
echo "START ALL DESIGNS"
date
echo "=================================="

########################################
# NANOBODY 2-12
########################################

for i in $(seq 2 12)
do
    echo
    echo "##################################"
    echo "RUN NANOBODY BATCH $i"
    echo "##################################"

    ./run_ppiflow_batch.sh nanobody $i 25
done

########################################
# ANTIBODY 1-12
########################################

for i in $(seq 1 12)
do
    echo
    echo "##################################"
    echo "RUN ANTIBODY BATCH $i"
    echo "##################################"

    ./run_ppiflow_batch.sh antibody $i 25
done

echo
echo "=================================="
echo "ALL JOBS FINISHED"
date
echo "=================================="
