#!/usr/bin/env bash
set -euo pipefail

NEXT=${1:?Usage: ./extend_both_selected_10ns.sh 20}

PREV=$((NEXT - 10))

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
MDP="$ROOT/mdp/md_extend_10ns.mdp"
GMX="/root/miniconda3/envs/gmx_cuda/bin/gmx"

run_one () {
    local RUN_DIR=$1
    local NAME=$2

    echo "======================================"
    echo "Running $NAME: ${PREV} ns -> ${NEXT} ns"
    echo "Directory: $RUN_DIR"
    echo "Start: $(date)"
    echo "GROMACS: $GMX"
    echo "======================================"

    cd "$RUN_DIR"

    if [ ! -f "md_${PREV}ns.gro" ]; then
        echo "ERROR: md_${PREV}ns.gro not found"
        exit 1
    fi

    if [ ! -f "md_${PREV}ns.cpt" ]; then
        echo "ERROR: md_${PREV}ns.cpt not found"
        exit 1
    fi

    # 重新跑这一段，先删除旧的/未完成的 NEXT 段输出
    rm -f md_${NEXT}ns.tpr \
          md_${NEXT}ns.xtc \
          md_${NEXT}ns.trr \
          md_${NEXT}ns.edr \
          md_${NEXT}ns.log \
          md_${NEXT}ns.gro \
          md_${NEXT}ns.cpt \
          md_${NEXT}ns_prev.cpt \
          md_${NEXT}ns_pullx.xvg \
          md_${NEXT}ns_pullf.xvg

    "$GMX" grompp \
      -f "$MDP" \
      -c "md_${PREV}ns.gro" \
      -t "md_${PREV}ns.cpt" \
      -p topol.top \
      -o "md_${NEXT}ns.tpr" \
      -maxwarn 2

    "$GMX" mdrun \
      -deffnm "md_${NEXT}ns" \
      -nb gpu \
      -pme gpu \
      -bonded gpu \
      -update gpu \
      -ntmpi 1 \
      -ntomp 16

    echo "DONE $NAME ${NEXT} ns at $(date)"
}

run_one "$ROOT/md_runs/nanobody_rank7_10ns" "nanobody_rank7"
run_one "$ROOT/md_runs/antibody_rank5_10ns" "antibody_rank5"

echo "======================================"
echo "Both systems finished to ${NEXT} ns"
echo "End: $(date)"
echo "======================================"
