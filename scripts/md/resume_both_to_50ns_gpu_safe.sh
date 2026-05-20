#!/usr/bin/env bash
set -euo pipefail

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
GMX="/root/miniconda3/envs/gmx_cuda/bin/gmx"
MDP="$ROOT/mdp/md_extend_10ns.mdp"
TARGET_STEP=5000000

GPU_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu -ntmpi 1 -ntomp 16"

last_step () {
    local LOG="$1"
    if [ ! -f "$LOG" ]; then
        echo 0
        return
    fi

    awk '
    /Step[[:space:]]+Time/ {
        getline
        if ($1 ~ /^[0-9]+$/) step=$1
    }
    END {
        if (step == "") print 0;
        else print step;
    }' "$LOG"
}

is_segment_complete () {
    local NEXT="$1"
    local STEP
    STEP=$(last_step "md_${NEXT}ns.log")

    if [ "$STEP" -ge "$TARGET_STEP" ] \
       && [ -f "md_${NEXT}ns.xtc" ] \
       && [ -f "md_${NEXT}ns.gro" ] \
       && [ -f "md_${NEXT}ns.cpt" ]; then
        return 0
    else
        return 1
    fi
}

postprocess_segment () {
    local NEXT="$1"

    echo "Post-processing md_${NEXT}ns"

    echo -e "Protein\nSystem" | "$GMX" trjconv \
      -s "md_${NEXT}ns.tpr" \
      -f "md_${NEXT}ns.xtc" \
      -o "md_${NEXT}ns_noPBC.xtc" \
      -pbc mol \
      -center

    echo -e "Backbone\nBackbone" | "$GMX" rms \
      -s "md_${NEXT}ns.tpr" \
      -f "md_${NEXT}ns_noPBC.xtc" \
      -o "rmsd_${NEXT}ns.xvg" \
      -tu ns

    echo -e "Backbone\nBackbone" | "$GMX" rms \
      -s em.tpr \
      -f "md_${NEXT}ns_noPBC.xtc" \
      -o "rmsd_${NEXT}ns_xtal.xvg" \
      -tu ns

    echo "Protein" | "$GMX" gyrate \
      -s "md_${NEXT}ns.tpr" \
      -f "md_${NEXT}ns_noPBC.xtc" \
      -o "gyrate_${NEXT}ns.xvg"
}

run_one_segment () {
    local RUN_DIR="$1"
    local NAME="$2"
    local NEXT="$3"
    local PREV=$((NEXT - 10))

    echo
    echo "======================================"
    echo "$NAME: ${PREV} ns -> ${NEXT} ns"
    echo "Directory: $RUN_DIR"
    echo "Start: $(date)"
    echo "======================================"

    cd "$RUN_DIR"

    local STEP
    STEP=$(last_step "md_${NEXT}ns.log")
    echo "Current md_${NEXT}ns step: $STEP / $TARGET_STEP"

    if is_segment_complete "$NEXT"; then
        echo "$NAME md_${NEXT}ns already complete. Skipping MD."

        if [ ! -f "rmsd_${NEXT}ns.xvg" ] || [ ! -f "gyrate_${NEXT}ns.xvg" ]; then
            echo "Analysis files missing, running post-processing only."
            postprocess_segment "$NEXT"
        fi

        return 0
    fi

    if [ ! -f "md_${PREV}ns.gro" ]; then
        echo "ERROR: md_${PREV}ns.gro not found in $RUN_DIR"
        exit 1
    fi

    if [ ! -f "md_${PREV}ns.cpt" ]; then
        echo "ERROR: md_${PREV}ns.cpt not found in $RUN_DIR"
        exit 1
    fi

    if [ ! -f "md_${NEXT}ns.tpr" ]; then
        echo "Generating md_${NEXT}ns.tpr from md_${PREV}ns.gro + md_${PREV}ns.cpt"

        "$GMX" grompp \
          -f "$MDP" \
          -c "md_${PREV}ns.gro" \
          -t "md_${PREV}ns.cpt" \
          -p topol.top \
          -o "md_${NEXT}ns.tpr" \
          -maxwarn 2
    else
        echo "Found existing md_${NEXT}ns.tpr"
    fi

    if [ -f "md_${NEXT}ns.cpt" ]; then
        echo "Continuing from checkpoint md_${NEXT}ns.cpt"

        "$GMX" mdrun \
          -deffnm "md_${NEXT}ns" \
          -cpi "md_${NEXT}ns.cpt" \
          -append \
          $GPU_FLAGS
    else
        echo "Starting new segment md_${NEXT}ns"

        "$GMX" mdrun \
          -deffnm "md_${NEXT}ns" \
          $GPU_FLAGS
    fi

    STEP=$(last_step "md_${NEXT}ns.log")
    echo "After mdrun, md_${NEXT}ns step: $STEP / $TARGET_STEP"

    if ! is_segment_complete "$NEXT"; then
        echo "ERROR: md_${NEXT}ns still not complete. Stop before post-processing."
        tail -n 80 "md_${NEXT}ns.log" || true
        exit 1
    fi

    postprocess_segment "$NEXT"

    echo "DONE $NAME to ${NEXT} ns at $(date)"
}

run_both_segment () {
    local NEXT="$1"

    run_one_segment "$ROOT/md_runs/nanobody_rank7_10ns" "nanobody_rank7" "$NEXT"
    run_one_segment "$ROOT/md_runs/antibody_rank5_10ns" "antibody_rank5" "$NEXT"
}

cd "$ROOT"

echo "======================================"
echo "START resume both systems to 50 ns"
echo "Start time: $(date)"
echo "======================================"

for NEXT in 20 30 40 50
do
    echo
    echo "######################################"
    echo "EXTENDING BOTH SYSTEMS TO ${NEXT} ns"
    echo "######################################"

    run_both_segment "$NEXT"
done

echo
echo "======================================"
echo "ALL DONE: both systems reached 50 ns"
echo "End time: $(date)"
echo "======================================"
