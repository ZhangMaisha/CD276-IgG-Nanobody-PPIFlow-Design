#!/usr/bin/env bash
set -u

NS=${1:-20}

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
GMX="/root/miniconda3/envs/gmx_cuda/bin/gmx"

analyze_one () {
    local RUN_DIR="$1"
    local NAME="$2"
    local TYPE="$3"

    echo
    echo "======================================"
    echo "Analyzing $NAME, ${NS} ns"
    echo "Directory: $RUN_DIR"
    echo "======================================"

    cd "$RUN_DIR" || exit 1

    local TPR="md_${NS}ns.tpr"
    local XTC="md_${NS}ns.xtc"
    local NOPBC="md_${NS}ns_noPBC.xtc"

    if [ ! -f "$TPR" ]; then
        echo "Missing $TPR, skip $NAME"
        return
    fi

    if [ ! -f "$XTC" ]; then
        echo "Missing $XTC, skip $NAME"
        return
    fi

    # 1. Remove PBC if missing
    if [ ! -f "$NOPBC" ]; then
        echo "Generating no-PBC trajectory..."
        echo -e "Protein\nSystem" | "$GMX" trjconv \
            -s "$TPR" \
            -f "$XTC" \
            -o "$NOPBC" \
            -pbc mol \
            -center
    fi

    # 2. RMSD
    echo "RMSD..."
    echo -e "Backbone\nBackbone" | "$GMX" rms \
        -s "$TPR" \
        -f "$NOPBC" \
        -o "rmsd_${NS}ns.xvg" \
        -tu ns

    # 3. RMSD vs minimized structure
    if [ -f em.tpr ]; then
        echo "RMSD vs minimized structure..."
        echo -e "Backbone\nBackbone" | "$GMX" rms \
            -s em.tpr \
            -f "$NOPBC" \
            -o "rmsd_${NS}ns_xtal.xvg" \
            -tu ns
    fi

    # 4. Radius of gyration
    echo "Rg..."
    echo "Protein" | "$GMX" gyrate \
        -s "$TPR" \
        -f "$NOPBC" \
        -o "gyrate_${NS}ns.xvg"

    # 5. RMSF
    echo "RMSF..."
    echo "Backbone" | "$GMX" rmsf \
        -s "$TPR" \
        -f "$NOPBC" \
        -o "rmsf_${NS}ns.xvg" \
        -res

    # 6. SASA
    echo "SASA..."
    echo "Protein" | "$GMX" sasa \
        -s "$TPR" \
        -f "$NOPBC" \
        -o "sasa_${NS}ns.xvg" \
        -or "sasa_residue_${NS}ns.xvg"

    # 7. Interface contacts / minimum distance using GROMACS selection
    # IgG: antibody = chain A or B, antigen = chain C
    # Nanobody: binder = chain A, antigen = chain C

    if [ "$TYPE" = "antibody" ]; then
        BINDER='(chain "A" or chain "B")'
    else
        BINDER='chain "A"'
    fi

    ANTIGEN='chain "C"'

    echo "Minimum distance between binder and CD276..."
    "$GMX" pairdist \
        -s "$TPR" \
        -f "$NOPBC" \
        -ref "$BINDER" \
        -sel "$ANTIGEN" \
        -o "mindist_${NS}ns.xvg" \
        2>&1 | tee "mindist_${NS}ns.log"

    echo "Interface contacts within 0.45 nm..."
    "$GMX" select \
        -s "$TPR" \
        -f "$NOPBC" \
        -select "${BINDER} and within 0.45 of ${ANTIGEN}" \
        -os "contacts_${NS}ns.xvg" \
        2>&1 | tee "contacts_${NS}ns.log"

    echo "Hydrogen bonds between binder and CD276..."
    "$GMX" hbond \
        -f "$NOPBC" \
        -s "$TPR" \
        -r "$BINDER" \
        -t "$ANTIGEN" \
        -num "hbonds_${NS}ns.xvg" \
        2>&1 | tee "hbonds_${NS}ns.log"

    echo
    echo "Generated files for $NAME:"
    ls -lh \
        "rmsd_${NS}ns.xvg" \
        "rmsd_${NS}ns_xtal.xvg" \
        "gyrate_${NS}ns.xvg" \
        "rmsf_${NS}ns.xvg" \
        "sasa_${NS}ns.xvg" \
        "sasa_residue_${NS}ns.xvg" \
        "mindist_${NS}ns.xvg" \
        "contacts_${NS}ns.xvg" \
        "hbonds_${NS}ns.xvg" 2>/dev/null
}

analyze_one "$ROOT/md_runs/nanobody_rank7_10ns" "nanobody_rank7" "nanobody"
analyze_one "$ROOT/md_runs/antibody_rank5_10ns" "antibody_rank5" "antibody"

echo
echo "All available MD analyses finished for ${NS} ns."
