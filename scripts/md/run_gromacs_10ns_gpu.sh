#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=${1:?Usage: ./run_gromacs_10ns_gpu.sh md_runs/antibody_rank5_10ns}

cd "$RUN_DIR"

GMX=gmx

echo "Running GPU GROMACS MD in: $(pwd)"
$GMX mdrun -version | grep -Ei "GPU|CUDA|OpenCL|SYCL" || true

FORCE_FIELD="oplsaa"
WATER_MODEL="spce"

$GMX pdb2gmx \
    -f input.pdb \
    -o processed.gro \
    -p topol.top \
    -i posre.itp \
    -ff "$FORCE_FIELD" \
    -water "$WATER_MODEL" \
    -ignh

$GMX editconf \
    -f processed.gro \
    -o newbox.gro \
    -c \
    -d 1.0 \
    -bt cubic

$GMX solvate \
    -cp newbox.gro \
    -cs spc216.gro \
    -o solv.gro \
    -p topol.top

$GMX grompp \
    -f ../../mdp/ions.mdp \
    -c solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 2

echo "SOL" | $GMX genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p topol.top \
    -pname NA \
    -nname CL \
    -neutral

$GMX grompp \
    -f ../../mdp/minim.mdp \
    -c solv_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 2

$GMX mdrun -v -deffnm em

echo "Potential" | $GMX energy \
    -f em.edr \
    -o potential.xvg

$GMX grompp \
    -f ../../mdp/nvt_short.mdp \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr \
    -maxwarn 2

$GMX mdrun -deffnm nvt -nb gpu -pme gpu

echo "Temperature" | $GMX energy \
    -f nvt.edr \
    -o temperature.xvg

$GMX grompp \
    -f ../../mdp/npt_short.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt.tpr \
    -maxwarn 2

$GMX mdrun -deffnm npt -nb gpu -pme gpu

echo "Pressure" | $GMX energy \
    -f npt.edr \
    -o pressure.xvg

echo "Density" | $GMX energy \
    -f npt.edr \
    -o density.xvg

$GMX grompp \
    -f ../../mdp/md_10ns.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -o md_10ns.tpr \
    -maxwarn 2

$GMX mdrun \
    -deffnm md_10ns \
    -nb gpu \
    -pme gpu

echo -e "Protein\nSystem" | $GMX trjconv \
    -s md_10ns.tpr \
    -f md_10ns.xtc \
    -o md_10ns_noPBC.xtc \
    -pbc mol \
    -center

echo -e "Backbone\nBackbone" | $GMX rms \
    -s md_10ns.tpr \
    -f md_10ns_noPBC.xtc \
    -o rmsd.xvg \
    -tu ns

echo -e "Backbone\nBackbone" | $GMX rms \
    -s em.tpr \
    -f md_10ns_noPBC.xtc \
    -o rmsd_xtal.xvg \
    -tu ns

echo "Protein" | $GMX gyrate \
    -s md_10ns.tpr \
    -f md_10ns_noPBC.xtc \
    -o gyrate.xvg

echo "DONE antibody GPU 10 ns MD"
