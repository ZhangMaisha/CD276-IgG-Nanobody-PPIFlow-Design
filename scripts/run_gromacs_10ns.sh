#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=${1:?Usage: ./run_gromacs_10ns.sh md_runs/nanobody_rank7_10ns}

cd "$RUN_DIR"

echo "=========================================="
echo "Running GROMACS 10 ns MD in: $(pwd)"
echo "Input PDB:"
ls -lh input.pdb
echo "=========================================="

# 自动找 GROMACS 命令
if command -v gmx >/dev/null 2>&1; then
    GMX=gmx
elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX=gmx_mpi
elif [ -x /software/gromacs-2024.3/build/bin/gmx_mpi ]; then
    GMX=/software/gromacs-2024.3/build/bin/gmx_mpi
elif [ -x /usr/local/gromacs/bin/gmx ]; then
    GMX=/usr/local/gromacs/bin/gmx
else
    echo "ERROR: Cannot find gmx or gmx_mpi"
    exit 1
fi

echo "Using GROMACS command: $GMX"
$GMX --version | head -20 || true

# 你老师旧脚本像是 OPLS + spce 风格
FORCE_FIELD="oplsaa"
WATER_MODEL="spce"

echo "Force field: $FORCE_FIELD"
echo "Water model: $WATER_MODEL"

# 1. Generate topology
$GMX pdb2gmx \
    -f input.pdb \
    -o processed.gro \
    -p topol.top \
    -i posre.itp \
    -ff "$FORCE_FIELD" \
    -water "$WATER_MODEL" \
    -ignh

# 2. Define the box
$GMX editconf \
    -f processed.gro \
    -o newbox.gro \
    -c \
    -d 1.0 \
    -bt cubic

# 3. Solvate
$GMX solvate \
    -cp newbox.gro \
    -cs spc216.gro \
    -o solv.gro \
    -p topol.top

# 4. Add ions
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

# 5. Energy minimization
$GMX grompp \
    -f ../../mdp/minim.mdp \
    -c solv_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 2

$GMX mdrun -v -deffnm em

# Energy potential
echo "Potential" | $GMX energy \
    -f em.edr \
    -o potential.xvg

# 6. NVT equilibration
$GMX grompp \
    -f ../../mdp/nvt_short.mdp \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr \
    -maxwarn 2

$GMX mdrun -deffnm nvt

echo "Temperature" | $GMX energy \
    -f nvt.edr \
    -o temperature.xvg

# 7. NPT equilibration
$GMX grompp \
    -f ../../mdp/npt_short.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt.tpr \
    -maxwarn 2

$GMX mdrun -deffnm npt

echo "Pressure" | $GMX energy \
    -f npt.edr \
    -o pressure.xvg

echo "Density" | $GMX energy \
    -f npt.edr \
    -o density.xvg

# 8. Production MD, 10 ns
$GMX grompp \
    -f ../../mdp/md_10ns.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -o md_10ns.tpr \
    -maxwarn 2

$GMX mdrun -deffnm md_10ns

# 9. Remove PBC
echo -e "Protein\nSystem" | $GMX trjconv \
    -s md_10ns.tpr \
    -f md_10ns.xtc \
    -o md_10ns_noPBC.xtc \
    -pbc mol \
    -center

# 10. RMSD relative to production start
echo -e "Backbone\nBackbone" | $GMX rms \
    -s md_10ns.tpr \
    -f md_10ns_noPBC.xtc \
    -o rmsd.xvg \
    -tu ns

# 11. RMSD relative to minimized structure
echo -e "Backbone\nBackbone" | $GMX rms \
    -s em.tpr \
    -f md_10ns_noPBC.xtc \
    -o rmsd_xtal.xvg \
    -tu ns

# 12. Radius of gyration
echo "Protein" | $GMX gyrate \
    -s md_10ns.tpr \
    -f md_10ns_noPBC.xtc \
    -o gyrate.xvg

echo "=========================================="
echo "DONE MD in: $(pwd)"
echo "Key outputs:"
ls -lh md_10ns.xtc md_10ns.gro rmsd.xvg rmsd_xtal.xvg gyrate.xvg
echo "=========================================="
