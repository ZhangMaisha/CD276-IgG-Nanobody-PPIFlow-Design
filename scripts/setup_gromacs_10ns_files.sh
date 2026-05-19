#!/usr/bin/env bash
set -euo pipefail

mkdir -p mdp md_runs/antibody_rank5_10ns md_runs/nanobody_rank7_10ns

############################################
# ions.mdp
############################################
cat > mdp/ions.mdp <<'EOF'
; ions.mdp - used as input into grompp to generate ions.tpr
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = cutoff
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

############################################
# minim.mdp
############################################
cat > mdp/minim.mdp <<'EOF'
; minim.mdp - used as input into grompp to generate em.tpr
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

############################################
# nvt_short.mdp
############################################
cat > mdp/nvt_short.mdp <<'EOF'
title                   = Protein complex NVT equilibration
define                  = -DPOSRES
integrator              = md
nsteps                  = 16000      ; 32 ps with dt = 0.002 ps
dt                      = 0.002

nstxout                 = 500
nstvout                 = 500
nstenergy               = 500
nstlog                  = 500

continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300

pcoupl                  = no

pbc                     = xyz

gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
EOF

############################################
# npt_short.mdp
############################################
cat > mdp/npt_short.mdp <<'EOF'
title                   = Protein complex NPT equilibration
define                  = -DPOSRES
integrator              = md
nsteps                  = 20000      ; 40 ps with dt = 0.002 ps
dt                      = 0.002

nstxout                 = 500
nstvout                 = 500
nstenergy               = 500
nstlog                  = 500

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com

pbc                     = xyz
gen_vel                 = no
EOF

############################################
# md_10ns.mdp
############################################
cat > mdp/md_10ns.mdp <<'EOF'
title                   = Protein complex 10 ns production MD
integrator              = md
nsteps                  = 5000000    ; 10 ns with dt = 0.002 ps
dt                      = 0.002

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000
compressed-x-grps       = System

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0

coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

pbc                     = xyz
DispCorr                = EnerPres
gen_vel                 = no
EOF

############################################
# Main GROMACS runner
############################################
cat > run_gromacs_10ns.sh <<'EOF'
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
EOF

chmod +x run_gromacs_10ns.sh

############################################
# Run both systems
############################################
cat > run_both_selected_md.sh <<'EOF'
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
EOF

chmod +x run_both_selected_md.sh

############################################
# R plotting script adapted from your old file
############################################
cat > MD_plots_auto.R <<'EOF'
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript MD_plots_auto.R <md_run_dir>")
}

setwd(args[1])
dir.create("plots", showWarnings = FALSE)

read_xvg <- function(file, skip_default) {
  read.table(file, sep = "", header = FALSE, skip = skip_default,
             na.strings = "", stringsAsFactors = FALSE)
}

save_plot <- function(p, filename, w=7, h=5) {
  ggsave(file.path("plots", filename), p, width=w, height=h, dpi=300)
}

# Potential
potential <- read_xvg("potential.xvg", 24)
p1 <- ggplot(potential, aes(x = V1, y = V2)) +
  geom_line() +
  geom_point(size=0.5) +
  labs(x = "Energy minimization step",
       y = bquote("Potential energy (kJ "*mol^-1*")")) +
  ggtitle("Energy minimization") +
  theme_bw()
save_plot(p1, "potential.png")

# Temperature
temperature <- read_xvg("temperature.xvg", 24)
temperature$average10ps <- NA
if (nrow(temperature) >= 10) {
  temperature$average10ps[10:nrow(temperature)] <- sapply(
    10:nrow(temperature),
    function(x) mean(temperature$V2[(x-9):x])
  )
}
p2 <- ggplot(temperature, aes(x = V1, y = V2)) +
  geom_line() +
  geom_point(size=0.5) +
  geom_line(aes(y = average10ps, col = "Running average 10 ps")) +
  labs(x = "Time (ps)", y = "Temperature (K)") +
  ggtitle("Temperature, NVT equilibration") +
  theme_bw() +
  theme(legend.title = element_blank())
save_plot(p2, "temperature.png")

# Pressure
pressure <- read_xvg("pressure.xvg", 24)
pressure$average10ps <- NA
if (nrow(pressure) >= 10) {
  pressure$average10ps[10:nrow(pressure)] <- sapply(
    10:nrow(pressure),
    function(x) mean(pressure$V2[(x-9):x])
  )
}
p3 <- ggplot(pressure, aes(x = V1, y = V2)) +
  geom_line() +
  geom_point(size=0.5) +
  geom_line(aes(y = average10ps, col = "Running average 10 ps")) +
  labs(x = "Time (ps)", y = "Pressure (bar)") +
  ggtitle("Pressure, NPT equilibration") +
  theme_bw() +
  theme(legend.title = element_blank())
save_plot(p3, "pressure.png")

# Density
density <- read_xvg("density.xvg", 24)
density$average10ps <- NA
if (nrow(density) >= 10) {
  density$average10ps[10:nrow(density)] <- sapply(
    10:nrow(density),
    function(x) mean(density$V2[(x-9):x])
  )
}
p4 <- ggplot(density, aes(x = V1, y = V2)) +
  geom_line() +
  geom_point(size=0.5) +
  geom_line(aes(y = average10ps, col = "Running average 10 ps")) +
  labs(x = "Time (ps)", y = bquote("Density (kg "*m^-3*")")) +
  ggtitle("Density, NPT equilibration") +
  theme_bw() +
  theme(legend.title = element_blank())
save_plot(p4, "density.png")

# RMSD
rmsd_equilibrated <- read_xvg("rmsd.xvg", 18)
rmsd_xtal <- read_xvg("rmsd_xtal.xvg", 18)
rmsd <- rmsd_equilibrated
names(rmsd) <- c("time", "equilibrated")
rmsd$original <- rmsd_xtal$V2

p5 <- ggplot(rmsd, aes(x = time)) +
  geom_line(aes(y = equilibrated, col = "Ref: production start")) +
  geom_line(aes(y = original, col = "Ref: minimized structure")) +
  labs(x = "Time (ns)", y = "Backbone RMSD (nm)") +
  ggtitle("Backbone RMSD") +
  theme_bw() +
  theme(legend.title = element_blank())
save_plot(p5, "rmsd.png")

# Radius of gyration
gyration <- read_xvg("gyrate.xvg", 27)
p6 <- ggplot(gyration, aes(x = V1/1000, y = V2)) +
  geom_line() +
  geom_point(size=0.5) +
  labs(x = "Time (ns)", y = bquote(R[g]*" (nm)")) +
  ggtitle("Radius of gyration") +
  theme_bw()
save_plot(p6, "gyrate.png")

message("Plots saved in: ", file.path(getwd(), "plots"))
EOF

echo "All MD files created."
