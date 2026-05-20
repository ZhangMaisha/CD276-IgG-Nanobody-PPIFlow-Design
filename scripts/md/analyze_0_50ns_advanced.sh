#!/usr/bin/env bash
set -euo pipefail

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
GMX="/root/miniconda3/envs/gmx_cuda/bin/gmx"
TARGET_STEP=5000000

SEGMENTS=(10 20 30 40 50)

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

check_segment_complete () {
    local SEG="$1"

    local LOG="md_${SEG}ns.log"
    local STEP
    STEP=$(last_step "$LOG")

    echo "Checking md_${SEG}ns: last step = $STEP / $TARGET_STEP"

    if [ "$STEP" -lt "$TARGET_STEP" ]; then
        echo "ERROR: md_${SEG}ns is not complete."
        exit 1
    fi

    for f in "md_${SEG}ns.tpr" "md_${SEG}ns.xtc" "md_${SEG}ns.gro" "md_${SEG}ns.cpt"; do
        if [ ! -f "$f" ]; then
            echo "ERROR: missing $f"
            exit 1
        fi
    done
}

make_noPBC_segment () {
    local SEG="$1"

    if [ -f "md_${SEG}ns_noPBC.xtc" ]; then
        echo "md_${SEG}ns_noPBC.xtc already exists."
        return
    fi

    echo "Generating md_${SEG}ns_noPBC.xtc"

    printf "Protein\nSystem\n" | "$GMX" trjconv \
      -s "md_${SEG}ns.tpr" \
      -f "md_${SEG}ns.xtc" \
      -o "md_${SEG}ns_noPBC.xtc" \
      -pbc mol \
      -center
}

shift_segment_time () {
    local SEG="$1"
    local T0="$2"

    if [ "$SEG" = "10" ]; then
        return
    fi

    local OUT="md_${SEG}ns_noPBC_shifted.xtc"

    if [ -f "$OUT" ]; then
        echo "$OUT already exists."
        return
    fi

    echo "Shifting md_${SEG}ns_noPBC.xtc to start at ${T0} ps"

    printf "System\n" | "$GMX" trjconv \
      -s "md_${SEG}ns.tpr" \
      -f "md_${SEG}ns_noPBC.xtc" \
      -o "$OUT" \
      -t0 "$T0"
}

make_continuous_trajectory () {
    if [ -f "md_0_50ns_noPBC.xtc" ]; then
        echo "md_0_50ns_noPBC.xtc already exists."
        return
    fi

    echo "Concatenating segments into md_0_50ns_noPBC.xtc"

    "$GMX" trjcat \
      -f md_10ns_noPBC.xtc \
         md_20ns_noPBC_shifted.xtc \
         md_30ns_noPBC_shifted.xtc \
         md_40ns_noPBC_shifted.xtc \
         md_50ns_noPBC_shifted.xtc \
      -o md_0_50ns_noPBC.xtc \
      -cat
}

make_chain_index () {
python3 - <<'PY'
from pathlib import Path

pdb = Path("input.pdb")
gro = Path("processed.gro")
out = Path("index_chain_auto.ndx")

if not pdb.exists():
    raise FileNotFoundError("input.pdb not found")
if not gro.exists():
    raise FileNotFoundError("processed.gro not found")

chain_order = []
chain_res = {}

for line in pdb.read_text(errors="ignore").splitlines():
    if not line.startswith("ATOM"):
        continue
    if line[12:16].strip() != "CA":
        continue

    chain = line[21].strip()
    if chain == "":
        chain = "blank"

    resi = line[22:26].strip()
    resn = line[17:20].strip()

    if chain not in chain_res:
        chain_res[chain] = []
        chain_order.append(chain)

    chain_res[chain].append((resi, resn))

print("Chains from input.pdb:")
for ch in chain_order:
    print(ch, len(chain_res[ch]))

lines = gro.read_text(errors="ignore").splitlines()
atom_lines = lines[2:-1]

residue_atoms = []
current = None
atoms = []

for line in atom_lines:
    resnr = line[0:5].strip()
    resn = line[5:10].strip()
    atomnr = int(line[15:20].strip())
    key = (resnr, resn)

    if current is None:
        current = key

    if key != current:
        residue_atoms.append(atoms)
        atoms = []
        current = key

    atoms.append(atomnr)

if atoms:
    residue_atoms.append(atoms)

print("Residues in processed.gro:", len(residue_atoms))

chain_atoms = {}
idx = 0

for ch in chain_order:
    n = len(chain_res[ch])
    selected = residue_atoms[idx:idx+n]
    flat = []
    for res_atoms in selected:
        flat.extend(res_atoms)
    chain_atoms[ch] = flat
    idx += n

print("Assigned atom counts:")
for ch in chain_order:
    print(ch, len(chain_atoms[ch]))

# AF3 Server convention in your current files:
# Nanobody:
#   A = CD276 antigen
#   B = nanobody
# Antibody:
#   A = CD276 antigen
#   B/C = antibody heavy/light
if set(["A", "B", "C"]).issubset(set(chain_atoms.keys())):
    binder_chains = ["B", "C"]
    antigen_chain = "A"
elif set(["A", "B"]).issubset(set(chain_atoms.keys())):
    binder_chains = ["B"]
    antigen_chain = "A"
else:
    raise ValueError(f"Unexpected chains: {chain_order}")

binder_atoms = []
for ch in binder_chains:
    binder_atoms.extend(chain_atoms[ch])

antigen_atoms = chain_atoms[antigen_chain]

def write_group(f, name, atoms):
    f.write(f"[ {name} ]\n")
    for i in range(0, len(atoms), 15):
        f.write(" ".join(str(a) for a in atoms[i:i+15]) + "\n")
    f.write("\n")

with out.open("w") as f:
    for ch in chain_order:
        write_group(f, f"Chain_{ch}", chain_atoms[ch])
    write_group(f, "Binder", binder_atoms)
    write_group(f, "Antigen", antigen_atoms)
    write_group(f, "Complex", binder_atoms + antigen_atoms)

print("Wrote:", out)
print("Binder chains:", binder_chains)
print("Antigen chain:", antigen_chain)
PY
}

analyze_continuous () {
    echo "RMSD 0–50 ns"
    printf "Backbone\nBackbone\n" | "$GMX" rms \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -o rmsd_0_50ns.xvg \
      -tu ns

    echo "RMSD vs minimized structure 0–50 ns"
    printf "Backbone\nBackbone\n" | "$GMX" rms \
      -s em.tpr \
      -f md_0_50ns_noPBC.xtc \
      -o rmsd_0_50ns_xtal.xvg \
      -tu ns

    echo "Radius of gyration 0–50 ns"
    printf "Protein\n" | "$GMX" gyrate \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -o gyrate_0_50ns.xvg

    echo "RMSF all backbone 0–50 ns"
    printf "Backbone\n" | "$GMX" rmsf \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -o rmsf_0_50ns.xvg \
      -res

    echo "SASA 0–50 ns"
    printf "Protein\n" | "$GMX" sasa \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -o sasa_0_50ns.xvg \
      -or sasa_residue_0_50ns.xvg

    echo "Creating chain index"
    make_chain_index

    echo "Minimum distance and interface contacts 0–50 ns"
    printf "Binder\nAntigen\n" | "$GMX" mindist \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -n index_chain_auto.ndx \
      -od mindist_0_50ns.xvg \
      -on contacts_0_50ns.xvg \
      -d 0.45 \
      2>&1 | tee mindist_contacts_0_50ns.log

    echo "Hydrogen bonds 0–50 ns"
    printf "Binder\nAntigen\n" | "$GMX" hbond \
      -s md_10ns.tpr \
      -f md_0_50ns_noPBC.xtc \
      -n index_chain_auto.ndx \
      -num hbonds_0_50ns.xvg \
      2>&1 | tee hbonds_0_50ns.log || true

    echo "Chain-specific RMSF 0–50 ns"
    for GROUP in Chain_A Chain_B Chain_C Binder Antigen
    do
        if grep -q "\[ $GROUP \]" index_chain_auto.ndx; then
            printf "$GROUP\n" | "$GMX" rmsf \
              -s md_10ns.tpr \
              -f md_0_50ns_noPBC.xtc \
              -n index_chain_auto.ndx \
              -o "rmsf_0_50ns_${GROUP}.xvg" \
              -res
        fi
    done
}

analyze_one_system () {
    local RUN_DIR="$1"
    local NAME="$2"

    echo
    echo "============================================================"
    echo "Analyzing continuous 0–50 ns for $NAME"
    echo "Directory: $RUN_DIR"
    echo "Start: $(date)"
    echo "============================================================"

    cd "$RUN_DIR"

    echo "Checking completed MD segments..."
    for SEG in "${SEGMENTS[@]}"
    do
        check_segment_complete "$SEG"
    done

    echo "Generating noPBC trajectories if missing..."
    for SEG in "${SEGMENTS[@]}"
    do
        make_noPBC_segment "$SEG"
    done

    echo "Generating time-shifted trajectories..."
    shift_segment_time 20 10000
    shift_segment_time 30 20000
    shift_segment_time 40 30000
    shift_segment_time 50 40000

    make_continuous_trajectory
    analyze_continuous

    echo
    echo "Generated files for $NAME:"
    ls -lh \
      md_0_50ns_noPBC.xtc \
      rmsd_0_50ns.xvg \
      rmsd_0_50ns_xtal.xvg \
      gyrate_0_50ns.xvg \
      rmsf_0_50ns.xvg \
      sasa_0_50ns.xvg \
      sasa_residue_0_50ns.xvg \
      mindist_0_50ns.xvg \
      contacts_0_50ns.xvg 2>/dev/null || true

    echo "Finished $NAME at $(date)"
}

analyze_one_system "$ROOT/md_runs/nanobody_rank7_10ns" "nanobody_rank7"
analyze_one_system "$ROOT/md_runs/antibody_rank5_10ns" "antibody_rank5"

echo
echo "============================================================"
echo "All continuous 0–50 ns advanced analyses finished."
echo "End: $(date)"
echo "============================================================"
