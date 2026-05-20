#!/usr/bin/env bash
set -u

NS=${1:-10}
ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
GMX="/root/miniconda3/envs/gmx_cuda/bin/gmx"

make_index () {
python - <<'PY'
from pathlib import Path

pdb = Path("input.pdb")
gro = Path("processed.gro")
out = Path("index_chain_auto.ndx")

if not pdb.exists():
    raise FileNotFoundError("input.pdb not found")
if not gro.exists():
    raise FileNotFoundError("processed.gro not found")

# 1. Count residues per chain from original input.pdb
chain_order = []
chain_residues = {}

for line in pdb.read_text(errors="ignore").splitlines():
    if not line.startswith("ATOM"):
        continue
    ch = line[21].strip()
    resi = line[22:26].strip()
    icode = line[26].strip()
    resn = line[17:20].strip()
    key = (resi, icode, resn)
    if ch not in chain_residues:
        chain_residues[ch] = []
        chain_order.append(ch)
    if key not in chain_residues[ch]:
        chain_residues[ch].append(key)

chain_counts = {ch: len(chain_residues[ch]) for ch in chain_order}
print("Chain residue counts from input.pdb:", chain_counts)

if "C" not in chain_counts:
    raise ValueError("Cannot find antigen chain C in input.pdb")

# 2. Parse processed.gro atom lines and reconstruct residue order
lines = gro.read_text(errors="ignore").splitlines()
atom_lines = lines[2:-1]

residue_atoms = []
current_key = None
current_atoms = []

for line in atom_lines:
    resnr = line[0:5].strip()
    resn = line[5:10].strip()
    atomnr = int(line[15:20].strip())
    key = (resnr, resn)

    if current_key is None:
        current_key = key

    if key != current_key:
        residue_atoms.append(current_atoms)
        current_atoms = []
        current_key = key

    current_atoms.append(atomnr)

if current_atoms:
    residue_atoms.append(current_atoms)

print("Residues reconstructed from processed.gro:", len(residue_atoms))

# 3. Assign reconstructed residues to chains by original chain lengths
chain_atoms = {}
idx = 0
for ch in chain_order:
    n = chain_counts[ch]
    atoms = []
    for res_atoms in residue_atoms[idx:idx+n]:
        atoms.extend(res_atoms)
    chain_atoms[ch] = atoms
    idx += n

print("Assigned atom counts:", {k: len(v) for k, v in chain_atoms.items()})

# 4. Define binder and antigen groups
if "B" in chain_atoms:
    binder_atoms = chain_atoms.get("A", []) + chain_atoms.get("B", [])
else:
    binder_atoms = chain_atoms.get("A", [])

antigen_atoms = chain_atoms.get("C", [])

if not binder_atoms:
    raise ValueError("No binder atoms found")
if not antigen_atoms:
    raise ValueError("No antigen chain C atoms found")

def write_group(f, name, atoms):
    f.write(f"[ {name} ]\n")
    for i in range(0, len(atoms), 15):
        f.write(" ".join(str(a) for a in atoms[i:i+15]) + "\n")
    f.write("\n")

with out.open("w") as f:
    for ch, atoms in chain_atoms.items():
        write_group(f, f"Chain_{ch}", atoms)
    write_group(f, "Binder", binder_atoms)
    write_group(f, "Antigen_C", antigen_atoms)
    write_group(f, "Complex", binder_atoms + antigen_atoms)

print("Wrote:", out)
PY
}

analyze_one () {
    local RUN_DIR="$1"
    local NAME="$2"

    echo
    echo "======================================"
    echo "Analyzing $NAME, ${NS} ns"
    echo "Directory: $RUN_DIR"
    echo "======================================"

    cd "$RUN_DIR" || exit 1

    if [ ! -f "md_${NS}ns.tpr" ]; then
        echo "Missing md_${NS}ns.tpr, skip"
        return
    fi

    if [ ! -f "md_${NS}ns.xtc" ]; then
        echo "Missing md_${NS}ns.xtc, skip"
        return
    fi

    if [ ! -f "md_${NS}ns_noPBC.xtc" ]; then
        echo "Generating noPBC trajectory..."
        echo -e "Protein\nSystem" | "$GMX" trjconv \
          -s "md_${NS}ns.tpr" \
          -f "md_${NS}ns.xtc" \
          -o "md_${NS}ns_noPBC.xtc" \
          -pbc mol \
          -center
    fi

    echo "Creating index_chain_auto.ndx..."
    make_index

    echo "Minimum distance and interface contacts..."
    echo -e "Binder\nAntigen_C" | "$GMX" mindist \
      -s "md_${NS}ns.tpr" \
      -f "md_${NS}ns_noPBC.xtc" \
      -n index_chain_auto.ndx \
      -od "mindist_${NS}ns.xvg" \
      -on "contacts_${NS}ns.xvg" \
      -d 0.45 \
      2>&1 | tee "mindist_contacts_${NS}ns.log"

    echo "Hydrogen bonds..."
    echo -e "Binder\nAntigen_C" | "$GMX" hbond \
      -s "md_${NS}ns.tpr" \
      -f "md_${NS}ns_noPBC.xtc" \
      -n index_chain_auto.ndx \
      -num "hbonds_${NS}ns.xvg" \
      2>&1 | tee "hbonds_${NS}ns.log" || true

    echo "Generated:"
    ls -lh "mindist_${NS}ns.xvg" "contacts_${NS}ns.xvg" "hbonds_${NS}ns.xvg" 2>/dev/null || true
}

analyze_one "$ROOT/md_runs/nanobody_rank7_10ns" "nanobody_rank7"
analyze_one "$ROOT/md_runs/antibody_rank5_10ns" "antibody_rank5"

echo
echo "Interface analyses finished for ${NS} ns."
