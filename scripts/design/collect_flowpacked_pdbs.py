from pathlib import Path
import shutil
import csv

ROOT = Path("/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main")
SRC = ROOT / "output/abmpnn_flowpacker_designs"
OUT = ROOT / "output/final_flowpacked_designs"

if OUT.exists():
    shutil.rmtree(OUT)

(OUT / "antibody").mkdir(parents=True, exist_ok=True)
(OUT / "nanobody").mkdir(parents=True, exist_ok=True)

rows = []
counts = {"antibody": 0, "nanobody": 0}

def copy_group(pattern, group, prefix):
    files = sorted(SRC.glob(pattern))
    for f in files:
        counts[group] += 1
        dst = OUT / group / f"{prefix}_{counts[group]:04d}.pdb"
        shutil.copy2(f, dst)
        rows.append({
            "group": group,
            "new_file": dst.name,
            "source_path": str(f),
        })

copy_group("IgG_batch_*/stage1/flowpacker_output/run_1/*.pdb", "antibody", "IgG_FP")
copy_group("NB_batch_*/stage1/flowpacker_output/run_1/*.pdb", "nanobody", "NB_FP")

with (OUT / "flowpacked_file_map.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["group", "new_file", "source_path"])
    writer.writeheader()
    writer.writerows(rows)

print("Collected antibody:", counts["antibody"])
print("Collected nanobody:", counts["nanobody"])
print("Saved to:", OUT)
