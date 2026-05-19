from pathlib import Path
import pandas as pd
import shutil

ROOT = Path("/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main")
RES = ROOT / "output/prodigy_flowpacked_results"
OUT = ROOT / "output/af3_lt_minus9_candidates"

if OUT.exists():
    shutil.rmtree(OUT)

(OUT / "antibody").mkdir(parents=True, exist_ok=True)
(OUT / "nanobody").mkdir(parents=True, exist_ok=True)
(OUT / "tables").mkdir(parents=True, exist_ok=True)

def export_group(group):
    csv_path = RES / f"{group}_ranked_all.csv"
    df = pd.read_csv(csv_path)
    df["dg_kcal_mol"] = pd.to_numeric(df["dg_kcal_mol"], errors="coerce")

    # strict less than -9
    sub = df[df["dg_kcal_mol"] < -9].copy()
    sub = sub.sort_values("dg_kcal_mol", ascending=True)
    sub.insert(0, "af3_candidate_id", [f"{group}_{i:04d}" for i in range(1, len(sub) + 1)])

    out_csv = OUT / "tables" / f"{group}_lt_minus9.csv"
    sub.to_csv(out_csv, index=False)

    for _, row in sub.iterrows():
        src = Path(row["path"])
        if not src.exists():
            print("MISSING PDB:", src)
            continue

        rank = int(row["rank"])
        dg = float(row["dg_kcal_mol"])
        dst_name = f"rank{rank:03d}_{dg:.2f}_{src.name}"
        dst = OUT / group / dst_name
        shutil.copy2(src, dst)

    print(f"{group}: {len(sub)} candidates with ΔG < -9")
    print("Table:", out_csv)
    print("PDB dir:", OUT / group)

export_group("antibody")
export_group("nanobody")
