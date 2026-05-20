from pathlib import Path
import subprocess
import csv
import re
import pandas as pd

ROOT = Path("/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main")
PDB_DIR = ROOT / "output/final_flowpacked_designs"
OUT = ROOT / "output/prodigy_flowpacked_results"
OUT.mkdir(parents=True, exist_ok=True)

def parse_value(pattern, text):
    m = re.search(pattern, text, flags=re.IGNORECASE)
    if not m:
        return ""
    return m.group(1)

def parse_stdout(stdout):
    dg = parse_value(r"Predicted binding affinity.*?:\s*([-\d\.]+)", stdout)
    kd = parse_value(r"Predicted dissociation constant.*?:\s*([-\deE\.]+)", stdout)
    contacts = parse_value(r"No\. of intermolecular contacts:\s*([-\d\.]+)", stdout)
    charged_charged = parse_value(r"No\. of charged-charged contacts:\s*([-\d\.]+)", stdout)
    charged_polar = parse_value(r"No\. of charged-polar contacts:\s*([-\d\.]+)", stdout)
    charged_apolar = parse_value(r"No\. of charged-apolar contacts:\s*([-\d\.]+)", stdout)
    polar_polar = parse_value(r"No\. of polar-polar contacts:\s*([-\d\.]+)", stdout)
    apolar_polar = parse_value(r"No\. of apolar-polar contacts:\s*([-\d\.]+)", stdout)
    apolar_apolar = parse_value(r"No\. of apolar-apolar contacts:\s*([-\d\.]+)", stdout)

    return {
        "dg_kcal_mol": dg,
        "kd_M": kd,
        "intermolecular_contacts": contacts,
        "charged_charged_contacts": charged_charged,
        "charged_polar_contacts": charged_polar,
        "charged_apolar_contacts": charged_apolar,
        "polar_polar_contacts": polar_polar,
        "apolar_polar_contacts": apolar_polar,
        "apolar_apolar_contacts": apolar_apolar,
    }

def run_group(group, selection):
    pdbs = sorted((PDB_DIR / group).glob("*.pdb"))
    rows = []

    print(f"\nRunning Prodigy for {group}: {len(pdbs)} PDBs")

    for i, pdb in enumerate(pdbs, 1):
        print(f"{group} {i}/{len(pdbs)} {pdb.name}", flush=True)

        cmd = ["prodigy", str(pdb), "--selection"] + selection
        res = subprocess.run(cmd, capture_output=True, text=True)

        parsed = parse_stdout(res.stdout)

        row = {
            "file": pdb.name,
            "path": str(pdb),
            "ok": res.returncode == 0,
            "returncode": res.returncode,
            **parsed,
            "stdout": res.stdout,
            "stderr": res.stderr,
        }
        rows.append(row)

    score_csv = OUT / f"{group}_prodigy_scores_all.csv"
    with score_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    df = pd.DataFrame(rows)
    df["dg_kcal_mol"] = pd.to_numeric(df["dg_kcal_mol"], errors="coerce")
    df_ranked = df[df["ok"] == True].dropna(subset=["dg_kcal_mol"]).copy()
    df_ranked = df_ranked.sort_values("dg_kcal_mol", ascending=True)
    df_ranked.insert(0, "rank", range(1, len(df_ranked) + 1))

    ranked_csv = OUT / f"{group}_ranked_all.csv"
    df_ranked.to_csv(ranked_csv, index=False)

    print("Saved scores:", score_csv)
    print("Saved ranked:", ranked_csv)

run_group("antibody", ["A,B", "C"])
run_group("nanobody", ["A", "C"])
