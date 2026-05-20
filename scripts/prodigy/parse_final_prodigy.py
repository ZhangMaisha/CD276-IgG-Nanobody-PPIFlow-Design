from pathlib import Path
import re
import csv

ROOT = Path("/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main")
IN = ROOT / "output/prodigy_final_recheck"
OUT = IN / "final_AF3_prodigy_summary.csv"

files = {
    "nanobody_rank7_AF3": IN / "nanobody_rank7_AF3_prodigy.txt",
    "antibody_rank5_AF3": IN / "antibody_rank5_AF3_prodigy.txt",
}

def find(pattern, text):
    m = re.search(pattern, text, flags=re.I)
    return m.group(1) if m else ""

rows = []

for name, fp in files.items():
    text = fp.read_text(errors="ignore")

    rows.append({
        "candidate": name,
        "intermolecular_contacts": find(r"No\. of intermolecular contacts:\s*([-\d\.]+)", text),
        "charged_charged_contacts": find(r"No\. of charged-charged contacts:\s*([-\d\.]+)", text),
        "charged_polar_contacts": find(r"No\. of charged-polar contacts:\s*([-\d\.]+)", text),
        "charged_apolar_contacts": find(r"No\. of charged-apolar contacts:\s*([-\d\.]+)", text),
        "polar_polar_contacts": find(r"No\. of polar-polar contacts:\s*([-\d\.]+)", text),
        "apolar_polar_contacts": find(r"No\. of apolar-polar contacts:\s*([-\d\.]+)", text),
        "apolar_apolar_contacts": find(r"No\. of apolar-apolar contacts:\s*([-\d\.]+)", text),
        "apolar_NIS_percent": find(r"Percentage of apolar NIS residues:\s*([-\d\.]+)", text),
        "charged_NIS_percent": find(r"Percentage of charged NIS residues:\s*([-\d\.]+)", text),
        "dg_kcal_mol": find(r"Predicted binding affinity.*?:\s*([-\d\.]+)", text),
        "kd_M": find(r"Predicted dissociation constant.*?:\s*([-\deE\.]+)", text),
    })

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print("Saved:", OUT)
