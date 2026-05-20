#!/usr/bin/env bash
set -euo pipefail

TYPE=${1:?Need type: antibody or nanobody}
BATCH=${2:?Need batch number}

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
PYTHON="/root/miniconda3/envs/ppiflow_af3/bin/python"
PPI_CFG="$ROOT/tool/PPIFlow/configs/inference_nanobody.yaml"
FLOW_CFG="$ROOT/tools/flowpacker/config/inference/base.yaml"

mkdir -p "$ROOT/demo_cd276/fp_batches"
mkdir -p "$ROOT/output/abmpnn_flowpacker_designs"
mkdir -p "$ROOT/logs"

if [ "$TYPE" = "antibody" ]; then
    RAW="$ROOT/output/raw_designs/IgG_batch_${BATCH}"
    OUT="$ROOT/output/abmpnn_flowpacker_designs/IgG_batch_${BATCH}"
    NAME="IgG_batch_${BATCH}_abmpnn_flowpacker"
    GENTYPE="antibody"
    FRAMEWORK="$ROOT/example/target_and_framework_pdb/6nou_scfv_framework.pdb"
    CHAIN_LIST="A B"
    EXTRA_TASK='  light_chain: "B"
  cdr_length: "CDRH1,8-8,CDRH2,8-8,CDRH3,10-20,CDRL1,6-9,CDRL2,3-3,CDRL3,9-11"'
elif [ "$TYPE" = "nanobody" ]; then
    RAW="$ROOT/output/raw_designs/NB_batch_${BATCH}"
    OUT="$ROOT/output/abmpnn_flowpacker_designs/NB_batch_${BATCH}"
    NAME="NB_batch_${BATCH}_abmpnn_flowpacker"
    GENTYPE="nanobody"
    FRAMEWORK="$ROOT/example/target_and_framework_pdb/7eow_nanobody_framework.pdb"
    CHAIN_LIST="A"
    EXTRA_TASK='  cdr_length: "CDRH1,8-8,CDRH2,8-8,CDRH3,9-21"'
else
    echo "TYPE must be antibody or nanobody"
    exit 1
fi

RAW_PDB_DIR="$RAW/stage1/ppiflow_output"

if [ ! -d "$RAW_PDB_DIR" ]; then
    echo "ERROR: raw PPIFlow directory not found: $RAW_PDB_DIR"
    exit 1
fi

RAW_COUNT=$(find "$RAW_PDB_DIR" -name "*.pdb" | wc -l)
echo "Raw PDB count in $RAW_PDB_DIR = $RAW_COUNT"

if [ "$RAW_COUNT" -eq 0 ]; then
    echo "ERROR: no raw PDB found"
    exit 1
fi

TASK="$ROOT/demo_cd276/fp_batches/${NAME}.yaml"
STEPS="$ROOT/demo_cd276/fp_batches/${NAME}_steps.yaml"

echo "=========================================="
echo "TYPE  = $TYPE"
echo "BATCH = $BATCH"
echo "RAW   = $RAW"
echo "OUT   = $OUT"
echo "START = $(date)"
echo "=========================================="

rm -rf "$OUT"
mkdir -p "$OUT/stage1/ppiflow_output"
mkdir -p "$OUT/stage1/mpnn_pdbs"

# Copy raw PPIFlow PDBs. Do not use cleaned_designs for pipeline continuation.
cp "$RAW_PDB_DIR"/*.pdb "$OUT/stage1/ppiflow_output/"

# Because /teams does not support symlink, manually copy to mpnn_pdbs.
# Use lowercase names only, because AbMPNN/Flowpacker often writes lowercase link_name.
for f in "$RAW_PDB_DIR"/*.pdb; do
    base=$(basename "$f")
    lower=$(echo "$base" | tr '[:upper:]' '[:lower:]')
    cp "$f" "$OUT/stage1/mpnn_pdbs/$lower"
done

echo "Prepared mpnn_pdbs:"
find "$OUT/stage1/mpnn_pdbs" -name "*.pdb" | wc -l
find "$OUT/stage1/mpnn_pdbs" -name "*.pdb" | head

cat > "$TASK" <<YAML
task:
  name: "$NAME"
  gentype: "$GENTYPE"
  antigen_pdb: "$ROOT/input/CD276_C.pdb"
  antigen_chain: "C"
  specified_hotspots: "C242,C243,C244,C245,C246,C247,C248,C282,C283,C284,C285"
  framework_pdb: "$FRAMEWORK"
  heavy_chain: "A"
$EXTRA_TASK
  samples_per_target: $RAW_COUNT
  output_base_dir: "$OUT"

steps:
  PPIFlowStep: false
  MPNNStep_stage1: false
  AbMPNNStep_stage1: true
  FlowpackerStep_stage1: true
  AF3scoreStep_stage1: false
  FilterStep_stage1: false

  RosettaFixStep: false
  PartialStep: false
  MPNNStep_stage2: false
  AbMPNNStep_stage2: false
  FlowpackerStep_stage2: false
  AF3scoreStep_stage2: false
  FilterStep_stage2: false
  ReFoldStep: false
  DockQStep: false
  RosettaRelaxStep: false
  RankStep: false
  ReportStep: false
YAML

cat > "$STEPS" <<YAML
python: "$PYTHON"

PPIFlowStep:
  model_weights: "ckpts/antibody.ckpt"
  config: "$PPI_CFG"

AbMPNNStep_stage1:
  chain_list: "$CHAIN_LIST"
  num_seq_per_target: 1
  sampling_temp: 0.1
  batch_size: 1
  omit_AAs: "C"

FlowpackerStep_stage1:
  config_path: "$FLOW_CFG"
  use_gt_masks: true
YAML

START=$(date +%s)

python "$ROOT/pipeline.py" \
  --task "$TASK" \
  --steps "$STEPS" \
  --stage 1 \
  2>&1 | tee "$ROOT/logs/${NAME}.log"

END=$(date +%s)
ELAPSED=$((END - START))

echo
echo "=========================================="
echo "DONE $TYPE batch $BATCH"
echo "Elapsed seconds: $ELAPSED"
echo "Elapsed hours: $(awk "BEGIN {printf \"%.2f\", $ELAPSED/3600}")"
echo "Flowpacker PDB count:"
find "$OUT/stage1/flowpacker_output/run_1" -name "*.pdb" 2>/dev/null | wc -l
echo "END = $(date)"
echo "=========================================="
