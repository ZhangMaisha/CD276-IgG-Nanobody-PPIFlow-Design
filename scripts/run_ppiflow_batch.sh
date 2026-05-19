#!/usr/bin/env bash
set -euo pipefail

############################################
# Usage:
# ./run_ppiflow_batch.sh antibody 1 25
# ./run_ppiflow_batch.sh nanobody 3 25
############################################

TYPE=${1:?Need type: antibody or nanobody}
BATCH=${2:?Need batch number}
N=${3:-25}

ROOT="/teams/YingChiLab_1702378116/MaishaZhang/PPIFlow-main"
PYTHON="/root/miniconda3/envs/ppiflow_af3/bin/python"
CFG="$ROOT/tool/PPIFlow/configs/inference_nanobody.yaml"

mkdir -p "$ROOT/demo_cd276/batches"
mkdir -p "$ROOT/output/raw_designs"
mkdir -p "$ROOT/logs"

############################################
# TYPE SETTINGS
############################################

if [ "$TYPE" = "antibody" ]; then

    NAME="CD276_IgG_batch${BATCH}"

    OUT="$ROOT/output/raw_designs/IgG_batch_${BATCH}"

    FRAMEWORK="$ROOT/example/target_and_framework_pdb/6nou_scfv_framework.pdb"

    MODEL="ckpts/antibody.ckpt"

    EXTRA_TASK=$(cat <<EOF
  light_chain: "B"
  cdr_length: "CDRH1,8-8,CDRH2,8-8,CDRH3,10-20,CDRL1,6-9,CDRL2,3-3,CDRL3,9-11"
EOF
)

elif [ "$TYPE" = "nanobody" ]; then

    NAME="CD276_NB_batch${BATCH}"

    OUT="$ROOT/output/raw_designs/NB_batch_${BATCH}"

    FRAMEWORK="$ROOT/example/target_and_framework_pdb/7eow_nanobody_framework.pdb"

    MODEL="ckpts/nanobody.ckpt"

    EXTRA_TASK=$(cat <<EOF
  cdr_length: "CDRH1,8-8,CDRH2,8-8,CDRH3,9-21"
EOF
)

else
    echo "TYPE must be antibody or nanobody"
    exit 1
fi

############################################
# CLEAN OUTPUT
############################################

rm -rf "$OUT"

TASK="$ROOT/demo_cd276/batches/${NAME}.yaml"
STEPS="$ROOT/demo_cd276/batches/${NAME}_steps.yaml"

############################################
# TASK YAML
############################################

cat > "$TASK" <<YAML
task:
  name: "$NAME"
  gentype: "$TYPE"
  antigen_pdb: "$ROOT/input/CD276_C.pdb"
  antigen_chain: "C"
  specified_hotspots: "C242,C243,C244,C245,C246,C247,C248,C282,C283,C284,C285"
  framework_pdb: "$FRAMEWORK"
  heavy_chain: "A"
$EXTRA_TASK
  samples_per_target: $N
  output_base_dir: "$OUT"

steps:
  PPIFlowStep: true
  MPNNStep_stage1: false
  AbMPNNStep_stage1: false
  FlowpackerStep_stage1: false
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

############################################
# STEPS YAML
############################################

cat > "$STEPS" <<YAML
python: "$PYTHON"

PPIFlowStep:
  model_weights: "$MODEL"
  config: "$CFG"
YAML

############################################
# RUN
############################################

echo "===================================="
echo "TYPE: $TYPE"
echo "BATCH: $BATCH"
echo "N: $N"
echo "OUTPUT: $OUT"
echo "START:"
date
echo "===================================="

START=$(date +%s)

python "$ROOT/pipeline.py" \
  --task "$TASK" \
  --steps "$STEPS" \
  --stage 1 \
  2>&1 | tee "$ROOT/logs/${NAME}.log"

END=$(date +%s)
ELAPSED=$((END - START))

echo
echo "===================================="
echo "DONE"
echo "Elapsed sec: $ELAPSED"
echo "Elapsed hr : $(awk "BEGIN {printf \"%.2f\", $ELAPSED/3600}")"
echo "PDB count:"
find "$OUT/stage1/ppiflow_output" -name "*.pdb" | wc -l
echo "END:"
date
echo "===================================="
