#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$HERE/../../.." && pwd)
LOG_ROOT="$PROJECT_ROOT/MethylVI/20260816/Logs/methscan_qc_methods"
CPU=${MVI_DSUB_CPU:-60}
MEM=${MVI_DSUB_MEM:-184320MB}
PROFILE=${MVI_PROFILE:-100k}
mkdir -p "$LOG_ROOT"

for method in scrublet doubletfinder; do
  dsub \
    -n "methylvi_${method}_${PROFILE}" \
    -R "cpu=${CPU};mem=${MEM}" \
    --cwd "$PROJECT_ROOT" \
    -oo "$LOG_ROOT/methylvi_${method}_${PROFILE}.%J.out" \
    -eo "$LOG_ROOT/methylvi_${method}_${PROFILE}.%J.err" \
    env MVI_THREADS="$CPU" MVI_MEMORY_GB=180 MVI_ACCELERATOR=cpu \
    bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh \
      "$method" full "$PROFILE"
done

echo "Both Methscan-QC MethylVI jobs submitted. Query them with: djob"
