#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$HERE/../../.." && pwd)
LOG_ROOT="$PROJECT_ROOT/MethylVI/20260816/Logs/methscan_vmr_methods"
CPU=${MVI_VMR_DSUB_CPU:-60}
MEM=${MVI_VMR_DSUB_MEM:-184320MB}
mkdir -p "$LOG_ROOT"

for method in scrublet doubletfinder; do
    dsub \
        -n "methylvi_vmr_${method}" \
        -R "cpu=${CPU};mem=${MEM}" \
        --cwd "$PROJECT_ROOT" \
        -oo "$LOG_ROOT/methylvi_vmr_${method}.%J.out" \
        -eo "$LOG_ROOT/methylvi_vmr_${method}.%J.err" \
        env MVI_THREADS="$CPU" MVI_MEMORY_GB=180 MVI_ACCELERATOR=cpu \
        bash MethylVI/20260816/Scripts/18_run_methscan_vmr_method.sh \
            "$method" full
done

echo "Both joint-VMR MethylVI jobs submitted. Query them with: djob"
