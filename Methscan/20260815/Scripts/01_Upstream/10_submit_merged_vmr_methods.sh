#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$HERE/../../../.." && pwd)
LOG_ROOT="$PROJECT_ROOT/Methscan/20260815/Logs/merged_vmr_methods"
CPU=${METHSCAN_MERGED_CPU:-60}
MEM=${METHSCAN_MERGED_MEM:-184320MB}
mkdir -p "$LOG_ROOT"

for method in scrublet doubletfinder; do
    dsub \
        -n "methscan_merged_${method}" \
        -R "cpu=${CPU};mem=${MEM}" \
        --cwd "$PROJECT_ROOT" \
        -oo "$LOG_ROOT/methscan_merged_${method}.%J.out" \
        -eo "$LOG_ROOT/methscan_merged_${method}.%J.err" \
        bash Methscan/20260815/Scripts/01_Upstream/09_run_merged_vmr_pipeline.sh \
            "$method" all "$CPU"
done

echo "Both joint 10-sample Methscan VMR jobs submitted. Query them with: djob"
