#!/usr/bin/env bash

# Scanpy current-batch paths derived from the repository-wide configuration.

SCANPY_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCLC_PROJECT_CONFIG="${SCLC_PROJECT_CONFIG:-$(cd "${SCANPY_SCRIPT_DIR}/../../.." && pwd)/project_config.sh}"
[[ -s "$SCLC_PROJECT_CONFIG" ]] || {
    echo "ERROR: project configuration missing: $SCLC_PROJECT_CONFIG" >&2
    return 1 2>/dev/null || exit 1
}
# shellcheck disable=SC1090
source "$SCLC_PROJECT_CONFIG"

: "${SCANPY_PYTHON:=${SCLC_CONDA_ROOT}/envs/scanpy310/bin/python}"
export SCANPY_SCRIPT_DIR SCANPY_PYTHON
