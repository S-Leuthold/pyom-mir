#!/bin/bash
#===============================================================================
# Script:  transfer-to-sybil.sh
# Purpose: Transfer PyOM analysis files to Sybil HPC scratch
# Author:  Sam Leuthold
# Created: 2026-02-04
#===============================================================================

SYBIL_USER="samleuth"
SYBIL_HOST="sybil.nrel.colostate.edu"
SYBIL_BASE="/scratch/samleuth/pyom-mir"

## ---------------------------------------------------------------------------
## Resolve paths
## ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "==============================================================================="
echo "Transfer PyOM-MIR files to Sybil"
echo "==============================================================================="
echo "Local project: ${PROJECT_DIR}"
echo "Target:        ${SYBIL_USER}@${SYBIL_HOST}:${SYBIL_BASE}"
echo

## ---------------------------------------------------------------------------
## Validate local files
## ---------------------------------------------------------------------------

echo "Checking local files..."

SCRIPTS=(
    "${SCRIPT_DIR}/run-pyom-analysis.R"
    "${SCRIPT_DIR}/launch-pyom-analysis.sh"
    "${SCRIPT_DIR}/config-pyom-analysis.yml"
)

DATA_FILE="${PROJECT_DIR}/data/processed/horizons_data.rds"

ALL_OK=true

for file in "${SCRIPTS[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $(basename "$file")"
    else
        echo "  ✗ $(basename "$file") NOT FOUND"
        ALL_OK=false
    fi
done

if [ -f "$DATA_FILE" ]; then
    SIZE=$(du -h "$DATA_FILE" | cut -f1)
    echo "  ✓ horizons_data.rds (${SIZE})"
else
    echo "  ✗ horizons_data.rds NOT FOUND at ${DATA_FILE}"
    ALL_OK=false
fi

if [ "$ALL_OK" = false ]; then
    echo
    echo "Missing files — aborting."
    exit 1
fi
echo

## ---------------------------------------------------------------------------
## Create remote directory structure
## ---------------------------------------------------------------------------

echo "Creating remote directories..."
ssh "${SYBIL_USER}@${SYBIL_HOST}" "mkdir -p ${SYBIL_BASE}/{data,results}"
echo "  ✓ ${SYBIL_BASE}/"
echo "  ✓ ${SYBIL_BASE}/data/"
echo "  ✓ ${SYBIL_BASE}/results/"
echo

## ---------------------------------------------------------------------------
## Transfer files
## ---------------------------------------------------------------------------

echo "Transferring scripts..."
scp "${SCRIPTS[@]}" "${SYBIL_USER}@${SYBIL_HOST}:${SYBIL_BASE}/"

echo "Transferring data..."
scp "${DATA_FILE}" "${SYBIL_USER}@${SYBIL_HOST}:${SYBIL_BASE}/data/"

if [ $? -eq 0 ]; then
    echo
    echo "✓ All files transferred!"
    echo
    echo "Next steps:"
    echo "  1. ssh ${SYBIL_USER}@${SYBIL_HOST}"
    echo "  2. cd ${SYBIL_BASE}"
    echo "  3. chmod +x launch-pyom-analysis.sh"
    echo "  4. ./launch-pyom-analysis.sh"
    echo
    echo "==============================================================================="
else
    echo
    echo "✗ Transfer failed!"
    exit 1
fi
