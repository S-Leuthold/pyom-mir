#!/bin/bash
#===============================================================================
# Script:  launch-pyom-analysis.sh
# Purpose: Launch PyOM-MIR analysis on Sybil HPC (nohup background)
# Author:  Sam Leuthold
# Created: 2026-02-04
#===============================================================================

## ---------------------------------------------------------------------------
## Thread controls — pin everything to 1 thread
## ---------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLAS_NUM_THREADS=1
export LAPACK_NUM_THREADS=1

## Allow parallelly to see all cores
export R_PARALLELLY_MAXWORKERS_LOCALHOST=100
export MC_CORES=60

## ---------------------------------------------------------------------------
## Validate environment
## ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

echo "==============================================================================="
echo "PyOM-MIR SPA-C Analysis — Sybil HPC"
echo "==============================================================================="
echo "Start time:        $(date)"
echo "Working directory: $(pwd)"
echo

echo "Checking required files..."
for file in "run-pyom-analysis.R" "config-pyom-analysis.yml" "data/horizons_data.rds"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file NOT FOUND"
        exit 1
    fi
done
echo

## ---------------------------------------------------------------------------
## Launch
## ---------------------------------------------------------------------------

LOG_FILE="pyom_$(date +%Y%m%d_%H%M%S).log"

echo "Launching production run..."
echo "Log file: ${LOG_FILE}"
echo "==============================================================================="

nohup Rscript run-pyom-analysis.R > "${LOG_FILE}" 2>&1 &
PID=$!

echo
echo "Started with PID: ${PID}"
echo "${PID}" > .pyom_pid
echo
echo "Monitor:  tail -f ${LOG_FILE}"
echo "Status:   ps -p ${PID}"
echo "Kill:     kill ${PID}"
echo
echo "==============================================================================="
echo "Launched at $(date)"
echo "==============================================================================="
