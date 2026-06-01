#!/bin/bash

# Run GADGET4 simulation pipeline with fixed potential.
# This script runs the full pipeline for a GADGET4 simulation with a fixed potential:
# 1) Initial conditions generation
# 2) GADGET parameter file (param.txt) update with specific values for the simulation
# 3) GADGET4 execution with MPI
#
# Usage: ./run_pipeline.sh <run_name> <mass_ratio> <num_particles>

# Exit immediately if:
#  - a command fails (-e)
#  - an undefined variable is used (-u)
#  - any command in a pipeline fails (-o pipefail)
set -euo pipefail

if [ $# -lt 3 ]; then
    echo "❌ Usage: ./run_pipeline.sh <run_name> <mass_ratio> <num_particles>"
    exit 1
fi

RUN_NAME="$1"
MASS_RATIO="$2"
NUM_PARTICLES="$3"

if ! [[ "$NUM_PARTICLES" =~ ^[0-9]+$ ]] || [ "$NUM_PARTICLES" -lt 2 ]; then
    echo "❌ num_particles must be an integer >= 2 (got: $NUM_PARTICLES)"
    exit 1
fi

rm -rf ./output

python3 initial_conditions.py "$MASS_RATIO" "$NUM_PARTICLES"


python3 - "$MASS_RATIO" "param.txt" <<'PY'
import sys, os, shutil, time

# Read command-line arguments
mr = float(sys.argv[1])
param_path = sys.argv[2]

# Parameters to overwrite in param.txt
params = {
    'TimeMax': 5e-9,
    'TimeBetSnapshot': 1e-12,
    'TimeBetStatistics': 1,
    'ErrTolIntAccuracy': 0.001,
    'MaxSizeTimestep': 1e-12,
    'SofteningComovingClass0': 1e-13,
    'SofteningMaxPhysClass0': 1e-13,
}

# Create backup directory and save a timestamped copy of param.txt
bak_dir = 'backups'
os.makedirs(bak_dir, exist_ok=True)
bak = os.path.join(bak_dir, f'param.txt.bak.{time.time():.0f}')
shutil.copyfile(param_path, bak)

# Read parameter file
with open(param_path) as f:
    lines = f.readlines()

# Replace relevant parameters while preserving all other lines
new_lines = []
for line in lines:
    stripped = line.lstrip()
    for k, v in params.items():
        if stripped.startswith(k):
            line = f"{k.ljust(24)}{v:.3e}\n"
            break
    new_lines.append(line)

# Write updated parameters back to param.txt
with open(param_path, 'w') as f:
    f.writelines(new_lines)
PY

nice -n 19 mpirun -np 32 --use-hwthread-cpus ./Gadget4 param.txt

exit 0
