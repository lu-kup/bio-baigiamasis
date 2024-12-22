#!/bin/bash

#SBATCH -p main
#SBATCH -n2

set -x
date
srun python3 run_models.py
date
