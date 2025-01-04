#!/bin/bash

#SBATCH -p main
#SBATCH -n1

set -x
date
srun python3 run_kprototypes_expression.py
date
