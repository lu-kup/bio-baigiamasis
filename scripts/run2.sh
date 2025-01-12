#!/bin/bash

#SBATCH -p main
#SBATCH -n1

set -x
date
srun python3 run_kmeans_expression.py TT_S1
date
