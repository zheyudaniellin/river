#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=10:00:00
#SBATCH --output=set_river.out
#SBATCH --partition=standard
#SBATCH --account=StarFormationTheory

python setup_river3.py

# 9 hours for nwalker=8, steps=300, thread=16
# 11.5 hours for nwalker=6, steps=400, thread=16
