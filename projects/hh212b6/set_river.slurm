#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=10:00:00
#SBATCH --output=set_river.out
#SBATCH --partition=standard
#SBATCH --account=StarFormationTheory

python setup_river2.py
