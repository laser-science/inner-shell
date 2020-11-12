#!/bin/bash
#SBATCH --job-name=Lithium
#SBATCH --ntasks=1
#SBATCH --account=phys
#SBATCH --time=8:00:00
#SBATCH --partition=phys

./Li.out
