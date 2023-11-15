#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH -J Group_converge
#SBATCH --output=logs/Group/converge.log
#SBATCH -p kipac

config=configs/Group/config.txt
python3 convergence_snaps.py ${config} -1
