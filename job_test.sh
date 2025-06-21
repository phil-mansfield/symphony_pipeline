#!/bin/bash

#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH -J Group_converge
#SBATCH --output=test.loh
#SBATCH -p kipac

echo ":)"
