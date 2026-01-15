#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH -J HR_include
#SBATCH --output=log.tag.include_MWest.oe
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

python3 write_includes.py ~/ZoomIns MWest -1
