#!/bin/sh
#SBATCH --array=31-44
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J MilkyWay_gal
#SBATCH --output=logs/MilkyWay/log.gal.MilkyWay_%A_%a.oe
#SBATCH --mem-per-cpu=32G
#SBATCH -p kipac

config=configs/MilkyWay/config.txt
python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} um &&
   echo "done"
