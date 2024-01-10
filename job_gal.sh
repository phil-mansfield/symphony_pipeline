#!/bin/sh
#SBATCH --array=0-4
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J r=0.025_MilkyWayHR_gal
#SBATCH --output=logs/MilkyWay/log.gal.MilkyWayHR_%A_%a.oe
#SBATCH --mem-per-cpu=24G
#SBATCH -p kipac

config=configs/MilkyWayHR/config.txt
python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.025 &&
   echo "done"
