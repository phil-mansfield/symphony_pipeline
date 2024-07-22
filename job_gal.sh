#!/bin/sh
#SBATCH --array=0-44
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J MilkyWay_gal
#SBATCH --output=logs/MilkyWay/log.gal.MilkyWayHR_%A_%a.oe
#SBATCH --mem-per-cpu=24G
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

date

config=configs/MilkyWay/config.txt
python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} um &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.005 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.008 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.015 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.025 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.05 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.1 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.2 &&
    python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} r=0.4 &&


date
