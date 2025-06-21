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
python3 galaxy_properties.py ${config} ${SLURM_ARRAY_TASK_ID} um um_fit r=0.005 r=0.008 r=0.015 r=0.025 r=0.05 r=0.1 r=0.2 r=0.4

date
