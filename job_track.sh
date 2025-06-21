#!/bin/bash

#SBATCH --array=0-53
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=28G
#SBATCH -J MWest_track
#SBATCH --output=logs/EDEN/log.track.EDEN_8K_%A_%a.oe
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

config=configs/EDEN_MilkyWay_8K/config.txt
suffix=fid4
snap_range=0:235

python3 find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&
   python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} --snap_range=${snap_range} --reset &&
   python3 convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   python3 convergence_snaps.py ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
