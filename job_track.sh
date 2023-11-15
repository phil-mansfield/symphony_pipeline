#!/bin/bash

#SBATCH --array=0-0
#SBATCH --time=2-00:00:0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH -J DMO_track
#SBATCH --output=logs/Cluster/log.track.Cluster_%A_%a.oe
#SBATCH -p kipac

config=configs/Cluster/config.txt
suffix=fid3
snap_range=0:200

python3 find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&
   python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} --snap_range=${snap_range} --reset &&
   python3 convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   python3 convergence_snaps.py ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
