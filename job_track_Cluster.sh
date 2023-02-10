#!/bin/bash
#SBATCH --array=0-95
#SBATCH --time=0:48:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/sdf/home/t/tdacunha/logs/track_logs/log.tag.Cluster.%j.oe
#SBATCH -p kipac
#SBATCH -J Cluster


config=~/symphony_pipeline/configs/Cluster/config.txt
suffix=fid
snap_range=0:199

python3 ~/symphony_pipeline/find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&
#python3 ~/symphony_pipeline/print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --reset --suffix=${suffix} --snap_range=${snap_range} &&
#python3 ~/symphony_pipeline/print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} --snap_range=${snap_range} &&
   #python3 ~/symphony_pipeline/convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   echo "done"
