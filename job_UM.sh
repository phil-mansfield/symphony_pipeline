#!/bin/sh
#SBATCH --array=0-48
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J Cluster_UM
#SBATCH --output=logs/Cluster/log.tag.Cluster_um_%A_%a.oe
#SBATCH --mem-per-cpu=32G
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

#python3 write_um_file.py configs/SymphonyMilkyWay/config.txt ${SLURM_ARRAY_TASK_ID} &&
#    echo "done"
python3 print_UM.py configs/Cluster/config.txt ${SLURM_ARRAY_TASK_ID} &&
   python3 write_um_file.py configs/Cluster/config.txt ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
