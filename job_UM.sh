#!/bin/sh
#SBATCH --array=0-48
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J Group_UM
#SBATCH --output=logs/Group/log.tag.Group_um_%A_%a.oe
#SBATCH --mem-per-cpu=32G
#SBATCH -p kipac

#python3 write_um_file.py configs/Group/config.txt ${SLURM_ARRAY_TASK_ID} &&
#    echo "done"
python3 print_UM.py configs/Group/config.txt ${SLURM_ARRAY_TASK_ID} &&
   python3 write_um_file.py configs/Group/config.txt ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
