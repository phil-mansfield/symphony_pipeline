#!/bin/sh
#SBATCH --array=1,36
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J EDEN_tag
#SBATCH --output=logs/EDEN/log.tag.EDEN_MilkyWay_8K_%A_%a.oe
#SBATCH --mem-per-cpu=32G
#SBATCH -p kipac

config=configs/EDEN_MilkyWay_8K/config.txt
go run scale_factor_table.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run tag_particles.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run xv.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
