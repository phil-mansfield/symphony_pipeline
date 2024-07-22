#!/bin/sh
#SBATCH --array=0-4
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J 8K_tag
#SBATCH --output=logs/EDEN/log.tag.EDEN_16K_%A_%a.oe
#SBATCH --mem-per-cpu=24G
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

config=configs/EDEN_MilkyWay_16K/config.txt
go run scale_factor_table.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run tag_particles.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run xv.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
