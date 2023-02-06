#!/bin/sh
#SBATCH --partition=kipac
#SBATCH --array=4-95
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J Cluster
#SBATCH --output=/sdf/home/t/tdacunha/logs/log.tag.Cluster.%j.oe
#SBATCH --mem-per-cpu=32G

echo "started"
config=~/symphony_pipeline/configs/Cluster/config.txt
go run ~/symphony_pipeline/tag_particles.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run ~/symphony_pipeline/xv.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
