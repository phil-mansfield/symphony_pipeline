#!/bin/sh
#SBATCH --partition=kipac
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J Cluster_make_trees
#SBATCH --output=/sdf/home/t/tdacunha/logs/log.Cluster_tree.%j.oe
#SBATCH --mem-per-cpu=32G

config=~/symphony_pipeline/configs/Cluster/config.txt
index=-1
go run ~/symphony_pipeline/write_binary_tree.go ${config} ${index} && go run ~/symphony_pipeline/write_tree_header.go ${config} ${index} && go run ~/symphony_pipeline/write_subhalos_file.go ${config} ${index}
