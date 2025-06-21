#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=16
#SBATCH -J MWest_include
#SBATCH --output=log.tag.All_include.oe
#SBATCH --mem=32G
#SBATCH --partition roma 
#SBATCH --account kipac:kipac

python3 write_includes.py configs/MilkyWay/config.txt -1 &&
    python3 write_includes.py configs/MWest/config.txt -1 &&
    python3 write_includes.py configs/MilkyWayHR/config.txt -1 &&
   echo "done"
