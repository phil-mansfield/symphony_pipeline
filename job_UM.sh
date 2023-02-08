#!/bin/sh
#SBATCH --partition=kipac
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J MilkyWay_print_UM
#SBATCH --output=logs/log.MilkyWay.%j.oe
#SBATCH --mem-per-cpu=32G

#python3 print_UM.py configs/MilkyWay/config.txt -1 &&
   python3 write_um_file.py configs/MilkyWay/config.txt -1
