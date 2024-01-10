#!/bin/bash
#SBATCH --job-name=UM023
#SBATCH --output=UM023.out
#SBATCH --error=UM023.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -p kipac

PERL5LIB=/sdf/home/y/ycwang/UMmw/universemachine/perl
export PERL5LIB
LD_LIBRARY_PATH=/sdf/home/y/ycwang/GSL/lib
export LD_LIBRARY_PATH

mkdir groupcat
/sdf/home/y/ycwang/UMmw/universemachine/split_halo_trees_phase1 hlists/hlist_1.00000.list 125 1 1 1 trees/forests.list > boxlist.txt
perl /sdf/home/y/ycwang/UMmw/universemachine/scripts/parallel_split_halo_trees.pl 125 0.286 0.7 outputs/scales.txt /sdf/group/kipac/u/ycwang/MWmass_new/join_halos/splines boxlist.txt 2 trees/tree_0_0_0.dat 

cd /sdf/home/y/ycwang/UMmw/universemachine/scripts
LD_LIBRARY_PATH=/sdf/home/y/ycwang/GSL/lib
export LD_LIBRARY_PATH
perl make_sf_catalog.pl zoom_cfg/Halo_023.cfg bestfit.dat 1 1
cd ../
./print_sm_catalog /sdf/group/kipac/u/ycwang/MWmass_new/Halo023/output/rockstar/groupcat/sfr_catalog_1.000000.bin scripts/zoom_cfg/Halo_023.cfg
cd /sdf/group/kipac/u/ycwang/MWmass_new/Halo023/output/rockstar/groupcat
mv sfr_catalog_1.000000.bin.txt sfr_catalog_1.000000.txt
