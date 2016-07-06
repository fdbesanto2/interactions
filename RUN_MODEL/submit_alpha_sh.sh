#!/bin/bash
#PBS -q qfbb
#PBS -N est_ab
#PBS -l nodes=12:ppn=1
#PBS -l walltime=01:00:00
#PBS -j oe

module load bioinformatics/R/3.2.5
cd $PBS_O_WORKDIR
Rscript ./parallel_prg_sh.r
