#!/bin/bash -l
#PBS -N research_shaom
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
time python run_run_wrapper.py