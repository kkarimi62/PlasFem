#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/PlasFem/src

module load openmpi/2.1.3-gnu
module load mpich/3.2.1-gnu
module load gcc/5.3.0
 					    export OMP_NUM_THREADS=8
$EXEC_DIR/fem_run -var OUT_PATH /scratch/${SLURM_JOB_ID} < in.txt
