#!/bin/bash
#SBATCH --partition=Nebula
#SBATCH --job-name=default
#SBATCH --output=%j.out
#SBATCH --error=%j.err 
#SBATCH --ntasks=1
#SBATCH --time=40:00:00
#SBATCH --mem=8GB

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# Code starts here -----------------------------------------


python knn.py test.npy


cd $SLURM_SUBMIT_DIR
echo "Hello World! I ran on compute node $(/bin/hostname -s)"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
