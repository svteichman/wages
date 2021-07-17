#!/bin/bash
#SBATCH --job-name worker_firm_var      # Set a name for your job.
#SBATCH --partition largemem     # Slurm partition to use
#SBATCH --ntasks 4          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 3-12:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=15000     # Memory limit for each tasks (in MB)
#SBATCH -o outfile_%j.out    # File to which STDOUT will be written
#SBATCH -e errfile_%j.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=whj666@uw.edu # Email to which notifications will be sent
module load R
Rscript run_model7-13-21.R 
