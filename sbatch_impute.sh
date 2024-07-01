#!/bin/bash
#SBATCH --job-name=IMPUTE               # Job name
#SBATCH --output=%j_output.log          # Output file (with job ID)
#SBATCH --error=%j_error.log            # Error file (with job ID)
#SBATCH --ntasks=1                      # Number of tasks (usually set to 1 for serial jobs)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=16G                       # Memory per node (e.g., 16GB)
#SBATCH --partition=singhlab-gpu        # Partition name

export PATH=$PATH:/hpc/home/jm688/tools/IgImputation

cd /hpc/home/jm688/projects/AbMAP_JM/diseases_heavy
# format_to_fasta.py -i /hpc/home/jm688/projects/AbMAP_JM/diseases_heavy \
#                 -o /hpc/home/jm688/projects/AbMAP_JM/diseases_heavy/FASTA \
#                 -l DEN CMV COVID

# wc -l FASTA/*

for disease in CMV  COVID  DEN
do
echo ${disease}" STARTS"
IgImputation.sh -i FASTA/${disease}.fasta \
                -o IMPUTED/${disease}/${disease}
done