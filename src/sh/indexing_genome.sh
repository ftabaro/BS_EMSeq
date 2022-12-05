#!/bin/bash 

#SBATCH --time=00-05:00:00
#SBATCH --mem=20000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output slurm_%x_%A_%a.out
#SBATCH --error slurm_%x_%A_%a.err

## run Bismark genome preparation step
## genome fasta has already been copied to /scratch/sformich/genome/Bismark folder in order to speed up IO computation
/g/boulard/sformich/conda/envs/smake_EMSeq/bin/bismark_genome_preparation --path_to_aligner /g/boulard/sformich/conda/envs/smake_EMSeq/bin/ /scratch/sformich/genome/Bismark/Mus_musculus
