#!/bin/bash -l
#SBATCH -t 60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=cull0084@umn.edu 
#SBATCH --job-name jg_snake
#SBATCH -o ./Jobs/%j.jg_snake.out
#SBATCH -e ./Jobs/%j.jg_snake.err
#SBATCH -p small,amdsmall

set -e

conda activate snake73
cd $SLURM_SUBMIT_DIR


snakemake -s joint_genotyping.smk \
    --use-singularity \
    --profile slurm.go_joint/ \
    --configfile canfam4_config.yaml \
    --keep-going

