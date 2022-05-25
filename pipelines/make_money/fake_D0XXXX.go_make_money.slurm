#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cull0084@umn.edu
#SBATCH --job-name QH_D0XXXX.one_wags.slurm
#SBATCH -o slurm_logs/%j.QH_D0XXXX.one_wags.out
#SBATCH -e slurm_logs/%j.QH_D0XXXX.one_wags.err
#SBATCH -A fried255

set -e

conda activate snake73
cd $SLURM_SUBMIT_DIR


FQ_DIR=/home/fried255/fried255/temp/lclark
PROC_DIR=/scratch.global/friedlab_MONEY_TEST
POP_DIR=/panfs/jay/groups/0/fried255/fried255/working/pipeline/canfam4/20220401

# extract reference dict from container
singularity exec --bind $PWD /home/fried255/cull0084/.sif/wags.sif \
    cp /home/refgen/dog/canfam4/canFam4.dict $PWD


snakemake -s go_make_money.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$FQ_DIR,$PROC_DIR,$POP_DIR" \
    --profile slurm.go_money \
    --configfile canfam4_money.yaml \
    --keep-going

