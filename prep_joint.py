#!/usr/bin/env python3

import os
import sys
import csv
import gzip
import wget
import yaml
import glob
import shutil
import string
import pathlib
import textwrap
import argparse
import fileinput
import pandas as pd
from datetime import datetime
from collections import defaultdict


refs = [
    "canfam3","canfam4","UU_Cfam_GSD_1.0_ROSY",
    "goldenPath",
    "tiger"
]

# NOTE - should also have checks that each sample and gvcf mapping are unique
def validate_mapping(f):
    ''' quick check that the sample-gvcf mapping contains two tab-delimited columns'''
    with open(f, 'r') as infile:
        for line in infile:
            cols = line.split('\t')
            col_length = len(cols)
            if col_length == 2:
                continue
            else:
                sys.exit(textwrap.dedent('''\n
                    ill-formatted input (-g/--gvcf) - ensure each row contains
                    two columns, sample followed by path to gvcf, separated
                    by a tab
                '''))

def main():

    # load known config
    try:
        with open(config) as f:
            doc = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"{config} does not exist - ensure correct path")
        
    # update config with mapping cohort, and temp dirs    
    doc['joint_cohort'] = os.path.basename(gvcfs)
    doc['date']        = datetime.today().strftime('%Y%m%d')
    if os.path.isfile(gvcfs):
        if validate_mapping(gvcfs):
            doc['joint_cohort'] = gvcfs
    else:
        sys.exit("mapping file (--gvcfs) does not exist - ensure correct path")
    doc['tmp_dir']['sites_only_gather_vcf'] = os.path.join(outdir, '.sites_gather')
    doc['tmp_dir']['unfilt_gather_vcf']     = os.path.join(outdir, '.unfilt_gather')

    # get values from config
    bucket    = doc['bucket']
    alias     = doc['alias']
    profile   = doc['profile']
    ref       = doc['ref']
    ref_fasta = doc['ref_fasta']
    ref_dict  = doc['ref_dict']
    ref_gtf   = doc['ref_gtf']

    # dump modified config to outdir
    os.makedirs(outdir, exist_ok=True)
    config_out = os.path.join(outdir, os.path.basename(config))
    with open(config_out,'w') as f:
        yaml.dump(doc, f, sort_keys=False)

    # prepare outdir with pipeline inputs and profile
    # create jobs dir if not exist
    jobs = os.path.join(outdir, f"{profile}_logs")
    if not os.path.exists(jobs):
        os.makedirs(jobs)

    # input templates
    pipeline  = "joint_genotyping"
    rules     = "rules"
    snake_n   = "many_wags.smk"
    profile_n = f"{profile}.go_wags"

    # copy snakefile, rules, config, and profile to working dir
    smk = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "pipelines",
        pipeline,
        f"inputs/{remote}/{recal}",
        snake_n
    )
    
    rules = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "pipelines",
        pipeline,
        f"inputs/{remote}/{recal}",
        rules
    )
    
    profile_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        f"profiles/{profile}",
        profile_n
    )
    
    src = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "pipelines",
        pipeline,
        "src"
    )
    
    input_names = [snake_n,"rules",profile_n,"src"]
    dst_files = [os.path.join(outdir,i) for i in input_names]
    
    for i in zip([smk,rules,profile_dir,src],dst_files):
        if not os.path.exists(i[1]):
            if os.path.isfile(i[0]):
                shutil.copy(i[0],i[1])
            else:
                shutil.copytree(i[0],i[1])
                # modify profile profile-submit.py for user-supplied parition
                if f"{profile}.go_wags" in i[0]:
                    job_sub = os.path.join(i[1],f"{profile}-submit.py")
                    with fileinput.FileInput(job_sub,inplace=True,backup=".bak") as file:
                        for line in file:
                            line = line.replace("DUMMY_PAR",partition)
                            line = line.replace("DUMMY_ACC",account)
                            print(line,end='')

    # copy gvcfs to output
    dst_gvcf = os.path.join(outdir, os.path.basename(gvcfs))
    if not os.path.isfile(dst_gvcf):
        shutil.copy(gvcfs, dst_gvcf)

    # submission destination
    job_name = snake_n.split('.')[0]
    gvcf_base = os.path.splitext(os.path.basename(gvcfs))[0]
    submiss = os.path.join(outdir, f"{gvcf_base}_{ref}.{job_name}.{profile}")
     
    # sbatch directives 
    header = (
        "#!/bin/bash -l\n"
        "#SBATCH -t 72:00:00\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks-per-node=1\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --mem=12gb\n"
        "#SBATCH --mail-type=ALL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --job-name {gvcf_base}_{ref}.{job_name}\n"
        f"#SBATCH -o slurm_logs/%j.{gvcf_base}_{ref}.{job_name}.out\n"
        f"#SBATCH -e slurm_logs/%j.{gvcf_base}_{ref}.{job_name}.err\n"
        f"#SBATCH -A {account}\n"
        f"#SBATCH -p {partition}\n"
    )             

    # job submission body 
    with open(submiss, "w") as f:
        print(header, file=f)
        print("set -e\n",file=f)
        print(f"conda activate {snake_env}",file=f)
        print("cd $SLURM_SUBMIT_DIR\n",file=f)

        if ref not in refs:
            print(f"FASTA_DIR={ref_fasta}\n", end="", file=f)
            print(f"DICT_DIR={ref_dict}\n", end="", file=f)
            print(f"GTF_DIR={ref_gtf}\n", end="", file=f)

            print(
                textwrap.dedent(
                    f"""
                    snakemake -s {snake_n} \\
                        --use-singularity \\
                        --singularity-args "-B $PWD,$FASTA_DIR,$DICT_DIR,$GTF_DIR" \\
                        --profile {profile_n} \\
                        --configfile {os.path.basename(config)} \\
                        --keep-going
                    """
                ),file=f
            ) 
        else:     # NOTE - THE GVCF DIR IS HARDCODED BELOW AND NEEDS TO BE FIXED BY CHECKING THE GVCF FILE       
            print(
                textwrap.dedent(
                    f"""
                    GVCF_DIR=/scratch.global/friedlab_JOINT_EVAL/fetched_gvcfs/gvcfs/
                    snakemake -s {snake_n} \\
                        --use-singularity \\
                        --singularity-args "-B $PWD,$GVCF_DIR" \\
                        --profile {profile_n} \\
                        --configfile {os.path.basename(config)} \\
                        --keep-going
                    """
                ),file=f
            ) 
            
        if "s3" in remote:
            print(f"# save {profile} err/out logs",end="",file=f)
            print(
                textwrap.dedent(
                    f"""
                    mc cp ./{profile}_logs/* \\
                        {alias}/{bucket}/wgs/pipeline/{ref}/{profile}_logs/
                    """
                ),file=f
            )
    
    num_gvcfs = sum(1 for line in open(gvcfs))
    print(f"many wags setup to process {num_gvcfs} samples at {outdir}")
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="wags",
        add_help=False,
        description=(
            "       ___           ___           ___           ___      \n"
            " MANY /__/\         /  /\         /  /\         /  /\     \n"
            "     _\_ \:\       /  /::\       /  /:/_       /  /:/_    \n"
            "    /__/\ \:\     /  /:/\:\     /  /:/ /\     /  /:/ /\   \n"
            "   _\_ \:\ \:\   /  /:/~/::\   /  /:/_/::\   /  /:/ /::\  \n"
            "  /__/\ \:\ \:\ /__/:/ /:/\:\ /__/:/__\/\:\ /__/:/ /:/\:\ \n"
            "  \  \:\ \:\/:/ \  \:\/:/__\/ \  \:\ /~~/:/ \  \:\/:/~/:/ \n"
            "   \  \:\ \::/   \  \::/       \  \:\  /:/   \  \::/ /:/  \n"
            "    \  \:\/:/     \  \:\        \  \:\/:/     \__\/ /:/   \n"
            "     \  \::/       \  \:\        \  \::/        /__/:/    \n"
            "      \__\/         \__\/         \__\/         \__\/   \n\n"
            "*many* wags generates all required input to joint call GVCFs \n"
            "into a single VCF following GATK best practices."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-c", "--config",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to config file generated by config_joint.py"
    )
    required.add_argument(
        "-g", "--gvcfs",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent(f'''\
            path to file with mapping of sample name to
            gvcf in tab delimited format with no header
        ''')
    )
    required.add_argument(
        "-s", "--snake-env",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="conda environment with snakemake"
    )
    required.add_argument(
        "-p", "--partition",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="default partition(s) to use (e.g. 'par1' or 'par1,par2'"
    )
    required.add_argument(
        "-e", "--email",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="email address for job logs"
    )
    required.add_argument(
        "-a", "--account",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="default scheduler account"
    )
    required.add_argument(
        "-o", "--out",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to out dir"
    )
    optional.add_argument(
        "--vqsr-snps",
        action='store_true',
        help=textwrap.dedent('''\
            recalibrate snps with known sites.
            if false, hard filtering will be applied
            to all snps [default: false]
        ''')
    )
    optional.add_argument(
        "--vqsr-indels",
        action='store_true',
        help=textwrap.dedent('''\
            recalibrate indels with known sites.
            if false, hard filtering will be applied
            to all non-snps [default: false]
        ''')
    )
    optional.add_argument(
        "--profile",
        default="slurm",
        help="HPC job scheduler [default: slurm]",
    )
    optional.add_argument(
        "--remote",
        default="local",
        type=str.lower,
        help="save outputs to remote: S3, SFTP [default: local]",
    )
    optional.add_argument(
        "--sif",
        default=os.path.join(os.path.expanduser("~"),".sif/wags.sif"),
        help="path to singularity image file [default: ~/.sif/wags.sif]"
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )

    args = parser.parse_args()
    config     = args.config
    gvcfs      = args.gvcfs
    snps       = args.vqsr_snps
    indels     = args.vqsr_indels
    outdir     = os.path.realpath(os.path.expanduser(args.out))
    snake_env  = args.snake_env
    partition  = args.partition
    email      = args.email
    account    = args.account
    sif        = args.sif
    profile    = args.profile
    remote     = args.remote.lower()

    # based on snps/indels, define pipeline
    recal = "vqsr_none"
    if snps and indels:
        recal = "vqsr_all"
    elif snps and not indels:
        recal = "vqsr_snps"
    
    # if non default sif location
    if sif != os.path.join(os.path.expanduser("~"),".sif/wags.sif"):
        sif = os.path.realpath(os.path.expanduser(sif))
       #if "~" in sif:
       #    sif = os.path.expanduser(sif)
       #else:
       #    sif = os.path.abspath(sif)
    
    # confirm image exitst
    if os.path.isfile(sif):
        print("wags image found!")
    else:
        sif_dir = os.path.join(os.path.expanduser("~"),".sif")
        print(
            f"wags image not found at {os.path.realpath(sif)} -> downloading to {sif_dir}"
        )
        url = "https://s3.msi.umn.edu/wags/wags.sif"
        os.makedirs(sif_dir,exist_ok=True)
        wget.download(url,sif_dir)
        sys.exit("\nrerun without --sif or point to directory containing wags.sif")

    main()
        
        
