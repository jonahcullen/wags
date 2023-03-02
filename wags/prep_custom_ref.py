#!/usr/bin/env python3

import os
import yaml
import gzip
import shutil
import textwrap
import argparse
import fileinput
from pathlib import Path

def main():
    global profile
    # prepare outdir
    ref_home = os.path.join(outdir,species,ref)
    ref_res  = os.path.join(ref_home,"resources")
    ref_fa   = os.path.join(ref_home,os.path.basename(fasta))

    # generate custom reference and resource directories
    os.makedirs(ref_home,exist_ok=True)
    
    # create jobs dir if not exist
    jobs = os.path.join(ref_home,f"{profile}_logs")
    os.makedirs(jobs,exist_ok=True)
    
    # custom known sites
    d_sites = {}
    if sites:
        os.makedirs(ref_res,exist_ok=True)
        with open(os.path.realpath(os.path.expanduser(sites)),'r') as f:
            for line in f:
                name,vcf = line.strip().split(',')
                assert os.path.isfile(vcf), f"{name} was not found, check path is correct"
                # copy known site to outdir (default ~/.wags/SPECIES/REF/resources)
                vcf_out = os.path.join(ref_res,os.path.basename(vcf))
                shutil.copy(vcf,vcf_out)
                # add updated location to dictionary
                d_sites[name] = vcf_out
   
    # copy reference fasta to outdir
    if not os.path.isfile(ref_fa):
        shutil.copy(fasta,ref_home)
    
    # copy and profile
    prep_dir = os.path.dirname(os.path.realpath(__file__))
    profile_dir = os.path.join(
        str(Path(prep_dir).parents[0]),
        f"profiles/{profile}/{profile}.go_wags"
    )
    
    shutil.copytree(profile_dir,os.path.join(ref_home,f"{profile}.go_wags"),dirs_exist_ok=True)
    
    job_sub = os.path.join(ref_home,f"{profile}.go_wags/{profile}-submit.py")
    with fileinput.FileInput(job_sub,inplace=True,backup=".bak") as file:
        for line in file:
            line = line.replace("DUMMY_PAR",partition)
            line = line.replace("DUMMY_ACC",account)
            print(line,end='')
    
    # custom config to modify
    config = os.path.join(
        str(Path(prep_dir).parents[0]),
        "pipelines/prep_ref/configs/custom_config.yaml"
    )

    # modify config file
    with open(config) as f:
        doc = yaml.safe_load(f)
    # update sif and other cli args
    doc['sif']      = sif
    doc['ref']      = ref
    doc['species']  = species
    doc['ref_dict'] = f"{os.path.splitext(ref_fa)[0]}.dict"
    doc['ref_fasta'] = ref_fa
    # add known site resources if provided
    doc['known_sites'] = d_sites

    # dump
    with open(os.path.join(ref_home,f"{ref}_config.yaml"),'w') as out:
        yaml.dump(doc,out,sort_keys=False)

    # copy snake
    smk = os.path.join(
        str(Path(prep_dir).parents[0]),
        "pipelines/prep_ref/prep_wags.smk",
    )
    shutil.copy(smk,ref_home)

    # submission destination
    job_name = "prep_wags"
    submiss = os.path.join(ref_home,f"{species}_{ref}.{job_name}.{profile}")
    
    # SBATCH directives 
    header = (
        "#!/bin/bash -l\n"
        "#SBATCH -t 6:00:00\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks-per-node=1\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --mem=12gb\n"
        "#SBATCH --mail-type=ALL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --job-name {species}_{ref}.{job_name}.slurm\n"
        f"#SBATCH -o slurm_logs/%j.{species}_{ref}.{job_name}.out\n"
        f"#SBATCH -e slurm_logs/%j.{species}_{ref}.{job_name}.err\n"
        f"#SBATCH -A {account}\n"
    )             
 
    # job submission body 
    with open(submiss, "w") as f:
        print(header, file=f)
        print("set -e\n",file=f)
        print(f"conda activate {snake_env}",file=f)
        print("cd $SLURM_SUBMIT_DIR\n",file=f)

        print(
            textwrap.dedent(
                f"""

                REF_DIR={ref_home}

                snakemake -s prep_wags.smk \\
                    --use-singularity \\
                    --singularity-args "-B $PWD,$REF_DIR" \\
                    --profile {profile}.go_wags \\
                    --configfile {ref}_config.yaml \\
                    --keep-going
                """
            ),file=f
        ) 

    print(
        textwrap.dedent(
            f"""
            submit {species}_{ref}.{job_name}.{profile} from {ref_home} 
            to prep {ref} for wags
        """)
    )

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="wags (custom ref)",
        add_help=False,
        description=(
            "      ___           ___           ___           ___                    \n"
            "     /__/\         /  /\         /  /\         /  /\ *                 \n"
            "    _\_ \:\       /  /::\       /  /:/_       /  /:/_                  \n"
            "   /__/\ \:\     /  /:/\:\     /  /:/ /\     /  /:/ /\                 \n"
            "  _\_ \:\ \:\   /  /:/~/::\   /  /:/_/::\   /  /:/ /::\                \n"
            " /__/\ \:\ \:\ /__/:/ /:/\:\ /__/:/__\/\:\ /__/:/ /:/\:\               \n"
            " \  \:\ \:\/:/ \  \:\/:/__\/ \  \:\ /~~/:/ \  \:\/:/~/:/               \n"
            "  \  \:\ \::/   \  \::/       \  \:\  /:/   \  \::/ /:/                \n"
            "   \  \:\/:/     \  \:\        \  \:\/:/     \__\/ /:/                 \n"
            "    \  \::/       \  \:\        \  \::/        /__/:/   *prepare custom\n"
            "     \__\/         \__\/         \__\/         \__\/     reference     \n\n"
            "Prepare custom reference input for generating GVCFs from FASTQs by \n"
            "providing the reference FASTA.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-r", "--ref",
       #nargs="?",
        metavar="\b",
        help="reference name (e.g. equcab3, canfam4)",
    )
    required.add_argument(
        "-n", "--species",
       #nargs="?",
        metavar="\b",
        help="species name (e.g. horse, dog)",
    )
    required.add_argument(
        "-f", "--fasta",
       #nargs="?",
        metavar="\b",
        help="path to reference fasta"
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
    optional.add_argument(
        "-o", "--out",
        default=os.path.join(os.path.expanduser("~"),".wags/"),
        metavar="",
        help="path to custom reference out dir [default: ~/.wags]"
    )
    optional.add_argument(
        "--sites",
        help=textwrap.dedent('''\
            comma-separated file containing names (col 1) and
            paths to resource VCFs (and indices) (col 2) to be 
            used with --ref custom and --bqsr
        ''')
    )
    optional.add_argument(
        "--profile",
        default="slurm",
        help="HPC job scheduler [default: slurm]",
    )
    optional.add_argument(
        "--sif",
        default=os.path.join(os.path.expanduser("~"),".sif/wags.sif"),
        help="location of container image [default: ~/.sif/wags.sif]"
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )
    
    args      = parser.parse_args()
    ref       = args.ref
    species   = args.species
    fasta     = os.path.realpath(os.path.expanduser(args.fasta))
    snake_env = args.snake_env
    partition = args.partition
    email     = args.email
    account   = args.account
    outdir    = args.out
    sites     = args.sites
    profile   = args.profile
    sif       = args.sif

    # assert fasta exists
    assert os.path.isfile(fasta), "fasta not found, check path is correct"

    # gunzip if necessary
    if fasta.endswith(".gz"):
        uncomp_fasta = os.path.join(outdir,os.path.splitext(fasta)[0])
        if not os.path.isfile(uncomp_fasta):
            with gzip.open(fasta,'r') as f_in, open(uncomp_fasta,'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        fasta = uncomp_fasta 

    # if non default outdir
    if outdir != os.path.join(os.path.expanduser("~"),".wags/"):
        outdir = os.path.realpath(os.path.expanduser(outdir))
 
    # if non default sif location
    if sif != os.path.join(os.path.expanduser("~"),".sif/wags.sif"):
        sif = os.path.realpath(os.path.expanduser(sif))
    
    # confirm image exitst
    if os.path.isfile(sif):
        print("wags image found!")
    else:
        sif_dir = os.path.join(os.path.expanduser("~"),".sif")
        print(
            f"wags image not found at {os.path.abspath(sif)} -> downloading to {sif_dir}"
        )
        url = "https://s3.msi.umn.edu/wags/wags.sif"
        os.makedirs(sif_dir,exist_ok=True)
        wget.download(url,sif_dir)
        sys.exit("\nrerun without --sif or point to directory containing wags.sif")
    
    main()
