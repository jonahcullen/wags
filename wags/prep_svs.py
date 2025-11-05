#!/usr/bin/env python3

import os
import sys
import csv
import gzip
import yaml
import glob
import shutil
import string
import tempfile
import textwrap
import argparse
import fileinput
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import defaultdict

refs = [
    "canfam3","canfam4","UU_Cfam_GSD_1.0_ROSY",
    "goldenPath","Arabian","Shire","Thoroughbred",
    "tiger",
    "Fca126_mat1.0",
    "alpaca",
    "ARS-UI_Ramb_v3.0",
    "ARS1.2"
]

# get config dir
prep_dir = Path(__file__).resolve().parent.parent
configs_dir = prep_dir / "pipelines" / "many_svs" / "configs"
# put available configs into a dictionary
config_d = {}
for species_dir in configs_dir.iterdir():
    if species_dir.is_dir():
        for config_f in species_dir.iterdir():
            if config_f.suffix == ".yaml":
                species = config_f.stem.replace("_config", "")
                config_d[species] = species_dir.name
                
def gather_crams(cram_name):
    Path(f"{outdir}/sample_lists/").mkdir(parents=True, exist_ok=True)
    cram_file = f"{outdir}/sample_lists/crams.txt"
    crams = 0
    with open(cram_file, 'w') as cf:
        for i in glob.glob(f"{cram_name}/*.cram"):
            file_name = i.split('/')[-1]
            sample_prefix = file_name.split('.')[0]
            cf.write(f"{sample_prefix}\n")
            crams += 1
    print(f"found {crams} cram files") 

def gather_svars(svar_name):
    svars = ['delly', 'gridss', 'manta', 'smoove']
    list_dir = f"{outdir}/sample_lists/"
    for sv in svars:
        svs = 0
        with open(f"{list_dir}/{sv}.txt", 'w') as sf:
            for i in glob.glob(f"{svar_name}/{sv}/*.gz"):
                file_name = i.split('/')[-1]
                sample_prefix = file_name.split('.')[0]
                sf.write(f"{sample_prefix}\n")
                svs += 1
        print(f"found {svs} {sv} files")
 
def main():
    
    # load known config
    prep_path = Path(__file__).resolve()
    config_path = prep_path.parent.parent / f"pipelines/many_svs/configs/{config_d[ref]}/{ref}_config.yaml"

    try:
        with open(config_path) as f:
            doc = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"{config} does not exist - ensure correct path")
        
    # update config with user info, mapping cohort, interval options, and temp dirs
    
    Path(outdir).mkdir(parents=True, exist_ok=True)
    cram_name = Path(crams).resolve()
    if Path(cram_name).exists():
        gather_crams(cram_name)
    else:
        sys.exit("cannot find cram directory - ensure correct path")

    svar_name = Path(svars).resolve()
    if Path(svars).exists():
        gather_svars(svar_name)
    else:
        sys.exit("cannot find svar directory - ensure correct path")

    doc['date'] = datetime.today().strftime('%Y%m%d')
    
    # sample files
    doc['cram_samples']    = f"{outdir}/sample_lists/crams.txt"
    doc['delly_samples']   = f"{outdir}/sample_lists/delly.txt"
    doc['gridss_samples']  = f"{outdir}/sample_lists/gridss.txt"
    doc['manta_samples']   = f"{outdir}/sample_lists/manta.txt"
    doc['smoove_samples']  = f"{outdir}/sample_lists/smoove.txt"
    doc['cram_dir']        = str(cram_name)
    doc['svar_dir']        = str(svar_name)

    #survivor / sv details
    doc['sv_length']  = sv_length 
    doc['sv_callers'] = sv_callers
    doc['sv_overlap'] = sv_overlap
    
    # other
    doc['sif']                = sif
    doc['bucket']             = bucket
    doc['alias']              = alias
    doc['profile']            = profile
    
    # get values from config
    ref_fasta = doc['ref_fasta']
    ref_gtf   = doc['ref_gtf']

    # dump modified config to outdir
    Path(outdir).mkdir(parents=True, exist_ok=True)
    config_out = str(Path(outdir) / f"{ref}_config.yaml")
    with open(config_out,'w') as f:
        yaml.dump(doc, f, sort_keys=False)

    # prepare outdir with pipeline inputs and profile
    jobs = Path(outdir) / f"{profile}_logs"
    jobs.mkdir(parents=True, exist_ok=True)

    # input templates
    pipeline  = "many_svs"
    rules     = "rules"
    snake_n   = "many_svs.smk"
    profile_n = f"{profile}.go_wags"

    # copy snakefile, rules, config, and profile to working dir
    prep_dir = Path(__file__).resolve().parent

    smk = (
        prep_dir.parents[0]
        / "pipelines/many_svs"
        / f"inputs/{remote}"
        / snake_n
    ) 
    
    rules = (
        prep_dir.parents[0]
        / "pipelines/many_svs"
        / f"inputs/{remote}"
        / rules
    )
    
    profile_dir = (
        prep_dir.parents[0]
        / f"profiles/{profile}"
        / profile_n
    )
    
    src = (
        prep_dir.parents[0]
        / "pipelines/many_svs"
        / "src"
    )
    
    input_names = [snake_n, "rules", profile_n, "src"]
    dst_files = [str(Path(outdir) / i) for i in input_names]
    
    for i in zip([smk, rules, profile_dir, src], dst_files):
        if not Path(i[1]).exists():
            if Path(i[0]).is_file():
                shutil.copy(i[0],i[1])
            else:
                shutil.copytree(i[0],i[1])
                # modify profile profile-submit.py for user-supplied parition
                if f"{profile}.go_wags" in str(i[0]):
                    if profile == 'lsf':
                        job_sub = str(Path(i[1]) / f"{profile}_submit.py")
                    else:
                        job_sub = str(Path(i[1]) / f"{profile}-submit.py")
                    with fileinput.FileInput(job_sub, inplace=True, backup=".bak") as file:
                        for line in file:
                            line = line.replace("DUMMY_PAR",partition)
                            line = line.replace("DUMMY_ACC",account)
                            print(line,end='')

    # submission destination
    job_name = snake_n.split('.')[0]
    submiss = os.path.join(outdir, f"svs_{ref}.{job_name}.{profile}")
     
    # sbatch directives 
    default_header = (
        "#!/bin/bash -l\n"
        f"#SBATCH -t {walltime}:00:00\n"
        "#SBATCH --nodes=1\n"
        "#SBATCH --ntasks-per-node=1\n"
        "#SBATCH --cpus-per-task=1\n"
        "#SBATCH --mem=12gb\n"
        "#SBATCH --mail-type=ALL\n"
        f"#SBATCH --mail-user={email}\n"
        f"#SBATCH --job-name svs_{ref}.{job_name}\n"
        f"#SBATCH -o slurm_logs/%j.svs_{ref}.{job_name}.out\n"
        f"#SBATCH -e slurm_logs/%j.svs_{ref}.{job_name}.err\n"
        f"#SBATCH -A {account}\n"
        f"#SBATCH -p {partition}\n"
    )             
    lsf_header = (
        "#!/bin/bash -l\n"
        f"#BSUB -W {walltime}:00\n"
        "#BSUB -n 1\n"
        "#BSUB -R span[hosts=1]\n"
        "#BSUB -R rusage[mem=12GB]\n"
        f"#BSUB -J svs_{ref}.{job_name}\n"
        f"#BSUB -o {profile}_logs/%J.svs_{ref}.{job_name}.out\n"
        f"#BSUB -e {profile}_logs/%J.svs_{ref}.{job_name}.err\n\n"
        "module load apptainer"
        
    )

    # job submission body 

    tmp_dir = tempfile.gettempdir()
    with open(submiss, "w") as f:
        if profile == 'lsf':
            print(lsf_header, file=f)
        else:
            print(default_header, file=f)
        print("set -e\n",file=f)
        print(f"conda activate {snake_env}",file=f)
        if profile != 'lsf':
            print("cd $SLURM_SUBMIT_DIR\n",file=f)

        print(
            textwrap.dedent(
                f"""
                CRAM_DIR={cram_name}
                SVAR_DIR={svar_name}
                TMP_DIR={tmp_dir}
                snakemake -s {snake_n} \\
                    --use-singularity \\
                    --singularity-args "-B $CRAM_DIR,$SVAR_DIR,$TMP_DIR" \\
                    --profile {profile_n} \\
                    --configfile {ref}_config.yaml \\
                    --keep-going \\
                    --rerun-incomplete
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
    
    # list avilable species and configs and exit
    parser.add_argument(
        "--configs",
        action="store_true",
        help="list all available species and assemblies (configs)"
    )

    args, unknown = parser.parse_known_args()

    if args.configs:
        print("available species and assemblies (configs):")
        for k,v in sorted(config_d.items(), key=lambda x: (x[1], x[0])):
            print(f"  {k:20} ({v})")
        sys.exit(0)

    # define required only if --configs was not used
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-r", "--ref",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent(f'''\
            select reference to use: 
                {", ".join(refs)}
            if using custom reference, ensure provided name
            is exact match to name (--ref) used with 
            prep_custom_ref.py
        ''')
    )
    required.add_argument(
        "-c", "--crams",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent(f"""\
            path to cram directory
        """)
    )
    required.add_argument(
        "-sv", "--svars",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent(f"""\
            path to sv directory
        """)
    )
    required.add_argument(
        "-b", "--bucket",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent('''\
            results directory or bucket name if using
            --remote s3
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
        "-i", "--sif",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to singularity image file"
    )
    required.add_argument(
        "-p", "--partition",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="default partition(s) to use (e.g. 'par1' or 'par1,par2')"
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
        "--walltime",
        default='72',
        help=textwrap.dedent('''\
            wall time (hours) for main job runner
            [default: 48]
        ''')
    )
    optional.add_argument(
        "--sv-overlap-size",
        default=100,
        help=textwrap.dedent('''\
            maximum bp between SV annotations to be merged
            [default: 100]
        ''')
    ) 
    optional.add_argument(
        "--sv-callers",
        default=2,
        help=textwrap.dedent('''\
            number of SV callers needed to keep an SV annotation 
            excluding MANTA insertions [default: 2]
        ''')
    )
    optional.add_argument(
        "--sv-length",
        default=50,
        help="minimum length needed to keep an SV [default: 50]",
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
        "--alias",
        default="s3",
        help="minio client S3 storage alias [default: s3]"
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )

    args        = parser.parse_args()
    ref         = args.ref
    crams       = args.crams
    svars       = args.svars 
    bucket      = args.bucket
    outdir      = os.path.realpath(os.path.expanduser(args.out))
    snake_env   = args.snake_env
    partition   = args.partition
    email       = args.email
    account     = args.account
    walltime    = args.walltime
    sif         = str(Path(args.sif).resolve())
    profile     = args.profile
    remote      = args.remote.lower()
    alias       = args.alias
    sv_length   = args.sv_length
    sv_callers  = args.sv_callers
    sv_overlap  = args.sv_overlap_size 
   
    # confirm image exitst
    if os.path.isfile(sif):
        print("wags image found!")
    else:
        print("wags image not found - confirm location")
        sys.exit(1)
    
    # check if user ref available
    if ref not in config_d:
        avail_refs = sorted(config_d.keys())
        print(f"\nERROR: reference '{ref}' not found in available configs")
        # get possible matches
        poss_matches = [i for i in avail_refs if i.startswith(ref[:3])]
        # and handle
        if poss_matches:
            print(f"did you mean one of these: {', '.join(poss_matches)}")
        else:
            print("no similar ref found")

        # show how to list all available configs
        print("\nto see available species and assemblies (configs), run the following command:")
        print("    python prep_many.py --configs\n")

    main()
        
        
