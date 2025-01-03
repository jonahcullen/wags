#!/usr/bin/env python3

import os
import sys
import csv
import gzip
import yaml
import glob
import shutil
import string
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
configs_dir = prep_dir / "pipelines" / "many_wags" / "configs"
# put available configs into a dictionary
config_d = {}
for species_dir in configs_dir.iterdir():
    if species_dir.is_dir():
        for config_f in species_dir.iterdir():
            if config_f.suffix == ".yaml":
                species = config_f.stem.replace("_config", "")
                config_d[species] = species_dir.name

# https://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings
def common_prefix(strings):
    '''
    find the longest string that is a prefix of all the strings.
    '''
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    # ensure prefex ends with directory and not split file name
    if os.path.isfile(prefix) or os.path.isdir(prefix):
        return prefix
    else:
        return os.path.dirname(prefix)

# NOTE - should also have checks that each sample and gvcf mapping are unique
def validate_mapping(f):
    '''
    quick check that the sample-gvcf mapping contains two tab-delimited columns
    '''
    prefs = []
    with open(f, 'r') as infile:
        for line in infile:
            cols = line.split('\t')
            col_length = len(cols)
            if col_length == 2:
                prefs.append(cols[1])
                continue
            else:
                sys.exit(textwrap.dedent('''\n
                    ill-formatted input (-g/--gvcf) - ensure each row contains
                    two columns, sample followed by path to gvcf, separated
                    by a tab
                '''))
        return common_prefix(prefs)

def main():
    
    # load known config
    prep_path = Path(__file__).resolve()
<<<<<<< Updated upstream:wags/prep_many.py
    config_path = prep_path.parent.parent / f"pipelines/many_wags/configs/{config_d[ref]}/{ref}_config.yaml"
=======
    #config_path = prep_path.parent.parent / f"pipelines/many_wags/configs/{config}_config.yaml"
>>>>>>> Stashed changes:wags/prep_joint.py
    try:
        with open(config) as f:
            doc = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"{config} does not exist - ensure correct path")
        
    # update config with user info, mapping cohort, interval options, and temp dirs
    doc['sif'] = sif
    gvcf_name = Path(gvcfs).name
    doc['joint_cohort'] = gvcf_name
    doc['date'] = datetime.today().strftime('%Y%m%d')
    if Path(gvcfs).is_file():
        if validate_mapping(gvcfs):
            doc['joint_cohort'] = str(Path(outdir) / gvcf_name)
           #doc['joint_cohort'] = os.path.join(outdir, os.path.basename(gvcfs))
    else:
        sys.exit("mapping file (--gvcfs) does not exist - ensure correct path")
    doc['tmp_dir']['sites_only_gather_vcf'] = str(Path(outdir) / ".sites_gather")
    doc['tmp_dir']['unfilt_gather_vcf']     = str(Path(outdir) / ".unfilt_gather")

    # interval optional updates
    doc['anchor_type']     = anchor_type
    doc['nrun_length']     = int(nrun_length)
    doc['interval_length'] = int(ival_length)
    doc['scatter_size']    = int(scat_size)
    
    # other
    doc['bucket']             = bucket
    doc['alias']              = alias
    doc['profile']            = profile
    doc['snp_filter_level']   = snp_val
    doc['indel_filter_level'] = indel_val
    if phase is not None:
        doc['phasing']  = true
        doc['link_map'] = phase

    # get values from config
    ref_fasta = doc['ref_fasta']
    ref_dict  = doc['ref_dict']
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
    pipeline  = "many_wags"
    rules     = "rules"
    snake_n   = "many_wags.smk"
    profile_n = f"{profile}.go_wags"

    # copy snakefile, rules, config, and profile to working dir
    prep_dir = Path(__file__).resolve().parent

    smk = (
        prep_dir.parents[0]
        / "pipelines/many_wags"
        / f"inputs/{remote}/{recal}"
        / snake_n
    ) 
    
    rules = (
        prep_dir.parents[0]
        / "pipelines/many_wags"
        / f"inputs/{remote}/{recal}"
        / rules
    )
    
    profile_dir = (
        prep_dir.parents[0]
        / f"profiles/{profile}"
        / profile_n
    )
    
    src = (
        prep_dir.parents[0]
        / "pipelines/many_wags"
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

    # copy gvcfs to output
    dst_gvcf = Path(outdir) / gvcf_name
    if not dst_gvcf.is_file():
        shutil.copy(gvcfs, dst_gvcf)

    # submission destination
    job_name = snake_n.split('.')[0]
    gvcf_base = Path(gvcfs).stem
    submiss = os.path.join(outdir, f"{gvcf_base}_{ref}.{job_name}.{profile}")
     
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
        f"#SBATCH --job-name {gvcf_base}_{ref}.{job_name}\n"
        f"#SBATCH -o slurm_logs/%j.{gvcf_base}_{ref}.{job_name}.out\n"
        f"#SBATCH -e slurm_logs/%j.{gvcf_base}_{ref}.{job_name}.err\n"
        f"#SBATCH -A {account}\n"
        f"#SBATCH -p {partition}\n"
    )             
    lsf_header = (
        "#!/bin/bash -l\n"
        f"#BSUB -W {walltime}:00\n"
        "#BSUB -n 1\n"
        "#BSUB -R span[hosts=1]\n"
        "#BSUB -R rusage[mem=12GB]\n"
        f"#BSUB -J {gvcf_base}_{ref}.{job_name}\n"
        f"#BSUB -o {profile}_logs/%J.{gvcf_base}_{ref}.{job_name}.out\n"
        f"#BSUB -e {profile}_logs/%J.{gvcf_base}_{ref}.{job_name}.err\n"
        f"#BSUB -q {partition}\n"
        f"#BSUB -B -N -u {email}\n"
    )

    # job submission body 
    with open(submiss, "w") as f:
        if profile == 'lsf':
            print(lsf_header, file=f)
        else:
            print(default_header, file=f)
        print("set -e\n",file=f)
        print(f"conda activate {snake_env}",file=f)
        if profile != 'lsf':
            print("cd $SLURM_SUBMIT_DIR\n",file=f)
        
        if "chrom" in anchor_type:
            print("# extract reference dict from container",end="",file=f)
            print(
                textwrap.dedent(
                    f"""
                    singularity exec --bind $PWD {sif} \\
                        cp /home/refgen/{species}/{ref}/{ref_dict} $PWD
                    """
                ),file=f
            )

       #if ref not in refs:
       #    print(f"FASTA_DIR={Path(ref_fasta).parent}\n", end="", file=f)
       #    print(f"DICT_DIR={Path(ref_dict).parent}\n", end="", file=f)
       #    print(f"GTF_DIR={Path(ref_gtf).parent}\n", end="", file=f)

       #    print(
       #        textwrap.dedent(
       #            f"""
       #            snakemake -s {snake_n} \\
       #                --use-singularity \\
       #                --singularity-args "-B $PWD,$FASTA_DIR,$DICT_DIR,$GTF_DIR" \\
       #                --profile {profile_n} \\
       #                --configfile {os.path.basename(ref)} \\
       #                --keep-going
       #            """
       #        ),file=f
       #    ) 
       #else:     # NOTE - THE GVCF DIR IS HARDCODED BELOW AND NEEDS TO BE FIXED BY CHECKING THE GVCF FILE       
        print(
            textwrap.dedent(
                f"""
                GVCF_DIR={validate_mapping(gvcfs)}
                snakemake -s {snake_n} \\
                    --use-singularity \\
                    --singularity-args "-B $PWD,$GVCF_DIR" \\
                    --profile {profile_n} \\
                    --configfile {ref}_config.yaml \\
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
        "-g", "--gvcfs",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent(f"""\
            path to file with mapping of sample name to
            gvcf in tab delimited format with no header
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
        "--anchor-type",
        default='nruns',
        help=textwrap.dedent('''\
            anchor type for generating intervals.
            'intergenic' requires an annotation file
            (e.g. gff/gtf) to be included during
            config_joint.py setup [options: nruns, intergenic, chroms]
        ''')
    )
    optional.add_argument(
        "--nrun-length",
        default='50',
        help=textwrap.dedent('''\
            maximum number of contiguous missing
            bases to tolerate. only relevant for
            anchor-type nruns
            allowed [default: 50]
        ''')
    )
    optional.add_argument(
        "--interval-length",
        default='10000000',
        help=textwrap.dedent('''\
            desired invterval target lengths.
            this does not guarantee a uniform
            length and is largely dependent on
            interval design choice (--anchor-type)
            and genome completeness [default: 10000000]
        ''')
    )
    optional.add_argument(
        "--scatter-size",
        default='50',
        help=textwrap.dedent('''\
            maximum number of intervals for vep-base
            annotation [default: 50]
        ''')
    )
    optional.add_argument(
        "--snp-tranche-level",
        default="99.0",
        help="truth sensitivity level for snp filtering (vqsr only) [default: 99.0]",
    )
    optional.add_argument(
        "--indel-tranche-level",
        default="99.0",
        help="truth sensitivity level for indel filtering (vqsr only) [default: 99.0]",
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
        "--phase",
        metavar="\b",
        help="path to linkage map to phase final vcf [default: false]"
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )

    args        = parser.parse_args()
    ref         = args.ref
    gvcfs       = args.gvcfs
    bucket      = args.bucket
    snps        = args.vqsr_snps
    indels      = args.vqsr_indels
    outdir      = os.path.realpath(os.path.expanduser(args.out))
    snake_env   = args.snake_env
    partition   = args.partition
    email       = args.email
    account     = args.account
    walltime    = args.walltime
    sif         = args.sif
    profile     = args.profile
    remote      = args.remote.lower()
    alias       = args.alias
    anchor_type = args.anchor_type.lower()
    nrun_length = args.nrun_length
    scat_size   = args.scatter_size
    ival_length = args.interval_length
    snp_val     = args.snp_tranche_level
    indel_val   = args.indel_tranche_level
    phase       = args.phase

    # based on snps/indels, define pipeline
    recal = "vqsr_none"
    if snps and indels:
        recal = "vqsr_all"
    elif snps and not indels:
        recal = "vqsr_snps"
    
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
        
        
