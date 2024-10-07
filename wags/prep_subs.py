#!/usr/bin/env python3

import os
import re
import sys
import csv
import bz2
import gzip
import lzma
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

class RawFile(object):
    def __init__(self,filename):
        self.filename = filename
        if filename.endswith('.gz'):
            if self.is_gzipped(filename):
                self.handle = gzip.open(filename, 'rt')
            else:
                self.handle = open(filename,'r')
        elif filename.endswith('bz2'):
            self.handle = bz2.open(filename,'rt')
        elif filename.endswith('xz'):
            self.handle = lzma.open(filenaem,'rt')
        else:
            self.handle = open(filename,'r')

    def is_gzipped(self, filename):
        try:
            with gzip.open(filename, 'rt') as f_in:
                f_in.read(1)
            return True
        except gzip.BadGzipFile:
            return False

    def __enter__(self):
        return self.handle

    def __exit__(self,dtype,value,traceback):
        self.handle.close()

refs = [
    "canfam3","canfam4","UU_Cfam_GSD_1.0_ROSY",
    "goldenPath","Arabian","Shire","Thoroughbred",
    "tiger",
    "Fca126_mat1.0",
    "alpaca",
    "vicpac32",
    "ARS-UI_Ramb_v3.0",
    "ARS1.2"
]

# get config dir
prep_dir = Path(__file__).resolve().parent.parent
configs_dir = prep_dir / "pipelines" / "one_wag" / "configs"
# put available configs into a dictionary
config_d = {}
for species_dir in configs_dir.iterdir():
    if species_dir.is_dir():
        for config_f in species_dir.iterdir():
            if config_f.suffix == ".yaml":
                species = config_f.stem.replace("_config", "")
                config_d[species] = species_dir.name


replace_map = {
    "_2.fq.gz": "_1.fq.gz",
    "_R2.fq.gz": "_R1.fq.gz",
    "_2.fastq.gz": "_1.fastq.gz",
    "_R2.fastq.gz": "_R1.fastq.gz",
    "_R2_001.fastq.gz": "_R1_001.fastq.gz",
    "_R2_002.fastq.gz": "_R1_002.fastq.gz",
    "_R2_001.part_001.fastq.gz": "_R1_001.part_001.fastq.gz",
    "_R2_001.part_002.fastq.gz": "_R1_001.part_002.fastq.gz",
    "_R2_001.part_003.fastq.gz": "_R1_001.part_003.fastq.gz",
    "_R2_001.part_004.fastq.gz": "_R1_001.part_004.fastq.gz",
    "_R2.part_001.fastq.gz": "_R1.part_001.fastq.gz",
    "_R2.part_002.fastq.gz": "_R1.part_002.fastq.gz",
    "_R2.part_003.fastq.gz": "_R1.part_003.fastq.gz",
    "_R2.part_004.fastq.gz": "_R1.part_004.fastq.gz",
}

def extract_pu(s):
    """ extracts platform unit from fastq header """
   #with gzip.open(s,"rt") as f:
    with RawFile(str(s)) as f:
        head = f.readline().strip()
        if head.startswith("@SRR"): # handle SRR fastqs from NCBI
            return f"{head.split('.')[0][1:]}.NA.NA" # sample name in place of flowcell
        elif head.startswith("@HWI"): # older MiSeq (?) fastqs...
            head = head.split(":")
            if head[-1]:
                return f"{head[2]}.{head[3]}.{head[-1]}"
            else:
                return f"{head[2]}.{head[3]}.NA"
        elif ":" not in head:
            if head.startswith("@ERR"): # specific fix...unsure if consistent for all ERR
                return f"{head.split('.')[0][1:]}.NA.NA" # sample name in place of flowcell
        else: # or "standard" illumina headers...
            head = head.split(":")
            ycoord = head[-1].split(" ",1)[0] # modify headers with additional info (eg length=150)
            return f"{head[2]}.{head[3]}.{ycoord}"

def find_r2s(name, fq_dir):
    """ find R2 files based on provided substring """
    r2_patterns = [
        "_R2_001.fastq.gz", "_R2_002.fastq.gz", "_R2.fastq.gz", 
        "_2_001.fastq.gz", "_2_002.fastq.gz", "_2.fastq.gz",
        "_R2_001.fq.gz", "_R2_002.fq.gz", "_R2.fq.gz", "_2_001.fq.gz", 
        "_2_002.fq.gz", "_2.fq.gz",
        "_R2_001.part_001.fastq.gz", "_R2_001.part_002.fastq.gz",
        "_R2_001.part_003.fastq.gz", "_R2_001.part_004.fastq.gz",
        "_R2.part_001.fastq.gz", "_R2.part_002.fastq.gz",
        "_R2.part_003.fastq.gz", "_R2.part_004.fastq.gz"
    ]
    r1_patterns = [
        "_R1_001.fastq.gz", "_R1_002.fastq.gz", "_R1.fastq.gz", 
        "_1_001.fastq.gz", "_1_002.fastq.gz", "_1.fastq.gz",
        "_R1_001.fq.gz", "_R1_002.fq.gz", "_R1.fq.gz", "_1_001.fq.gz", 
        "_1_002.fq.gz", "_1.fq.gz",
        "_R1_001.part_001.fastq.gz", "_R1_001.part_002.fastq.gz",
        "_R1_001.part_003.fastq.gz", "_R1_001.part_004.fastq.gz",
        "_R1.part_001.fastq.gz", "_R1.part_002.fastq.gz",
        "_R1.part_003.fastq.gz", "_R1.part_004.fastq.gz"
    ]

    matched = []

    for f in Path(fq_dir).rglob(f"{name}*"):
        is_r2 = any(str(f).endswith(indicator) for indicator in r2_patterns)
        is_not_r1 = not any(str(f).endswith(indicator) for indicator in r1_patterns)
        
        if is_r2 and is_not_r1:
            matched.append(f)

    matched = [f for f in matched if not f.name.endswith('.md5')]

    return matched


def main():
    global profile
   
    d = defaultdict(dict)

    cols = ["breed","sample_name","readgroup_name",
            "fastq_1","fastq_2","library_name",
            "platform_unit","flowcell","run_date",
            "platform_name","sequencing_center"]

    with open(sample_meta) as ids:
        reader = csv.reader(ids, delimiter=',')
        next(reader, None)
        for line in reader:
            
            tmp = find_r2s(line[-1], fq_dir)
            # print the number of found fastq pairs/sample
            pairs = "pair" if len(tmp) == 1 else "pairs"
            print(f"found {len(tmp)} FASTQ {pairs} for {line[-1]}")
            sample_input = []  
            for i,v in enumerate(sorted(tmp)):
                
                platform_unit = extract_pu(v)
                cdate = datetime.fromtimestamp(os.path.getctime(v)).strftime('%Y-%m-%dT%H:%M:%S')                   
                
                # get first read
                v1 = None
                for R2,R1 in replace_map.items():
                    if R2 in os.path.basename(v):
                        v1 = os.path.join(
                            os.path.dirname(v),
                            os.path.basename(v).replace(R2, R1)
                        )
                        break

                sample_input.append(
                    [   
                        line[1],
                        line[0],
                        f"{line[0]}_{string.ascii_uppercase[i]}",
                        v1,
                        v,
                        f"illumina-{line[0]}",
                        platform_unit,
                        platform_unit.split(".")[0],
                        cdate,
                        "illumina",
                        "unknown"
                    ]
                )
            
            d[line[0]]['work_dir'] = os.path.join(outdir,line[1],line[0],ref)
            d[line[0]]['breed']    = line[1]
            d[line[0]]['df']       = pd.DataFrame(sample_input, columns=cols)
            
    # copy pipeline input and submission file
    for k,v in d.items():
        # sample_name and breed
        sample_name = k
        breed = v['breed']        

        os.makedirs(v['work_dir'], exist_ok=True)

        # create jobs dir if not exist
        jobs = os.path.join(v['work_dir'],f"{profile}_logs")
        if not os.path.exists(jobs):
            os.makedirs(jobs)
   
        # write sample input to working directory
        with open(os.path.join(v['work_dir'],"input.tsv"),"w") as out:
            v['df'].to_csv(out,sep='\t',index=False)
       
        # input templates
        pipeline  = "one_wag"
        rules     = "rules"
        snake_n   = "one_wag.smk"
        config_n  = f"{ref}_config.yaml"
        profile_n = f"{profile}.go_wags"
        # switch to money templates if true - NEEDS TO BE UPDATED STILL
        if money:
            pipeline  = "only_wag"
            rules     = "rules"
            snake_n   = "only_wag.smk"
            config_n  = f"{ref}_only_wag.yaml"
 
        # copy snakefile, rules, config, and profile to working dir
       #prep_dir = os.path.dirname(os.path.realpath(__file__))
        smk = os.path.join(
           #str(Path(prep_dir).parents[0]),
            prep_dir,
            "pipelines",
            pipeline,
            f"inputs/{remote}/{bqsr}",
            snake_n
        )
        
        rules = os.path.join(
            prep_dir,
            "pipelines",
            pipeline,
            f"inputs/{remote}/{bqsr}",
            rules
        )
        config = os.path.join(
            prep_dir,
            "pipelines",
            pipeline,
            "configs",
            config_d[ref],
            config_n
        )
        # selected reference not in container, presumed prep_custom_ref.py already
        # executed
        # NEEDS TO BE FIXED ALONG WITH prep_custom_ref.py TO FORCE THE CONFIG DIR
        # AS SPECIES AND CONFIG ARE PLACED WITH OTHERS
        if ref not in refs:
            tmp = glob.glob(f"{ref_dir}/**/{ref}_config.yaml",recursive=True)
            assert tmp, f"config not found for {ref}, ensure prep_custom_ref.py ran successfully and check ref_dir"
            config = tmp[0]
    
        profile_dir = os.path.join(
            prep_dir,
            f"profiles/{profile}",
            profile_n
        )
        
        src = os.path.join(
            prep_dir,
            "pipelines",
            pipeline,
            "src"
        )

        # modify config file
        with open(config) as f:
            doc = yaml.safe_load(f)
        
        # get config values
        species  = doc['species']
        ref_dict = doc['ref_dict']
        fasta    = doc['ref_fasta']
        
        # update sif and other cli args
        doc['profile'] = profile
        doc['sif']     = sif
        doc['alias']   = alias
        doc['bucket']  = bucket
        doc['left_align'] = left_align
        doc['tmp_dir']['sort_tmp']  = os.path.join(
            outdir,".sort",breed,sample_name,".tmp"
        )
        doc['tmp_dir']['fastq_tmp'] = os.path.join(
            outdir,".fastq",breed,sample_name,".tmp"
        )
        
        # for private variant analysis, add pop.vcf and common.vcf path
        if money:
            doc['pop_vcf']     = pop
            doc['common_vcf']  = common
            doc['allele_freq'] = float(allele_freq)
            doc['tmp_dir']['select_variants_to_table'] = os.path.join(
                outdir,".select",breed,sample_name,".tmp"
            )
            doc['tmp_dir']['unfilt_gather_vcf'] = os.path.join(
                outdir,".gather",breed,sample_name,".tmp"
            )

        # interval optional updates
        doc['nrun_length']  = int(nrun_len)
        doc['scatter_size'] = int(scat_size)
        
        # dump
        with open(os.path.join(v['work_dir'],config_n),'w') as out:
            yaml.dump(doc,out,sort_keys=False)
   
        input_names = [snake_n,"rules",profile_n,"src"]
        dst_files = [os.path.join(v['work_dir'],i) for i in input_names]
        
        for i in zip([smk,rules,profile_dir,src],dst_files):
            if not os.path.exists(i[1]):
                if os.path.isfile(i[0]):
                    shutil.copy(i[0],i[1])
                else:
                    shutil.copytree(i[0],i[1])
                    # modify profile profile-submit.py for user-supplied parition
                    if f"{profile}.go_wags" in i[0]:
                        if profile == 'lsf':
                            job_sub = os.path.join(i[1],f"{profile}_submit.py")
                        else:
                            job_sub = os.path.join(i[1],f"{profile}-submit.py")
                        with fileinput.FileInput(job_sub,inplace=True,backup=".bak") as file:
                            for line in file:
                                line = line.replace("DUMMY_PAR",partition)
                                line = line.replace("DUMMY_ACC",account)
                                print(line,end='')
                    
        # submission destination
        job_name = snake_n.split('.')[0]
        submiss = os.path.join(v['work_dir'],f"{breed}_{sample_name}.{job_name}.{profile}")
        
        # SBATCH directives 
        default_header = (
            "#!/bin/bash -l\n"
            f"#SBATCH -t {walltime}:00:00\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks-per-node=1\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --mem=12gb\n"
            "#SBATCH --mail-type=ALL\n"
            f"#SBATCH --mail-user={email}\n"
            f"#SBATCH --job-name {v['breed']}_{k}.{job_name}.{profile}\n"
            f"#SBATCH -o slurm_logs/%j.{v['breed']}_{k}.{job_name}.out\n"
            f"#SBATCH -e slurm_logs/%j.{v['breed']}_{k}.{job_name}.err\n"
            f"#SBATCH -A {account}\n"
            f"#SBATCH -p {partition}\n"
        )             
        lsf_header = (
            "#!/bin/bash -l\n"
            f"#BSUB -W {walltime}:00\n"
            "#BSUB -n 1\n"
            "#BSUB -R span[hosts=1]\n"
            "#BSUB -R rusage[mem=12GB]\n"
            f"#BSUB -J {v['breed']}_{k}.{job_name}.{profile}\n"
            f"#BSUB -o {profile}_logs/%J.{v['breed']}_{k}.{job_name}.out\n"
            f"#BSUB -e {profile}_logs/%J.{v['breed']}_{k}.{job_name}.err\n"
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

            if ref not in refs:
                print(f"REF_DIR={ref_dir}",end="",file=f)
            if money:
                print(f"\nPOP_VCF={os.path.dirname(pop)}",end="",file=f)

            print(
                textwrap.dedent(
                    f"""
                    TMP_DIR=/share/stern/mwvandew/tmp
                    FQ_DIR={fq_dir}
                    PROC_DIR={outdir} 
                    """
                ),file=f
            )

            # if ref in container, extract dictionary or generate from fasta 
            # if does not exist
            if ref in refs:
                print("# extract reference dict from container",end="",file=f)
                print(
                    textwrap.dedent(
                        f"""
                        singularity exec --bind $PWD {sif} \\
                            cp /home/refgen/{species}/{ref}/{ref_dict} $PWD
                        """
                    ),file=f
                )
                # added gap file from TK to use for T2T reference when public
                if "Thoroughbred" in ref:
                    print("# extract T2T gaps file from container",end="",file=f)
                    print(
                        textwrap.dedent(
                            f"""
                            singularity exec -B $PWD {sif} \\
                                cp /home/refgen/horse/Thoroughbred/resources/TB-T2T_gaps.csv $PWD
                            """
                        ),file=f
                    )

            print(
                textwrap.dedent(
                    f"""
                    snakemake -s {snake_n} \\
                        --use-singularity \\
                        --singularity-args "-B $TMP_DIR,$PWD,$REF_DIR,$POP_VCF,$FQ_DIR,$PROC_DIR" \\
                        --profile {profile_n} \\
                        --configfile {config_n} \\
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
                            {alias}/{bucket}/wgs/{breed}/{sample_name}/{ref}/{profile}_logs/
                        """
                    ),file=f
                )

    plural = "sample"
    if len(d) > 1:
        plural = "samples"
    print(f"{len(d)} {plural} setup for processing")
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="wags",
        add_help=False,
        description=(
            "       ___           ___           ___           ___      \n"
            " ONE  /__/\         /  /\         /  /\         /  /\     \n"
            "     _\_ \:\       /  /::\       /  /:/_       /  /:/_    \n"
            "    /__/\ \:\     /  /:/\:\     /  /:/ /\     /  /:/ /\   \n"
            "   _\_ \:\ \:\   /  /:/~/::\   /  /:/_/::\   /  /:/ /::\  \n"
            "  /__/\ \:\ \:\ /__/:/ /:/\:\ /__/:/__\/\:\ /__/:/ /:/\:\ \n"
            "  \  \:\ \:\/:/ \  \:\/:/__\/ \  \:\ /~~/:/ \  \:\/:/~/:/ \n"
            "   \  \:\ \::/   \  \::/       \  \:\  /:/   \  \::/ /:/  \n"
            "    \  \:\/:/     \  \:\        \  \:\/:/     \__\/ /:/   \n"
            "     \  \::/       \  \:\        \  \::/        /__/:/    \n"
            "      \__\/         \__\/         \__\/         \__\/   \n\n"
            "wags generates all required input to process FASTQs to GVCF \n"
            "following GATK best practices. For each sample (and associated \n"
            "FASTQ pair), wags outputs a directory structure organized by \n"
            "breed, wherein GATK pipeline input are contained by sample ID."
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
        "-m", "--meta",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="csv of meta data to associate sample names and fastqs"
    )
    required.add_argument(
        "-f", "--fastqs",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to fastq directory"
    )
    required.add_argument(
        "-r", "--ref",
        default="canfam4",
        metavar="",
        help=textwrap.dedent(f'''\
            select reference to use: 
                {", ".join(refs)}
            if using custom reference, ensure provided name
            is exact match to name (--ref) used with 
            prep_custom_ref.py
        ''')
    )
    required.add_argument(
        "-o", "--out",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to out dir"
    )
    required.add_argument(
        "-b", "--bucket",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help=textwrap.dedent('''\
            bucket (with --remote s3) or results 
            directory name. if using sftp, bucket is the
            host name and path to output director
            (e.g. hostname/path/to/dir)
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
    optional.add_argument(
        "--walltime",
        default='48',
        help=textwrap.dedent('''\
            wall time (hours) for main job runner
            [default: 48]
        ''')
    )
    optional.add_argument(
        "--ref-dir",
        default=os.path.join(os.path.expanduser("~"),".wags/"),
        help=textwrap.dedent('''\
            path to custom reference directory generated by
            prep_custom_ref.py - multiple assemblies from
            the same species must have different names 
            [default ~/.wags]
        '''),
    )
    optional.add_argument(
        "--nrun-length",
        default='50',
        help=textwrap.dedent('''\
            maximum number of contiguous missing
            bases to tolerate for generating intervals
            [default: 50]
        ''')
    )
    optional.add_argument(
        "--scatter-size",
        default='50',
        help=textwrap.dedent('''\
            maximum number of intervals across which
            haplotype calling and bqsr (if desired) occur
            [default: 50]
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
        "--alias",
        default="s3",
        help="minio client S3 storage alias [default: s3]"
    )
    optional.add_argument(
        "--no-bqsr",
        action="store_true",
        help="turn off bqsr [default: bqsr on]"
    )
    optional.add_argument(
        "--left-align",
        action="store_true",
        help="left align analysis-ready bam [default: off]"
    )
    optional.add_argument(
        "--sv",
        action="store_true",
        help="call structural variants [default: off]"
    )
    optional.add_argument(
        "--money",
        action="store_true",
        help="switch to setup private SNP pipeline"
    )
    optional.add_argument(
        "--pop",
        default=None,
        help="path to population variants vcf (only wags pipeline) [default: None]"
    )
    optional.add_argument(
        "--allele-freq",
        default='0.005',
        help=textwrap.dedent('''\
            threshold to define common variants compared
            to the population vcf (only relevant with 
            --money) [default: 0.005]
        ''')
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )

    args = parser.parse_args()
    sample_meta = os.path.realpath(os.path.expanduser(args.meta))
    fq_dir      = os.path.realpath(os.path.expanduser(args.fastqs))
    outdir      = os.path.realpath(os.path.expanduser(args.out))
    bucket      = args.bucket
    snake_env   = args.snake_env
    partition   = args.partition
    email       = args.email
    account     = args.account
    walltime    = args.walltime
    sif         = args.sif
    ref         = args.ref
    ref_dir     = os.path.realpath(os.path.expanduser(args.ref_dir))
    alias       = args.alias
    money       = args.money
    pop         = args.pop
    allele_freq = args.allele_freq
    nrun_len    = args.nrun_length
    scat_size   = args.scatter_size
    profile     = args.profile
    remote      = args.remote.lower()
    no_bqsr     = args.no_bqsr
    left_align  = args.left_align
    sv_call     = args.sv

    # QUICK FIX FOR goldenPath - NEED TO ADJUST CONTAINER TO BE horse/goldenpath
   #if "golden" in ref:
   #    ref = "goldenPath"

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
        print("    python prep_subs.py --species\n")

        sys.exit(1)

    # get bqsq or no bqsr
    bqsr = "bqsr"
    if no_bqsr:
        bqsr = "no_bqsr"

    if money:
       #bqsr = "recal"
        # require path to pop
        if pop is None:
            parser.error("--money requires --pop")
        else:
            pop = os.path.realpath(os.path.expanduser(args.pop))
            # generate path to common vars assuming pop vcf was generated by
            # the joint pipeline
            p = Path(pop)
            # create directory common_vars within pop vcf directory
            base = os.path.join(str(Path(*p.parts[:-1])),"common_vars")
            os.makedirs(base,exist_ok=True)
            # generate the common vcf name and set variable
            common_name = f"{p.parts[-1].rsplit('.',2)[0]}.af_nonmajor.vcf.gz"
            common = os.path.join(base,common_name)
            

    # if non default sif location
    if sif != os.path.join(os.path.expanduser("~"),".sif/wags.sif"):
        sif = os.path.realpath(os.path.expanduser(sif))
    
    # confirm image exitst
    if os.path.isfile(sif):
        print("wags image found!")
    else:
        print("wags image not found - confirm location")
        sys.exit(1)

    # check if scratch dir exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} created!")
        
    main()
        
        
