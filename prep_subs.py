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
#remotes = ["local","s3","sftp"]

def extract_pu(s):
    """ extracts platform unit from fastq header """
    with gzip.open(s,"rt") as f:
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
            return f"{head[2]}.{head[3]}.{head[-1]}"


#def main(dog_meta,outdir,fq_dir,ref,profile):
def main():
    global profile
   
    d = defaultdict(dict)

    cols = ["breed","sample_name","readgroup_name",
            "fastq_1","fastq_2","library_name",
            "platform_unit","flowcell","run_date",
            "platform_name","sequencing_center"]
    
    with open(dog_meta) as ids:
        reader = csv.reader(ids, delimiter=',')
        next(reader, None)
        for line in reader:
            
            # there is a problem with this strategy if the sample name
            # contains _R1 or _R2 
            # find associated fastqs
            tmp = []
            for root,dirs,files in os.walk(fq_dir):
                for f in files:
                    if f.startswith(line[-1]) and not f.endswith("md5") and "_R2" in f:
                        tmp.append(os.path.join(root,f))
            first = "_R1"           
            second = "_R2"
            # special circumstance where fastqs labelled as sample_{1,2}.fastq.gz
            # there is DEFINITELY a better to deal with finding fastqs...
            if not tmp:
                for root,dirs,files in os.walk(fq_dir):
                    for f in files:
                        if f.startswith(line[-1]) and f.endswith("_2.fastq.gz"): 
                            tmp.append(os.path.join(root,f))
                first = "_1.fastq.gz"
                second = "_2.fastq.gz"
            # add fastq information for each pair per sample
            dog_input = []  
            for i,v in enumerate(sorted(tmp)):
                
                platform_unit = extract_pu(v)
                cdate = datetime.fromtimestamp(os.path.getctime(v)).strftime('%Y-%m-%dT%H:%M:%S')                   
   
                # where second read is v, first read is
                v1 = os.path.join(
                    os.path.dirname(v),
                    os.path.basename(v).replace(second,first)
                )
 
                dog_input.append(
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
            
            d[line[0]]['work_dir']    = os.path.join(outdir,line[1],line[0],ref)
            d[line[0]]['breed']       = line[1]
            d[line[0]]['df']          = pd.DataFrame(dog_input, columns=cols)
            
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
   
        # write dog input to working directory
        with open(os.path.join(v['work_dir'],"input.tsv"),"w") as out:
            v['df'].to_csv(out,sep='\t',index=False)
       
        # input templates
        pipeline  = "process_wgs"
        rules     = "rules"
        snake_n   = "one_wags.smk"
        config_n  = f"{ref}_config.yaml"
        profile_n = f"{profile}.go_wags"
        # switch to money templates if true - NEEDS TO BE UPDATED STILL
        if money:
            pipeline  = "make_money"
            rules     = "rules"
            snake_n   = "go_make_money.smk"
            config_n  = f"{ref}_money.yaml"
           #profile_n = "slurm.go_money"
 
        # copy snakefile, rules, config, and profile to working dir
        smk = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "pipelines",
            pipeline,
            f"inputs/{remote}/{bqsr}",
            snake_n
        )
        
        rules = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "pipelines",
            pipeline,
            f"inputs/{remote}/{bqsr}",
            rules
        )
        
        config = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "pipelines",
            pipeline,
            "configs",
            config_n
        )
        # selected reference not in container, presumed prep_custom_ref.py already
        # executed
        if ref not in refs:
            tmp = glob.glob(f"{ref_dir}/**/{ref}_config.yaml",recursive=True)
            assert tmp, f"config not found for {ref}, ensure prep_custom_ref.py ran successfully and check ref_dir"
            config = tmp[0]
    
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
        # for private variant analysis, add pop.vcf and common.vcf path
        if money:
            doc['tmp_dir']['sort_tmp']              = os.path.join(outdir,".sort",breed,sample_name,".tmp")
            doc['tmp_dir']['sites_only_gather_vcf'] = os.path.join(outdir,".sort",breed,sample_name,".tmp")
            doc['pop_vcf']    = pop
            doc['common_vcf'] = common
        else:
            doc['sort_tmp']  = os.path.join(outdir,".sort",breed,sample_name,".tmp")
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
        header = (
            "#!/bin/bash -l\n"
            "#SBATCH -t 60:00:00\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks-per-node=1\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --mem=12gb\n"
            "#SBATCH --mail-type=ALL\n"
            f"#SBATCH --mail-user={email}\n"
            f"#SBATCH --job-name {v['breed']}_{k}.{job_name}.slurm\n"
            f"#SBATCH -o slurm_logs/%j.{v['breed']}_{k}.{job_name}.out\n"
            f"#SBATCH -e slurm_logs/%j.{v['breed']}_{k}.{job_name}.err\n"
            f"#SBATCH -A {account}\n"
        )             
     
        # job submission body 
        with open(submiss, "w") as f:
            print(header, file=f)
            print("set -e\n",file=f)
            print(f"conda activate {snake_env}",file=f)
            print("cd $SLURM_SUBMIT_DIR\n",file=f)

            if ref not in refs:
                print(f"REF_DIR={ref_dir}",end="",file=f)
            if money:
                print(f"\nPOP_VCF={os.path.dirname(pop)}",end="",file=f)

            print(
                textwrap.dedent(
                    f"""
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
 
            print(
                textwrap.dedent(
                    f"""
                    snakemake -s {snake_n} \\
                        --use-singularity \\
                        --singularity-args "-B $PWD,$REF_DIR,$POP_VCF,$FQ_DIR,$PROC_DIR" \\
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
                        mc cp --recursive ./{profile}_logs/ \\
                            {alias}/{bucket}/wgs/{breed}/{sample_name}/{ref}/
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
            "      ___           ___           ___           ___      \n"
            "     /__/\         /  /\         /  /\         /  /\     \n"
            "    _\_ \:\       /  /::\       /  /:/_       /  /:/_    \n"
            "   /__/\ \:\     /  /:/\:\     /  /:/ /\     /  /:/ /\   \n"
            "  _\_ \:\ \:\   /  /:/~/::\   /  /:/_/::\   /  /:/ /::\  \n"
            " /__/\ \:\ \:\ /__/:/ /:/\:\ /__/:/__\/\:\ /__/:/ /:/\:\ \n"
            " \  \:\ \:\/:/ \  \:\/:/__\/ \  \:\ /~~/:/ \  \:\/:/~/:/ \n"
            "  \  \:\ \::/   \  \::/       \  \:\  /:/   \  \::/ /:/  \n"
            "   \  \:\/:/     \  \:\        \  \:\/:/     \__\/ /:/   \n"
            "    \  \::/       \  \:\        \  \::/        /__/:/    \n"
            "     \__\/         \__\/         \__\/         \__\/   \n\n"
            "wags generates all required input to process FASTQs to gVCF \n"
            "following GATK best practices. For each sample (and associated \n"
            "FASTQ pair), wags outputs a directory structure organized by \n"
            "breed, wherein GATK pipeline input are contained by sample ID."
        ),
       #formatter_class=argparse.RawDescriptionHelpFormatter
        formatter_class=argparse.RawTextHelpFormatter
    )
    
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
            select reference to use: {", ".join(refs)}.
            if using custom reference, ensure provided name
            is exact matche to name (--ref) used with 
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
            bucket name. if using sftp, bucket is the
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
        "--ref-dir",
        default=os.path.join(os.path.expanduser("~"),".wags/"),
        help=textwrap.dedent('''\
            path to custom reference directory generated by
            prep_custom_ref.py - multiple references from
            the same species must have different names 
            [default ~/.wags]
        '''),
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
        help="left align analysis-ready bam  [default: off]"
    )
    optional.add_argument(
        "--money",
        action="store_true",
        help="switch to setup private SNP pipeline"
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )
   #optional.add_argument(
   #    "--common",
   #    default=None,
   #    help="path to common variants vcf (money pipeline) [default: None]"
   #)
    optional.add_argument(
        "--pop",
        default=None,
        help="path to population variants vcf (money pipeline) [default: None]"
    )

    args = parser.parse_args()
    dog_meta   = os.path.realpath(os.path.expanduser(args.meta))
    fq_dir     = os.path.realpath(os.path.expanduser(args.fastqs))
    outdir     = os.path.realpath(os.path.expanduser(args.out))
    bucket     = args.bucket
    snake_env  = args.snake_env
    partition  = args.partition
    email      = args.email
    account    = args.account
    sif        = args.sif
    ref        = args.ref
   #ref_dir    = os.path.realpath(os.path.expanduser(args.ref_dir))
    ref_dir    = os.path.realpath(os.path.expanduser(args.ref_dir))
    alias      = args.alias
    money      = args.money
   #common     = args.common
   #pop        = os.path.realpath(os.path.expanduser(args.pop))
    pop        = args.pop
    profile    = args.profile
    remote     = args.remote.lower()
    no_bqsr    = args.no_bqsr
    left_align = args.left_align

    # QUICK FIX FOR goldenPath - NEED TO ADJUST CONTAINER TO BE horse/goldenpath
    if "golden" in ref:
        ref = "goldenPath"

    # get bqsq or no bqsr
    bqsr = "bqsr"
    if no_bqsr:
        bqsr = "no_bqsr"

    if money:
        bqsr = "recal"
        # require path to pop
        if pop is None:
            parser.error("--money requires --pop")
        else:
            pop = os.path.realpath(os.path.expanduser(args.pop))
            # generate path to common vars assuming pop vcf was generated by
            # the joint pipeline
            p = pathlib.Path(pop)
            # create directory common_vars within pop vcf directory
            base = os.path.join(str(pathlib.Path(*p.parts[:-1])),"common_vars")
            os.makedirs(base,exist_ok=True)
            # generate the common vcf name and set variable
            common_name = f"{p.parts[-1].rsplit('.',2)[0]}.af_nonmajor.vcf.gz"
            common = os.path.join(base,common_name)
            

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

    # check if scratch dir exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} created!")
        
   #main(dog_meta,outdir,fq_dir,ref,profile)
    main()
        
        
