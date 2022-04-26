#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 20:43:49 2021

@author: jonahcullen
"""


import argparse
import os
import sys
import csv
from collections import defaultdict
import gzip
from datetime import datetime
import shutil
import string
import textwrap
import pandas as pd
import wget
import yaml


refs = ["canfam3", "canfam4"]

def extract_pu(s):
    """ extracts platform unit from fastq header """
    with gzip.open(s, "rt") as f:
        head = f.readline().strip()
        if head.startswith("@SRR"): # handle SRR fastqs from NCBI
            return f"{head.split('.')[0][1:]}.NA.NA"
        elif head.startswith("@HWI"): # older MiSeq (?) fastqs...
            head = head.split(":")
            if head[-1]:
                return f"{head[2]}.{head[3]}.{head[-1]}"
            else:
                return f"{head[2]}.{head[3]}.NA"
        else: # or "standard" illumina headers...
            head = head.split(":")
            return f"{head[2]}.{head[3]}.{head[-1]}"


def main(dog_meta, outdir, fq_dir, ref):
   
    d = defaultdict(dict)

    cols = ["breed","sample_name","readgroup_name",
            "fastq_1","fastq_2","library_name",
            "platform_unit","flowcell","run_date",
            "platform_name","sequencing_center"]
    
    with open(dog_meta) as ids:
        reader = csv.reader(ids, delimiter=',')
        next(reader, None)
        for line in reader:
            
            # find associated fastqs
            tmp = []
            for root,dirs,files in os.walk(fq_dir):
                for f in files:
                    if f.startswith(line[-1]) and not f.endswith("md5") and "_R2" in f:
                        tmp.append(os.path.join(root,f))
            
            # add fastq information for each pair per sample
            dog_input = []  
            for i,v in enumerate(sorted(tmp)):
                
                platform_unit = extract_pu(v)
                cdate = datetime.fromtimestamp(os.path.getctime(v)).strftime('%Y-%m-%dT%H:%M:%S')                   
    
                dog_input.append(
                    [   
                        line[1],
                        line[0],
                        f"{line[0]}_{string.ascii_uppercase[i]}",
                        v.replace("_R2","_R1"),
                        v,
                        f"illumina-{line[0]}",
                        platform_unit,
                        platform_unit.split(".")[0],
                        cdate,
                        "illumina",
                        "unknown"
                    ]
                )
            
            d[line[0]]["work_dir"]    = os.path.join(outdir,line[1],line[0],ref)
            d[line[0]]["breed"]       = line[1]
            d[line[0]]["df"]          = pd.DataFrame(dog_input, columns=cols)
            
            # save breed and sample name
            breed       = line[1]
            sample_name = line[0]
            
    # copy pipeline input and slurm file
    for k,v in d.items():
        os.makedirs(v["work_dir"], exist_ok=True)

        # create jobs dir if not exist
        jobs = os.path.join(v["work_dir"],"slurm_logs")
        if not os.path.exists(jobs):
            os.makedirs(jobs)
   
        # write dog input to working directory
        with open(os.path.join(v["work_dir"],"input.tsv"),"w") as out:
            v["df"].to_csv(out,sep='\t',index=False)
       
        # input templates
        pipeline  = "GoProcessWGS"
        rules     = "rules"
        snake_n   = "go_process_wgs.smk"
        config_n  = f"{ref}_config.yaml"
        profile_n = "slurm.go_wgs"
        # switch to money templates if true
        if money:
            pipeline  = "GoMakeMoney"
            rules     = "rules"
            snake_n   = "go_make_money.smk"
            config_n  = f"{ref}_money.yaml"
            profile_n = "slurm.go_money"
 
        # copy snakefile, rules, config, and profile to working dir
        smk = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "Pipelines",
            pipeline,
            snake_n
        )
        
        rules = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "Pipelines",
            pipeline,
            rules
        )
        
        config = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "Pipelines",
            pipeline,
            config_n
        )
    
        profile = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "Pipelines",
            pipeline,
            profile_n
        )
        
        src = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "Pipelines",
            pipeline,
            "src"
        )

        # modify config file
        with open(config) as f:
            doc = yaml.safe_load(f)
        # update sif location
        doc['sif']      = sif
        doc['sort_tmp'] = os.path.join(outdir,'.tmp') 
        doc['alias']    = alias
        doc['bucket']   = bucket
        # dump
        with open(os.path.join(v["work_dir"],config_n),'w') as out:
            yaml.dump(doc,out,sort_keys=False)
   
        input_names = [snake_n,"rules",profile_n,"src"]
        dst_files = [os.path.join(v["work_dir"],i) for i in input_names]
        
        for i in zip([smk,rules,profile,src],dst_files):
            if not os.path.exists(i[1]):
                if os.path.isfile(i[0]):
                    shutil.copy(i[0],i[1])
                else:
                    shutil.copytree(i[0],i[1])
                    
        # slurm destination
        job_name = snake_n.split('.')[0]
        slurm = os.path.join(v["work_dir"],f"{v['breed']}_{k}.{job_name}.slurm")
        
        # SBATCH directives 
        header = (
            "#!/bin/bash -l\n"
            "#SBATCH -t 60:00:00\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks-per-node=1\n"
            "#SBATCH --cpus-per-task=1\n"
            "#SBATCH --mem=12gb\n"
            "#SBATCH --mail-type=ALL\n"
            f"#SBATCH --mail-user={os.environ['USER']}@umn.edu\n"
            f"#SBATCH --job-name {v['breed']}_{k}.{job_name}.slurm\n"
            f"#SBATCH -o slurm_logs/%j.{v['breed']}_{k}.{job_name}.out\n"
            f"#SBATCH -e slurm_logs/%j.{v['breed']}_{k}.{job_name}.err\n"
        )             
     
        # slurm submission body 
        with open(slurm, "w") as f:
            print(header, file=f)
            print("set -e\n",file=f)
            print("conda activate snake532",file=f)
            print("cd $SLURM_SUBMIT_DIR\n",file=f)

            print(f"FQ_DIR={fq_dir}",file=f)
            print(f"PROC_DIR={outdir}\n",file=f)

            print("# extract reference dict from container",end="",file=f)
            print(
                textwrap.dedent(
                    f"""
                    singularity exec --bind $PWD {sif} \\
                        cp /home/refgen/dog/canfam{ref[-1]}/canFam{ref[-1]}.dict $PWD
                    """
                ),file=f
            ) 

            print("# begin pipeline",file=f)
            print("start=`date +%s`",file=f)
            print(
                textwrap.dedent(
                    f"""
                    snakemake -s {snake_n} \\
                        --use-singularity \\
                        --singularity-args "-B $PWD,$FQ_DIR,$PROC_DIR" \\
                        --profile {profile_n} \\
                        --configfile {config_n} \\
                        --keep-going\n
                    """
                ),file=f
            )
            print("end=`date +%s`",file=f)
            print("echo Runtime \(seconds\): $((end-start))\n",file=f)

            print("# save slurm err/out logs",end="",file=f)
            print(
                textwrap.dedent(
                    f"""
                    mc cp --recursive ./slurm_logs/ \\
                        s3Fried/{bucket}/wgs/{breed}/{sample_name}/canfam4/
                    """
                ),file=f
            )

    plural = "sample"
    if len(d) > 1:
        plural = "samples"
    print(f"{len(d)} {plural} setup for processing")
        




# dog_meta = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/dog_ids_convert.canfam3.csv"
# outdir = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_OUT"
# fq_dir = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/Fastq"
# fq_paths = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/fullpath_fq.canfam3.list"


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
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-m", "--meta",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="csv of meta data with fastq to UMN ID conversions"
    )
    required.add_argument(
        "-f", "--fastqs",
        default=argparse.SUPPRESS,
        metavar="\b",
        required=True,
        help="path to fastq directory"
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
        help="bucket name"
    )
    optional.add_argument(
        "-r", "--ref",
        nargs="?",
        const="canfam4",
        default="canfam4",
        choices=refs,
        help="select canfam reference to use - "+ \
            " or ".join(refs) + " [default: canfam4]",
        metavar=""
    )
    optional.add_argument(
        "--sif",
        default=os.path.join(os.path.expanduser("~"),".sif/wags.sif"),
        help="location of container image [default: ~/.sif/wags.sif]"
    )
    optional.add_argument(
        "--alias",
        default="s3",
        help="minio client S3 storage alias [default: s3]"
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

    args = parser.parse_args()
    dog_meta = os.path.abspath(args.meta)
    fq_dir   = os.path.abspath(args.fastqs)
    outdir   = os.path.abspath(args.out)
    bucket   = args.bucket
    sif      = args.sif
    ref      = args.ref.lower()
    alias    = args.alias
    money    = args.money
    
    # if non default sif location
    if sif != os.path.join(os.path.expanduser("~"),".sif/wags.sif"):
        if "~" in sif:
            sif = os.path.expanduser(sif)
        else:
            sif = os.path.abspath(sif)
    
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

    # sync job
   #sync = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/SyncDogs/sync.slurm"

    # check if scratch dir exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} created!")

        
    main(dog_meta,outdir,fq_dir,ref)    
        
        