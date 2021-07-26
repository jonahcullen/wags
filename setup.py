#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 20:43:49 2021

@author: jonahcullen
"""


import argparse
import os
import csv
from collections import defaultdict
import gzip
from datetime import datetime
import shutil
import string
import textwrap
import pandas as pd


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
                    if f.startswith(line[-1]) and "_R2" in f:
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
            
            d[line[0]]["work_dir"] = os.path.join(outdir,line[1],line[0],ref)
            d[line[0]]["breed"] = line[1]
            d[line[0]]["df"] = pd.DataFrame(dog_input, columns=cols)
            
            
    # copy pipeline input and slurm file
    for k,v in d.items():
        os.makedirs(v["work_dir"], exist_ok=True)

        # create jobs dir if not exist
        jobs = os.path.join(v["work_dir"],"Jobs")
        if not os.path.exists(jobs):
            os.makedirs(jobs)
   
        # write dog input to working directory
        with open(os.path.join(v["work_dir"],"input.tsv"),"w") as out:
            v["df"].to_csv(out,sep='\t',index=False)
       
        # input templates
        pipeline  = "GoProcessWGS"
        snake_n   = "go_process_wgs.smk"
        config_n  = f"{ref}_config.yaml"
        profile_n = "slurm.go_wgs"
        # switch to money templates if true
        if money:
            pipeline  = "GoMakeMoney"
            snake_n   = "go_make_money.smk"
            config_n  = f"{ref}_money.yaml"
            profile_n = "slurm.go_money"
 
        # copy snakefile, config, and profile to working dir
        smk = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "Pipelines",
                               pipeline,
                               snake_n)
        
        config = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "Pipelines",
                                  pipeline,
                                  config_n)
    
        profile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "Pipelines",
                                  pipeline,
                                  profile_n)
        
        src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "Pipelines",
                                  pipeline,
                                  "src")
        
       #input_names = ["go_process_wgs.smk",f"{ref}_config.yaml","slurm.go_wgs","src"]
        input_names = [snake_n,config_n,profile_n,"src"]
        dst_files = [os.path.join(v["work_dir"],i) for i in input_names]
        
        for i in zip([smk,config,profile,src],dst_files):
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
            "#SBATCH --mem=60gb\n"
            "#SBATCH --mail-type=ALL\n"
            f"#SBATCH --mail-user={os.environ['USER']}@umn.edu\n"
            f"#SBATCH --job-name {v['breed']}_{k}.{job_name}.slurm\n"
            f"#SBATCH -o Jobs/%j.{v['breed']}_{k}.{job_name}.out\n"
            f"#SBATCH -e Jobs/%j.{v['breed']}_{k}.{job_name}.err\n"
        )             
     
        # slurm submission body 
        with open(slurm, "w") as f:
            print(header, file=f)
            print("set -e\n",file=f)
            print("conda activate snake532",file=f)
            print("cd $SLURM_SUBMIT_DIR\n",file=f)

            print(
                textwrap.dedent(
                    f"""
                    snakemake -s {snake_n} \\
                        --profile {profile_n} \\
                        --configfile {config_n} \\
                        --keep-going
                    """
                ),
            file=f)

           #ref_dog_breed = os.path.join(base,ref.lower(),v["breed"],k)

           #print(f"mkdir -p {os.path.join(ref_dog_breed,'jobs')}\n" +
           #      f"cp -t {os.path.join(ref_dog_breed,'jobs')} *.{{slurm,smk,yaml,err,out}}\n",
           #      file=f)
           #
           #print(f"sbatch  --mail-user={os.environ['USER']}@umn.edu \\\n" +
           #      f'    --export=BREED={v["breed"]},DOG_ID={k},REF={ref.lower()},DOG_DIR="$SLURM_SUBMIT_DIR" \\\n' +
           #      f'    {sync}\n\nwait\n',file=f)
           #
           #print(
           #    textwrap.dedent(
           #        f"""
           #        sbatch --mail-user={os.environ['USER']}@umn.edu \\
           #            --export=BREED={v["breed"]},DOG_ID={k},REF={ref.lower()},DOG_DIR="$SLURM_SUBMIT_DIR" \\
           #            {sync}\n\nwait
           #        """
           #    ),
           #file=f)

    print(f"{len(d)} dogs setup for processing")
        




# dog_meta = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/dog_ids_convert.canfam3.csv"
# outdir = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_OUT"
# fq_dir = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/Fastq"
# fq_paths = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/RESCUER/TEST_DATA/fullpath_fq.canfam3.list"


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="Dogpile",
        add_help=False,
        description=(
            "Dogpile generates all required input to process FASTQs to gVCF \n"
            "following GATK best practices. For each sample (and associated \n"
            "FASTQ pair), Dogpile outputs a directory structure organized by \n"
            "breed, wherein GATK pipeline input are contained by sample ID."
        )
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
    fq_dir = os.path.abspath(args.fastqs)
    outdir = os.path.abspath(args.out)
    ref = args.ref.lower()
    money = args.money

    # base dir to save pipeline output in primary
    base = "/panfs/roc/groups/0/fried255/fried255/working"

    # sync job
   #sync = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/SyncDogs/sync.slurm"

    # check if scratch dir exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"{outdir} created!")

        
    main(dog_meta,outdir,fq_dir,ref)    
        
        
