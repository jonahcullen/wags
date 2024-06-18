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

def main():
    global species
    # if using a known ref - load and modify.
    # ref_fasta, ref_dict, and ref_gtf are already
    # avaiable within container and do not need
    # to be specified
    prep_dir = os.path.dirname(os.path.realpath(__file__))
    if ref in refs:
        config_n  = f"{ref}_config.yaml"
        config = os.path.join(
            str(Path(prep_dir).parents[0]),
            "pipelines/many_wags/configs",
            config_n
        )
        # load known config
        with open(config) as f:
            doc = yaml.safe_load(f)
        
        # get config values
        species  = doc['species']
    else:
        # build config from custom_config.yaml
        config_n = f"{ref}_config.yaml"
        config = os.path.join(
            str(Path(prep_dir).parents[0]),
            "pipelines/many_wags/configs/custom_config.yaml",
        )
        # load known config
        with open(config) as f:
            doc = yaml.safe_load(f)
        # update config with non-known ref 
        doc['species']   = species
        doc['ref']       = ref
        doc['ref_fasta'] = fasta
        doc['ref_dict']  = fasta_dict
        doc['ref_gtf']   = gtf
    
    # set directory for config file
    if not args.out:
        outdir = os.path.join(os.path.expanduser("~"), f".wags/{species}/{ref}/many_wags")
    else:
        outdir = os.path.realpath(os.path.expanduser(args.out))
    os.makedirs(outdir, exist_ok=True)

    # update known ref config with user cli args
   #doc['date']               = date
    doc['bucket']             = bucket
    doc['alias']              = alias
    doc['sif']                = sif
    doc['profile']            = profile
    doc['snp_filter_level']   = snp_val
    doc['indel_filter_level'] = indel_val
    if phase is not None:
        doc['phasing'] = true
        doc['link_map'] = phase

    # check if exists and dump
    config_out = os.path.join(outdir, config_n)
    if os.path.isfile(config_out):
        ow = input(f"config already exists at {config_out} - overwrite? [y/n] ")
        if ow.lower() == 'y':
            with open(config_out,'w') as f:
                yaml.dump(doc, f, sort_keys=False)
                print(f"{ref}_config.yaml saved to {outdir}")
        else:
            print("nothing to do - confirm ref name correct or try again")
    else:
        with open(config_out,'w') as f:
            yaml.dump(doc, f, sort_keys=False)
            print(f"{ref}_config.yaml saved to {outdir}")

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="wags",
        add_help=False,
        description=(
            "script to generate config file for joint genotpying \n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-r", "--ref",
       #default="canfam4",
       #metavar="",
        default=argparse.SUPPRESS,
        metavar="\b",
        help=textwrap.dedent(f'''\
            select reference to use: {", ".join(refs)}.
            included references have --fasta, --dict, and
            --gtf already included and do not need to be 
            set as arguments if desired. if using different
            reference, set --fasta, --fasta_dict, and
            --gtf arguments
        ''')
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
    optional.add_argument(
        "-o", "--out",
        metavar="\b",
        help="location of config out [default: ~/.wags/SPECIES/REF/many_wags/REF_config.yaml]"
    )
    optional.add_argument(
        "-s", "--species",
        metavar="\b",
        help="species to be joint genotyped"
    )
    optional.add_argument(
        "-f", "--fasta",
        metavar="\b",
        help="(if custom ref) path to reference fasta"
    )
    optional.add_argument(
        "-d", "--fasta_dict",
        metavar="\b",
        help="(if custom ref) path to reference dict"
    )
    optional.add_argument(
        "-g", "--gtf",
        metavar="\b",
        help="(if custom ref) path to reference gtf"
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
        "--snp-tranche-level",
        default="99.0",
        help="truth sensitivity level for snp filtering (vqsr only) [default: 99.0]",
    )
    optional.add_argument(
        "--indel-tranche-level",
        default="99.0",
        help="truth sensitivity level for snp filtering (vqsr only) [default: 99.0]",
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

    args = parser.parse_args()
    ref        = args.ref
    snp_val    = args.snp_tranche_level
    indel_val  = args.indel_tranche_level
    bucket     = args.bucket
    outdir     = args.out
    species    = args.species
    alias      = args.alias
    sif        = args.sif
    profile    = args.profile
    remote     = args.remote.lower()
    fasta      = args.fasta
    fasta_dict = args.fasta_dict
    gtf        = args.gtf
    phase      = args.phase
    
    # require species argument if not known
    if ref not in refs: 
        if species is None:
            parser.error(f"--species (-s) required if using ref other than {refs}")
 
    # if non default sif location
    if sif != os.path.join(os.path.expanduser("~"),".sif/wags.sif"):
        sif = os.path.realpath(os.path.expanduser(sif))
    
    # confirm image exitst
   #if os.path.isfile(sif):
   #    print("wags image found!")
   #else:
   #    sif_dir = os.path.join(os.path.expanduser("~"),".sif")
   #    print(
   #        f"wags image not found at {os.path.realpath(sif)} -> downloading to {sif_dir}"
   #    )
   #    url = "https://s3.msi.umn.edu/wags/wags.sif"
   #    os.makedirs(sif_dir,exist_ok=True)
   #    wget.download(url,sif_dir)
   #    sys.exit("\nrerun without --sif or point to directory containing wags.sif")

    main()
        
        
