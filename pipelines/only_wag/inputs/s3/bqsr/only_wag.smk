import os
import re
import textwrap
import pandas as pd

#######################################
# S3 with BQSR ########################
#######################################

localrules: sv_done,
            manifest_and_archive,
            multiqc,
            upload_fastqs,
            upload_pipe_and_logs

singularity: config['sif']
include: "src/utils.py"

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id,
    secret_access_key=s3_access_key
)

# read inputs
units = pd.read_table(config['units'],dtype=str).set_index('readgroup_name',drop=False)

# get breed and sample name from units
breed = units['breed'].values[0]
sample_name = units['sample_name'].values[0]

# get sequence group intervals with unmapped and hc caller intervals
sequence_grouping(config['bucket'],config['ref_dict'])
intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/with_unmap","{interval}.tsv"))

rule all:
    input:
        # structural variants
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/sv.done",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
        ),
        # money archive
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/{breed}_{sample_name}.{ref}.K9MM.tar.gz",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
        ),
        # multiqc
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
            
        ),
        # upload fastqs
        expand(
            "{bucket}/fastqc/{breed}/{sample_name}/{u.readgroup_name}.upload",
            u=units.itertuples(), 
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
        ),
        # upload pipeline, and logs
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/upload_pipe.done",
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
        ),

# rules to include based on user setup
include: "rules/bam.smk"
include: "rules/sv.smk"
include: "rules/gvcf.smk"
include: "rules/genotype.smk"
include: "rules/gather_and_select.smk"
include: "rules/hardfltr.smk"
include: "rules/gather_and_vep.smk"
include: "rules/common.smk"
include: "rules/compare.smk"
include: "rules/final_output.smk"
include: "rules/qc.smk"
include: "rules/save.smk"

