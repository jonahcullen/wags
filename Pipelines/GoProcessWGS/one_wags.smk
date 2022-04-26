import pandas as pd
import os

localrules: hc_intervals,
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

units = pd.read_table(config['units'],dtype=str).set_index('readgroup_name',drop=False)

# get breed and sample name from units
breed = units['breed'].values[0]
sample_name = units['sample_name'].values[0]

sequence_grouping(config['bucket'],config['ref_dict'])
# get sequence group intervals without unmapped, with unmapped, and hc caller intervals
intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/no_unmap","{interval}.tsv"))
unmap_intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/with_unmap","{interval}.tsv"))

rule all:
    input:
       # apply bqsr - DUE TO THE DIFFERENCE IN INTERVALS WITH OR WITHOUT UNMAPPED,
       # THIS IS NEEDED IN ADDITION TO THE FINAL MERGEGVCFS...FOR NOW...
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
            interval=unmap_intervals,
        ),
       # save fastqs
        expand(
            "{bucket}/fastqc/{breed}_{sample_name}/{u.readgroup_name}.upload",
            u=units.itertuples(), 
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
        ),
       # multiqc and save
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
            
        ),

# rules to include based on config
include: "rules/qc.smk"

if config['bqsr']:
    include: "rules/bam.smk"
else:
    include: "rules/bam.no_bqsr"

include: "rules/gvcf.smk"
include: "rules/save.smk"

