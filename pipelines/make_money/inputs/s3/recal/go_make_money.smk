import os
import textwrap
import pandas as pd

#######################################
# S3 with BQSR ########################
#######################################

localrules: manifest_and_archive,
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

# get sequence group intervals without unmapped, with unmapped, and hc caller intervals
sequence_grouping(config['bucket'],config['ref_dict'])
intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/no_unmap","{interval}.tsv"))
unmap_intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/with_unmap","{interval}.tsv"))

rule all:
    input:
        # apply bqsr - DUE TO THE DIFFERENCE IN INTERVALS WITH OR WITHOUT UNMAPPED,
        # THIS IS NEEDED IN ADDITION TO THE FINAL MERGEGVCFS...FOR NOW...
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam"
                if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{interval}.left_aligned.duplicates_marked.recalibrated.bam",
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
            interval=unmap_intervals,
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
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"],
       #    
       #),
        # upload fastqs
        expand(
            "{bucket}/fastqc/{breed}_{sample_name}/{u.readgroup_name}.upload",
            u=units.itertuples(), 
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
        ),
        # upload pipeline, and logs
        expand(
            "{bucket}/{breed}_{sample_name}.{ref}.done",
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
        ),

# rules to include based on user setup
include: "rules/bam.smk"
include: "rules/gvcf.smk"
include: "rules/genotype.smk"
include: "rules/fltr_sites_vcf.smk"
include: "rules/recal.smk"
include: "rules/gather_vep.smk"
include: "rules/compare.smk"
include: "rules/final_output.smk"
include: "rules/qc.smk"
include: "rules/save.smk"

