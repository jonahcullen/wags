import os
import pandas as pd
from pathlib import Path

localrules: scatter_intervals,
            generate_intervals

singularity: config['sif']
include: "src/utils.py"

# setup snakemake s3 remote provider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id,
    secret_access_key=s3_access_key
)

rule all:
    input:
        # vep'ed vcf
        expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
            bucket=config['bucket'],
            ref=config['ref'],
            date=config['date'],
        ),
        # multiqc report
        expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html",
            bucket=config['bucket'],
            ref=config['ref'],
            date=config['date'],
        ),
       ## variant table
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{date}.{ref}.table",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #),
       ## phasing
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_genotype.{ref}.snps.phased.vcf.gz",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #)


include: "rules/intervals.smk"
include: "rules/genotype.smk"
include: "rules/fltr_sites_vcf.smk"
include: "rules/recal.smk"
include: "rules/gather_vep.smk"
include: "rules/qc.smk"
include: "rules/table.smk"
include: "rules/phasing.smk"
