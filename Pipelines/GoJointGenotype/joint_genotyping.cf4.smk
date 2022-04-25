import pandas as pd
from pathlib import Path
import os

localrules: input_list,
            scatter_intervals,
            generate_intervals,

#configfile: "canfam4_config.yaml"
singularity: config['sif']


# setup snakemake s3 remote provider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id,
    secret_access_key=s3_access_key
)

units = pd.read_csv("gvcfs.list",sep="\t")


rule all:
    input:
        # vep'ed vcf
        S3.remote(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vep.vcf.gz",
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
            ),keep_local=True
        ),
        # multiqc reprot
        S3.remote(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html",
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
            ),keep_local=True
        ),
        # variant table
        S3.remote(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{date}.{ref}.table",
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
            ),keep_local=True
        )
        # phasing - no CF4 yet


include: "src/get_gvcfs.py"
include: "src/utils.py"
include: "rules/inputs.smk"
include: "rules/intervals.smk"
include: "rules/qc.smk"
include: "rules/genotype.smk"
include: "rules/fltr_sites_vcf.smk"
include: "rules/recal.smk"
include: "rules/gather_vep.smk"
include: "rules/table.smk"


