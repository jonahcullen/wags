import pandas as pd
from pathlib import Path
import os

localrules: input_list,
            scatter_intervals,
            generate_intervals,
           #plot_interval_lengths,
           #sites_only_gather_vcf

#configfile: "canfam4_config.yaml"
singularity: config['sif']

include: "src/get_gvcfs.py"
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

units = pd.read_csv("gvcfs.list",sep="\t")


rule all:
    input:
       #"friedlab/wgs/gvcfs/boxr/D00151/canfam4/gvcf/D00151.canfam4.done",
       #"friedlab/wgs/gvcfs/boxr/D00220/canfam4/gvcf/D00220.canfam4.done",
       #expand(
       #    "friedlab/wgs/{u.breed}/{u.dogid}/{ref}/gvcf/{u.dogid}.{ref}.done",
       #   #zip, breed=units['breed'].tolist(), dogid=units['dogid']
       #    u=units.itertuples(),
       #    ref=config['ref']
       #)
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/inputs.list",
       #    bucket = config['bucket'],
       #    ref = config['ref'],
       #    date = config['date']
       #)
       #expand("{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list",
       #    bucket = config['bucket'],
       #    ref = config['ref'],
       #    date = config['date']
       #)
       # interval plot - NOT WORKING - NEED TO FIX !!!!
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths.tiff",
       #        bucket=config['bucket'],
       #        ref=config['ref'],
       #        date=config['date']
       #    )
       #),
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/intervals/collapsed_lengths.csv",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date']
       #),
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz",
       #    bucket = config['bucket'],
       #    ref = config['ref'],
       #    date = config['date'],
       #    interval = [str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
       #),
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
       #    bucket = config['bucket'],
       #    ref = config['ref'],
       #    date = config['date'],
       #)
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/{vtype}.tranches",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #    vtype=["indels","snps"]
       #)
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/wags_{interval}/recal.{interval}.vcf.gz",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #   #interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
       #    interval=[str(i).zfill(4) for i in range(0,1+1)]
       #),
       # ONE FINAL OUTPUT - KEEP
        S3.remote(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz",
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
               #interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
            )
        ),
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/recal.{interval}.vep.vcf.gz",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #    interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
       #   #interval=[str(i).zfill(4) for i in range(0,1+1)]
       #)
       # TWO FINAL OUTPUT - KEEP
        S3.remote(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vep.vcf.gz",
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
               #interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
               #interval=[str(i).zfill(4) for i in range(0,1+1)]
            ),keep_local=True
        )
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html",
       #        bucket=config['bucket'],
       #        ref=config['ref'],
       #        date=config['date'],
       #    ),keep_local=True
       #)


include: "rules/inputs.smk"
include: "rules/intervals.smk"
include: "rules/qc.smk"
include: "rules/genotype.smk"
include: "rules/fltr_sites_vcf.smk"
include: "rules/recal.smk"
include: "rules/gather_vep.smk"

