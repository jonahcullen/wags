import os
import pandas as pd
from pathlib import Path

localrules: input_list,
            scatter_intervals,
            generate_intervals,

singularity: config['sif']
include: "src/utils.py"
# below is generate tsv of sample, breed, gvcf from s3 based on wags' fastq
# processing pipeline - will not work for local...should be rethought...
#include: "src/get_gvcfs.py"

units = pd.read_csv("gvcfs.list",sep="\t")


rule all:
    input:
        # vep'ed vcf
        expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vep.vcf.gz",
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
        # variant table
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


include: "rules/inputs.smk"
include: "rules/intervals.smk"
include: "rules/genotype.smk"
include: "rules/fltr_sites_vcf.smk"
include: "rules/recal.smk"
include: "rules/gather_vep.smk"
include: "rules/qc.smk"
include: "rules/table.smk"
include: "rules/phasing.smk"

