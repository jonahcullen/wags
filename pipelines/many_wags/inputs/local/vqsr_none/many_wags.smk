import os
import pandas as pd
from pathlib import Path


singularity: config['sif']
include: "src/utils.py"

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
        # variant table
        expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{ref}.{date}.table",
            bucket=config['bucket'],
            ref=config['ref'],
            date=config['date'],
        ),
       ## phasing
       #expand(
       #    "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz",
       #    bucket=config['bucket'],
       #    ref=config['ref'],
       #    date=config['date'],
       #)


include: "rules/intervals.smk"
include: "rules/genotype.smk"
include: "rules/gather_and_select.smk"
include: "rules/hardfltr.smk"
include: "rules/gather_and_vep.smk"

qc_rules = "rules/qc.smk" if config['coverage_sites'] else 'rules/qc.no_cov.smk'
include: qc_rules

include: "rules/table.smk"
include: "rules/phasing.smk"

