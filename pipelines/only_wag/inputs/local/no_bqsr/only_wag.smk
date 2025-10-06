import re
import os
import pandas as pd

#######################################
# local without BQSR ##################
#######################################

localrules: sv_done,
            multiqc,
            upload_fastqs,
            upload_pipe_and_logs

singularity: config['sif']
include: "src/utils.py"

# read inputs
units = pd.read_table(config['units'],dtype=str).set_index('readgroup_name',drop=False)

# get breed and sample name from units
breed = units['breed'].values[0]
sample_name = units['sample_name'].values[0]

# get sequence group intervals with unmapped and hc caller intervals
sequence_grouping(config['bucket'],config['ref_dict'])
intervals, = glob_wildcards(os.path.join(
    config['bucket'], "seq_group", "with_unmap", "{interval}.tsv"
))

rule all:
    input:
        # gvcf
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
        ),
        # structural variants
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/sv.done",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
        ),
        # multiqc
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
            
        ),
        #money archive
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/{breed}_{sample_name}.{ref}.MM.tar.gz",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
	)

# rules to include based on user setup
include: "rules/qc.smk"
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

