import pandas as pd
import os

#######################################
# local with BQSR #####################
#######################################

localrules: hc_intervals,
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
        # gvcf
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz",
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

# rules to include based on user setup
include: "rules/qc.smk"
include: "rules/bam.smk"
include: "rules/gvcf.smk"

