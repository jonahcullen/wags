import re
import os
import pysftp
import pandas as pd

#######################################
# SFTP without BQSR ###################
#######################################

localrules: multiqc,
            upload_fastqs,
            upload_pipe_and_logs

singularity: config['sif']
include: "src/utils.py"

from snakemake.remote.SFTP import RemoteProvider
server_user = os.environ.get('SERVER_USER')
server_pass = os.environ.get('SERVER_PASS')

cnopts = pysftp.CnOpts()
cnopts.hostkeys = None

SFTP = RemoteProvider(
    username=server_user,
    password=server_pass,
    cnopts=cnopts,
    mkdir_remote=True
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
beds, = glob_wildcards(os.path.join(f"{config['bucket']}/bed_group/","{bed}.bed"))
# NOTE THE INTERVALS ARE NOT CORRECT HERE
#
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

