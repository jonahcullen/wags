import re
import os
import pandas as pd

#######################################
# one_wag #############################
#######################################

# check if sv calling
sv_enabled = config.get("sv_call", False)

# check if using s3
using_s3 = config.get('remote', 'local').lower() == 's3'

# local rules
localrules: multiqc, sv_done, upload_fastqs, upload_pipe_and_logs

singularity: config['sif']
include: "src/utils.py"

if using_s3:
    # check for s3 credentials
    s3_key_id = os.environ.get('AWS_ACCESS_KEY')
    s3_access_key = os.environ.get('AWS_SECRET_KEY')

    if not s3_key_id:
        raise ValueError(
            "ERROR: S3 remote storage requested (--remote s3) but AWS_ACCESS_KEY variable not set. "
            "Set and export AWS_ACCESS_KEY or use --remote local"
        )
    if not s3_access_key:
        raise ValueError(
            "ERROR: S3 remote storage requested (--remote s3) but AWS_SECRET_KEY variable not set. "
            "Set and export AWS_SECRET_KEY or use --remote local"
        )

    try:
        from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
        
        S3 = S3RemoteProvider(
            endpoint_url='https://s3.msi.umn.edu',
            access_key_id=s3_key_id,
            secret_access_key=s3_access_key
        )
    except ImportError as e:
        raise ImportError(
            "ERROR: S3 remote storage requested but snakemake S3 provider could not be imported. "
            "Please install required dependencies or use --remote local. Error: {}".format(e)
        )

# read inputs
units = pd.read_table(config['units'],dtype=str).set_index('readgroup_name',drop=False)

# get breed and sample name from units
breed = units['breed'].values[0]
sample_name = units['sample_name'].values[0]

# get sequence group intervals with unmapped and hc caller intervals
sequence_grouping(config['bucket'],config['ref_dict'])
intervals, = glob_wildcards(os.path.join(config['bucket'], "seq_group", "with_unmap", "{interval}.tsv"))
beds, = glob_wildcards(os.path.join(config['bucket'], "bed_group", "{bed}.bed"))

core_targets = [
    # gvcf
    expand(
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz",
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
        ref=config['ref'],
    )
]

# add sv targets if requested at prep_subs
sv_targets = []
if sv_enabled:
    sv_targets = expand(
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/sv.done",
        bucket=config['bucket'],
        breed=breed,
        sample_name=sample_name,
        ref=config['ref'],
    )

# add upload targets if using s3
upload_targets = []
if using_s3:
    upload_targets.extend([
        # upload fastqs
        expand(
            "{bucket}/fastqc/{breed}/{sample_name}/{u.readgroup_name}.upload",
            u=units.itertuples(), 
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
        ),
        # upload pipeline, and logs
        expand(
            "{bucket}/{ref}/{breed}/{sample_name}/{ref}/pipe_and_logs.upload",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
        )
    ])

# flatten all targets to single list
all_targets = []
for target in core_targets:
    all_targets.extend(target)
all_targets.extend(sv_targets)
# flatten upload targets and extend
for target in upload_targets:
    all_targets.extend(target)

rule all:
    input:
        all_targets

# rules to include based on user setup
include: "rules/qc.smk"
include: "rules/bam.smk"
include: "rules/gvcf.smk"

if sv_enabled:
    include: "rules/sv.smk"
if using_s3:
    include: "rules/save.smk"

