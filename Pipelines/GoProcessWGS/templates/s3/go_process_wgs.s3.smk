import pandas as pd
import os

localrules: hc_intervals,
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

units = pd.read_table(config['units'],dtype=str).set_index('readgroup_name',drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

# get breed and sample name from units
breed = units['breed'].values[0]
#breed = units["breed"].unique()[0]
sample_name = units['sample_name'].values[0]

sequence_grouping(config['bucket'],config['ref_dict'])
# get sequence group intervals without unmapped, with unmapped, and hc caller intervals
intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/no_unmap","{interval}.tsv"))
unmap_intervals, = glob_wildcards(os.path.join(f"{config['bucket']}/seq_group/with_unmap","{interval}.tsv"))
#hc_intervals, = glob_wildcards(os.path.join(config['hc_intervals'],"{hc_interval}.interval_list"))

rule all:
    input:
       ## fastqc - FOR NOW JUST USING mino client in shell command due to 
       ## fastqc's output naming conventions
       #expand(
       #   #"{bucket}/fastqc/{u.sample_name}/{u.flowcell}/{u.readgroup_name}/qc.done",
       #    "{bucket}/wgs/{breed}/{u.sample_name}/fastqc/{u.sample_name}/{u.flowcell}/{u.readgroup_name}_R2_fastqc.zip",
       #    u=units.itertuples(), 
       #    bucket=config['bucket'],
       #    breed=breed
       #),
       #S3.remote(expand("{bucket}/wgs/{breed}/{u.sample_name}/fastqc/{u.readgroup_name}/"
       # fastp
       #S3.remote(
       #    expand(
       #        '{bucket}/private/fastq/trimmed/{u.tissue}/{u.tissue}_{u.sample}_{u.lane}_{read}_001.fastq.gz',
       #        u=units.itertuples(), bucket=config['bucket'], 
       #        read=['R1','R2']
       #    )
       #)
       ## fastq to ubam
       #expand(
       #    "results/fastqs_to_ubam/{u.sample_name}/{u.readgroup_name}.unmapped.bam",
       #    u=units.itertuples()
       #)
       # mark adapters
       #expand(
       #    "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.readgroup_name}.mark_adapt.unmapped.bam",
       #   #"results/fastqs_to_ubam/{u.sample_name}/{u.readgroup_name}.mark_adapt.metrics.txt"
       #    u=units.itertuples(),
       #    bucket=config['bucket'],
       #    breed=breed,
       #    ref=config['ref'],
       #)
       ## sam_to_fastq_and_bwa_mem
       #expand(
       #   #"{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.readgroup_name}.mark_adapt.unmapped.bam",
       #   #"results/sam_to_fastq_and_bwa_mem/{u.sample_name}/{u.readgroup_name}.{ref}_aligned.unmerged.bam",
       #    "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.readgroup_name}.{ref}_aligned.unmerged.bam",
       #    u=units.itertuples(),
       #    bucket=config['bucket'],
       #    breed=breed,
       #    ref=config['ref'],
       #)
       ## merge bams
       #expand(
       #   #"{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.readgroup_name}.mark_adapt.unmapped.bam",
       #   #"results/sam_to_fastq_and_bwa_mem/{u.sample_name}/{u.readgroup_name}.{ref}_aligned.unmerged.bam",
       #    "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.readgroup_name}.{ref}.merged.unsorted.bam",
       #    u=units.itertuples(),
       #    bucket=config['bucket'],
       #    breed=breed,
       #    ref=config['ref'],
       #)
       #expand("results/merge_bam_alignment/{u.sample_name}/{u.readgroup_name}.{ref}.merged.unsorted.bam",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## mark duplicates
       #expand(
       #   #"results/mark_duplicates/{u.sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config['ref'],
       #),
       ## sort and fix tags
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
       #   #"results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config['ref'],
       #),
       ## base recalibrator
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{interval}.recal_data.csv",
       #   #"results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config['ref'],
       #    interval=intervals
       #),
       ## gather bqsr reports
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.recal_data.csv",
       #   #"results/gather_bqsr_reports/{u.sample_name}.{ref}.recal_data.csv",
       #    bucket=config["bucket"],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"]
       #),
       # apply bqsr - DUE TO THE DIFFERENCE IN INTERVALS WITH OR WITHOUT UNMAPPED,
       # THIS IS NEEDED IN ADDITION TO THE FINAL MERGEGVCFS...FOR NOW...
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
            interval=unmap_intervals,
        ),
       ## gather bams
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.{ext}",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"],
       #    ext=["bam","bai","bam.md5"]
       #),
       # doc and flagstat - KEEP
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.{ext}",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"],
       #    ext=["flagstat"]
       #),
       # analyze covariates
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.analyze_cov.{ext}",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"],
       #    ext=["pdf","csv"]
       #),
       # haplotype caller
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{sample_name}.00{split}.g.vcf.gz",
       #   #"{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{hc_interval}/{sample_name}.{ref}.g.vcf.gz",
       #   #"results/haplotype_caller/{hc_interval}/{u.sample_name}.{ref}.g.vcf.gz",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config["ref"],
       #    split=list(map("{:02d}".format, list(range(0,config['scatter_size']))))
       #),
       ## merge gvcfs
       #expand(
       #    "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.{ext}",
       #    bucket=config['bucket'],
       #    breed=breed,
       #    sample_name=sample_name,
       #    ref=config['ref'],
       #    ext=["gz","gz.tbi"]
       #),
       # save fastqs
        expand(
            "{bucket}/fastqc/{breed}_{sample_name}/{u.readgroup_name}.upload",
            u=units.itertuples(), 
            bucket=config["bucket"],
            breed=breed,
            sample_name=sample_name,
        ),
       # multiqc and save
        expand(
            "{bucket}/{breed}_{sample_name}_{ref}.done",
           #"{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
            
        ),

include: "rules/qc.smk"
include: "rules/bam.smk"
include: "rules/gvcf.smk"
include: "rules/save.smk"



