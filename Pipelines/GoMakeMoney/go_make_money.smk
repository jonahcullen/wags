# CanFam3 money

import pandas as pd
import os

localrules: manifest_and_upload,
            save_jobs

include: "src/utils.py"

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id,
    secret_access_key=s3_access_key
)

units = pd.read_table(config["units"],dtype=str).set_index("readgroup_name",drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

# get breed from units
breed = units["breed"].unique()[0]

sequence_grouping(config["ref_dict"])
# get sequence group intervals without unmapped, with unmapped, and hc caller intervals
intervals, = glob_wildcards(os.path.join("results/seq_group/no_unmap","{interval}.tsv"))
unmap_intervals, = glob_wildcards(os.path.join("results/seq_group/with_unmap","{interval}.tsv"))
hc_intervals, = glob_wildcards(os.path.join(config["hc_intervals"],"{hc_interval}.interval_list"))

rule all:
    input:
       ## fastqc - FOR NOW JUST USING mino client in shell command due to 
       ## fastqc's output naming conventions
        expand(
            "{bucket}/fastqc/{u.sample_name}/{u.readgroup_name}/qc.done",
            u=units.itertuples(), bucket=config["bucket"]
        ),
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
       #expand("results/fastqs_to_ubam/{u.sample_name}/{u.readgroup_name}.unmapped.bam",
       #    u=units.itertuples()
       #)
       ## sam_to_fastq_and_bwa_mem
       #expand("results/sam_to_fastq_and_bwa_mem/{u.sample_name}/{u.readgroup_name}.aligned.unmerged.bam",
       #    u=units.itertuples()
       #),
       ## merge bams
       #expand("results/merge_bam_alignment/{u.sample_name}/{u.readgroup_name}.{ref}.merged.unsorted.bam",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## mark duplicates
       #expand("results/mark_duplicates/{u.sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## sort and fix tags
       #expand("results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
       #    sample_name=units["sample_name"].values[0],ref=config["ref"]
       #),
       ## base recalibrator
       #expand("results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
       #    sample_name=units["sample_name"].values[0],ref=config["ref"],
       #    interval=intervals
       #),
       ## gather bqsr reports
       #expand("results/gather_bqsr_reports/{u.sample_name}.{ref}.recal_data.csv",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## apply bqsr - DUE TO THE DIFFERENCE IN INTERVALS WITH OR WITHOUT UNMAPPED,
       ## THIS IS NEEDED IN ADDITION TO THE FINAL MERGEGVCFS...FOR NOW...
        expand("{bucket}/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
            sample_name=units["sample_name"].values[0],ref=config["ref"],
            interval=unmap_intervals, bucket=config["bucket"],
            breed=breed
        ),
       ## gather bams
       #expand("results/gather_bam_files/{u.sample_name}.{ref}.bam",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## gather bams
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/bam/{u.sample_name}.{ref}.{ext}",
       #        u=units.itertuples(), bucket=config['bucket'],
       #        ref=config["ref"], breed=breed, 
       #        ext=["bam","bai","bam.md5"]
       #    )
       #),
       ## doc and flagstat
       #expand("results/coverage_depth_and_flagstat/{u.sample_name}.{ref}.depthofcoverage.sample_summary",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## haplotype caller
       #expand("results/haplotype_caller/{hc_interval}/{u.sample_name}.{ref}.g.vcf.gz",
       #    u=units.itertuples(),ref=config["ref"],
       #    hc_interval=hc_intervals
       #),
       ## merge gvcfs
       #expand("results/merge_gvcfs/{u.sample_name}.{ref}.g.vcf.gz",
       #    u=units.itertuples(),ref=config["ref"]
       #)
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/gvcf/{u.sample_name}.{ref}.g.vcf.gz",
       #        u=units.itertuples(), bucket=config['bucket'],
       #        ref=config["ref"], breed=breed, 
       #       #ext=["bam","bai","bam.md5"]
       #    )
       #),
       ## sites only gather 
       #"results/sites_only_gather_vcf/gather.sites_only.vcf.gz",
       ## final gather 
       #expand("results/final_gather_vcfs/joint_genotype.{ref}.{ext}",
       #    ref=config["ref"],
       #    ext=["vcf.gz","vcf.gz.tbi"]
       #)
       ## final output
       #f"results/final_output/{units['sample_name'].values[0]}_{config['ref']}_tables.xlsx",
       #expand("results/final_output/{ref}/{sample_out}_{ref}_tables.xlsx",
       #    sample_out=units['sample_name'].values[0],
       #    ref=config["ref"]
       #)
       ## archive and gzip final out
       #expand("results/final_output/{ref}/{sample_out}_{ref}_K9MM.tar.gz",
       #    sample_out=units['sample_name'].values[0],
       #    ref=config["ref"]
       #)
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_K9MM.tar.gz",
       #        sample_out=units['sample_name'].values[0], 
       #        bucket=config['bucket'], ref=config["ref"], 
       #        breed=breed, 
       #    )
       #),
        expand(
            "{bucket}/save_jobs/{breed}_{sample_out}_{ref}_jobs.done",
            sample_out=units['sample_name'].values[0], 
            bucket=config['bucket'], ref=config["ref"], 
            breed=breed, 
        )


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        touch("{bucket}/fastqc/{sample_name}/{readgroup_name}/qc.done")
       #fqs = S3.remote('{bucket}/private/fastq/trimmed/{tissue}/{tissue}_{sample}_{lane}.fastp.html'),
       #qc  =
    params:
        fastqc = config["fastqc"],
        alias  = config["alias"]
    resources:
         time   = 360,
         mem_mb = 6000, 
         cpus   = 4
    run:
        # get flowcell
        flowcell = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"flowcell"].values[0]

        # output path to fastqc by flowcell
       #qc_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "fastqc",
       #        flowcell
       #)
        qc_outdir = os.path.join(
                "results",
                "fastqc",
                config["ref"],
                flowcell
        )
        
       #qc_name = os.path.basename(input.r1).split(".")[0] + "_fastqc"
        # output path to fastq by flowcell
       #fq_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "fastq",
       #        flowcell
       #
       #)

        shell(f'''
            set -e

            mkdir -p {qc_outdir}

            {{params.fastqc}} --outdir {qc_outdir} \
                {{input.r1}} {{input.r2}}

            mcli cp --recursive {qc_outdir} \
                {{params.alias}}/{config["bucket"]}/wgs/{breed}/{wildcards.sample_name}/fastqc/{flowcell}/

            mcli cp {{input.r1}} {{input.r2}} \
                {{params.alias}}/{config["bucket"]}/wgs/{breed}/{wildcards.sample_name}/fastq/{flowcell}/
        ''')

rule fastqs_to_ubam:
    input:
        unpack(get_fastq),
       #"results/fastqc/{sample_name}/{readgroup_name}/qc.done"
    output:
        ubam = "{bucket}/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    params:
        java_opt  = "-Xms6000m",
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
       #tmp_dir   = f"/scratch.global/friedlab_{os.environ['USER']}/{{readgroup_name}}.ubam.tmp"
        tmp_dir   = f"/dev/shm/friedlab_{os.environ['USER']}/{{readgroup_name}}_{config['ref']}.ubam.tmp"
    threads: 6
    resources:
         time   = 400,
         mem_mb = 240000,
    run:

        # get fastq meta data for logging in the bam
        library_name = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"library_name"].values[0]
        platform_unit = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"platform_unit"].values[0]
        run_date = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"run_date"].values[0]
        platform_name = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"platform_name"].values[0]
        sequencing_center = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"sequencing_center"].values[0]

        shell(f'''
            set -e

	        mkdir -p {{params.tmp_dir}}

            {{params.gatk}} --java-options {{params.java_opt}} \
                FastqToSam \
                --TMP_DIR {{params.tmp_dir}} \
                --FASTQ {{input.r1}} \
                --FASTQ2 {{input.r2}} \
                --OUTPUT {{output.ubam}} \
                --READ_GROUP_NAME {{wildcards.readgroup_name}} \
                --SAMPLE_NAME {{wildcards.sample_name}} \
                --LIBRARY_NAME {library_name} \
                --PLATFORM_UNIT {platform_unit} \
                --RUN_DATE {run_date} \
                --PLATFORM {platform_name} \
                --SEQUENCING_CENTER {sequencing_center}

	        rm -rf {params.tmp_dir}
        ''')

rule sam_to_fastq_and_bwa_mem:
    input:
        ubam = "{bucket}/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    output:
        bwa_log = "{bucket}/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.{ref}_aligned.unmerged.bwa.stderr.log",
        bam     = "{bucket}/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.{ref}_aligned.unmerged.bam"
    params:
        comp_level = 5,
        java_opt   = "-Xms3000m",
        gatk       = config["gatk"],
        picard     = config["picard"],
        samtools   = config["samtools"],
        bwa        = config["bwa"],
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        ref_fasta  = config["ref_fasta"],
    threads: 16
    resources:
         time   = 1440,
         mem_mb = 24000, 
    shell:
        '''
            set -o pipefail
            set -e

            # set the bash variable needed for the command-line
            bash_ref_fasta={params.ref_fasta}

		    java -Dsamjdk.compression_level={params.comp_level} {params.java_opt} -jar {params.picard} \
                SamToFastq \
			    INPUT={input.ubam} \
			    FASTQ=/dev/stdout \
			    INTERLEAVE=true \
			    NON_PF=true \
            | \
		    {params.bwa} {params.bwa_cl} $bash_ref_fasta /dev/stdin -  2> >(tee {output.bwa_log} >&2) \
            | \
		    {params.samtools} view -1 - > {output.bam}
        '''

rule merge_bam_alignment:
    input:
        ubam = "{bucket}/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam",
        bam  = "{bucket}/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.{ref}_aligned.unmerged.bam"
    output:
        merged_bam = "{bucket}/merge_bam_alignment/{sample_name}/{readgroup_name}.{ref}.merged.unsorted.bam"
    params:
        comp_level = 5,
        java_opt   = "-Xms3000m",
        gatk       = config["gatk"],
        bwa        = config["bwa"],
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        bwa_ver    = "0.7.17-r1188",
        ref_fasta  = config["ref_fasta"]
    threads: 4
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        # DOUBLE CHECK BAM HEADER - UNCLEAR IF BWA_VER WILL EVALUATE AS EXPECTED
        # BWA_VER=$(/panfs/roc/groups/0/fried255/shared/gatk4_workflow/tools/bwa-0.7.17/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //' 2>&1)
        '''
            set -e
            
            # set the bash variable needed for the command-line
            bash_ref_fasta={params.ref_fasta}
            
            # get bwa version
            # this had to removed as it keep causing an error - unsure why...TROUBLE SHOOT NEED

            {params.gatk} --java-options "-Dsamjdk.compression_level={params.comp_level} {params.java_opt}" \
                MergeBamAlignment \
                --VALIDATION_STRINGENCY SILENT \
                --EXPECTED_ORIENTATIONS FR \
                --ATTRIBUTES_TO_RETAIN X0 \
                --ALIGNED_BAM {input.bam} \
                --UNMAPPED_BAM {input.ubam} \
                --OUTPUT {output.merged_bam} \
                --REFERENCE_SEQUENCE {params.ref_fasta} \
                --PAIRED_RUN true \
                --SORT_ORDER "unsorted" \
                --IS_BISULFITE_SEQUENCE false \
                --ALIGNED_READS_ONLY false \
                --CLIP_ADAPTERS false \
                --MAX_RECORDS_IN_RAM 2000000 \
                --ADD_MATE_CIGAR true \
                --MAX_INSERTIONS_OR_DELETIONS -1 \
                --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
                --PROGRAM_RECORD_ID "bwamem" \
                --PROGRAM_GROUP_VERSION "{params.bwa_ver}" \
                --PROGRAM_GROUP_COMMAND_LINE "{params.bwa_cl}" \
                --PROGRAM_GROUP_NAME "bwamem" \
                --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
                --ALIGNER_PROPER_PAIR_FLAGS true \
                --UNMAP_CONTAMINANT_READS true
        '''

rule mark_duplicates:
    input:
        merged_bams = expand(
                        "{bucket}/merge_bam_alignment/{sample_name}/{readgroup_name}.{ref}.merged.unsorted.bam",
                        sample_name=units["sample_name"].values[0],
                        readgroup_name=list(units["readgroup_name"]),
                        ref=config["ref"], bucket=config["bucket"]
                      )
    output:
        dedup_bam = S3.remote("{bucket}/mark_duplicates/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam"),
       #metrics   = "results/mark_duplicates/{sample_name}.{ref}.duplicate_metrics"
        metrics   = S3.remote(f"{{bucket}}/wgs/{breed}/{{sample_name}}/{{ref}}/logs/{{sample_name}}.{{ref}}.duplicate_metrics"),
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m -Xmx16g",
        gatk       = config["gatk"],
        tmp_dir    = "/dev/shm/{sample_name}_{ref}.md.tmp"
    threads: 4
    resources:
         time   = 720,
         mem_mb = 240000
    run:
        # separate bams by --INPUT
        bams = " --INPUT ".join(map(str,input.merged_bams))

        # primary logs dir
       #log_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "logs"
       #)

        shell(f'''
            mkdir -p {{params.tmp_dir}}

            {{params.gatk}} --java-options "-Dsamjdk.compression_level={{params.comp_level}} {{params.java_opt}}" \
                MarkDuplicates \
                --TMP_DIR {{params.tmp_dir}} \
                --INPUT {bams} \
                --OUTPUT {{output.dedup_bam}} \
                --METRICS_FILE {{output.metrics}} \
                --VALIDATION_STRINGENCY SILENT \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --ASSUME_SORT_ORDER "queryname" \
                --CREATE_MD5_FILE true
        ''')

rule sort_and_fix_tags:
    input:
        dedup_bam = "{bucket}/mark_duplicates/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam"
    output:
        sorted_bam = "{bucket}/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        sorted_bai = "{bucket}/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai"
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m",
        gatk       = config["gatk"],
        ref_fasta  = config["ref_fasta"],
        tmp_dir    = f"/scratch.global/friedlab_{os.environ['USER']}/{{sample_name}}_{config['ref']}.sort.tmp"
    threads: 12
    resources:
         time   = 720,
         mem_mb = 32000
    shell:
        '''
            set -o pipefail
            
            mkdir -p {params.tmp_dir}

            {params.gatk} --java-options "-Dsamjdk.compression_level={params.comp_level} {params.java_opt}" \
                SortSam \
                --TMP_DIR {params.tmp_dir} \
                --INPUT {input.dedup_bam} \
                --OUTPUT /dev/stdout \
                --SORT_ORDER "coordinate" \
                --CREATE_INDEX false \
                --CREATE_MD5_FILE false \
            | \
            {params.gatk} --java-options "-Dsamjdk.compression_level={params.comp_level} {params.java_opt}" \
                SetNmMdAndUqTags \
                --INPUT /dev/stdin \
                --OUTPUT {output.sorted_bam} \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true \
                --REFERENCE_SEQUENCE {params.ref_fasta}

            rm -rf {params.tmp_dir}
       '''

rule base_recalibrator:
    input:
        sorted_bam = "{bucket}/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        interval   = ancient("results/seq_group/no_unmap/{interval}.tsv")
    output:
        recal_csv = "{bucket}/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv"
    params:
        java_opt         = "-Xms4000m",
        gatk             = config["gatk"],
        ref_fasta        = config["ref_fasta"],
        dbsnp_snp_vcf    = config["dbsnp_snp_vcf"],
        broad_snp_vcf    = config["broad_snp_vcf"],
        axelsson_snp_vcf = config["axelsson_snp_vcf"],
        dbsnp_indels_vcf = config["dbsnp_indels_vcf"],
    resources:
         time   = 120,
         mem_mb = 20000
    run:
        # separate content of interval by -L
        with open(input.interval,"r") as f:
            ival = f.read().strip().replace("\t", " -L ")

        shell(f'''
            {{params.gatk}} --java-options {{params.java_opt}} \
                BaseRecalibrator \
                -R {{params.ref_fasta}} \
                -I {{input.sorted_bam}} \
                --use-original-qualities \
                -O {{output.recal_csv}} \
                --known-sites {{params.dbsnp_snp_vcf}} \
                --known-sites {{params.broad_snp_vcf}} \
                --known-sites {{params.axelsson_snp_vcf}} \
                --known-sites {{params.dbsnp_indels_vcf}} \
                -L {ival}
        ''')

rule gather_bqsr_reports:
    input:
       #bqsr_reports = sorted(
       #            expand("results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
       #                sample_name=units["sample_name"].values[0],
       #                ref=config["ref"],
       #                interval=intervals
       #            )
       #)
        bqsr_reports = sorted(
                    expand("{bucket}/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
                        sample_name=units["sample_name"].values[0],
                        ref=config["ref"], bucket=config["bucket"],
                        interval=intervals
                        ), key=lambda item: int(os.path.basename(item).split(".")[-3].split("_")[1])
        )
    output:
       #report = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv"
        report = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.recal_data.csv"),
    params:
        java_opt = "-Xms3000m",
        gatk     = config["gatk"]
    resources:
         time   = 10,
         mem_mb = 4000
    run:
        # separate reports by -I
        reports = " -I ".join(map(str,input.bqsr_reports))

        # primary logs dir
       #log_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "logs"
       #)

        shell(f'''
            {{params.gatk}} --java-options {{params.java_opt}} \
                GatherBQSRReports \
                -I {reports} \
                -O {{output.report}}
        ''')
        
rule apply_bqsr:
    input:
        sorted_bam = "{bucket}/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
       #report     = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv",
        report     = S3.remote(f"{{bucket}}/wgs/{breed}/{{sample_name}}/{{ref}}/logs/{{sample_name}}.{{ref}}.recal_data.csv"),
        interval   = ancient("results/seq_group/with_unmap/{interval}.tsv")
    output:
        recal_bam = "{bucket}/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
        recal_bai = "{bucket}/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bai"
    params:
        java_opt  = "-Xms3000m",
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"]
    resources:
         time   = 180,
         mem_mb = 10000
    run:
        # separate content of interval by -L
        with open(input.interval,"r") as f:
            ival = f.read().strip().replace("\t", " -L ")

        shell(f'''
            {{params.gatk}} --java-options {{params.java_opt}} \
                ApplyBQSR \
                -R {{params.ref_fasta}} \
                -I {{input.sorted_bam}} \
                -O {{output.recal_bam}} \
                -L {ival} \
                -bqsr {{input.report}} \
                --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                --add-output-sam-program-record \
                --create-output-bam-md5 \
                --use-original-qualities
        ''')

rule gather_bam_files:
    input:
        recal_bams = sorted(
                expand("{bucket}/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
                    sample_name=units["sample_name"].values[0],
                    ref=config["ref"], bucket=config["bucket"],
                    interval=intervals
                    ), key=lambda item: int(os.path.basename(item).split(".")[2].split("_")[1])
        )
    output:
       #final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam",
       #final_bai = "results/gather_bam_files/{sample_name}.{ref}.bai",
       #final_md5 = "results/gather_bam_files/{sample_name}.{ref}.bam.md5"
        final_bam = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"),
        final_bai = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"),
        final_md5 = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam.md5"),
    params:
        comp_level = 5,
        java_opt   = "-Xms2000m",
        gatk       = config["gatk"]
    threads: 4
    resources:
         time   = 180,
         mem_mb = 16000
    run:
        # separate bams by -I
        bams = " -I ".join(map(str,input.recal_bams))
        
        # primary bam dir
       #bam_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "bam"
       #)

        shell(f'''
            {{params.gatk}} --java-options "-Dsamjdk.compression_level={{params.comp_level}} {{params.java_opt}}" \
                GatherBamFiles \
                --INPUT {bams} \
                --OUTPUT {{output.final_bam}} \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true
        ''')

rule coverage_depth_and_flagstat:
    input:
       #final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam"
        final_bam = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"),
        final_bai = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"),
    output:
       #doc_smry = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_summary",
       #doc_stat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_statistics",
       #flagstat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat"
        doc_smry = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.depthofcoverage.sample_summary"),
        doc_stat = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.depthofcoverage.sample_statistics"),
        flagstat = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.flagstat"),
    params:
        doc_base   = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.depthofcoverage"),
        comp_level = 5,
        java_opt   = "-Xmx32000m",
        gatk3      = config["gatk3"],
        ref_fasta  = config["ref_fasta"],
    threads: 8
    resources:
         time   = 360,
         mem_mb = 32000
    run:
        # primary bam dir
       #bam_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "bam"
       #)
        
        shell(f'''
            java -jar {{params.java_opt}} {{params.gatk3}} \
                -T DepthOfCoverage \
                -R {{params.ref_fasta}} \
                -omitBaseOutput \
                -omitLocusTable \
                -omitIntervals \
                -I {{input.final_bam}} \
                -o {{params.doc_base}} \
                -ct 5 \
                -ct 15 \
                -ct 30 \
                -nt 8
        
            java -jar {{params.java_opt}} {{params.gatk3}} \
                -T FlagStat \
                -R {{params.ref_fasta}} \
                -I {{input.final_bam}} \
                -o {{output.flagstat}} \
                -nct 8
        ''')

rule haplotype_caller:
    input:
       #final_bam = "{bucket}/gather_bam_files/{sample_name}.{ref}.bam",
        final_bam = S3.remote(f"{{bucket}}/wgs/{breed}/{{sample_name}}/{{ref}}/bam/{{sample_name}}.{{ref}}.bam"),
        final_bai = S3.remote(f"{{bucket}}/wgs/{breed}/{{sample_name}}/{{ref}}/bam/{{sample_name}}.{{ref}}.bai"),
        interval  = f"{os.path.join(config['hc_intervals'],'{hc_interval}')}.interval_list"
    output:
        hc_gvcf = "{bucket}/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz"
    params:
        java_opt  = "-Xmx10G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
    threads: 4
    resources:
         time   = 720,
         mem_mb = 8000
    shell:
        '''
            {params.gatk} --java-options "{params.java_opt}" \
                HaplotypeCaller \
                -R {params.ref_fasta} \
                -I {input.final_bam} \
                -L {input.interval} \
                -O {output.hc_gvcf} \
                -contamination 0 -ERC GVCF
        '''

rule merge_gvcfs:
    input:
        flagstat = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.flagstat"),
       #flagstat  = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat",
        hc_gvcfs = sorted(
            expand(
                "{bucket}/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz",
                sample_name=units["sample_name"].values[0],
                ref=config["ref"], bucket=config["bucket"],
                hc_interval=hc_intervals
                ), key=lambda item: int(item.split("/")[-2].split("-")[0])
        )
    output:
       #final_gvcf     = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz",
       #final_gvcf_tbi = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz.tbi"
        final_gvcf     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz"),
        final_gvcf_tbi = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz.tbi"),
    params:
        java_opt  = "-Xmx2000m",
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
    threads: 4
    resources:
         time   = 120,
         mem_mb = 4000
    run:
        # separate bams by -I
        gvcfs = " -I ".join(map(str,input.hc_gvcfs))
        
        # primary gvcf dir
       #gvcf_outdir = os.path.join(
       #        config["primary"],
       #        config["ref"],
       #        BREED,
       #        wildcards.sample_name,
       #        "gvcf"
       #)

        shell(f'''
            {params.gatk} --java-options {params.java_opt}  \
                MergeVcfs \
                --INPUT {gvcfs} \
                --OUTPUT {{output.final_gvcf}}
        ''')

rule genotype_gvcfs:
    input:
       #ival_db  = "results/import_gvcfs/{interval}",
       #interval = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/CanFam4/Intervals/{interval}.interval_list"
        final_gvcf = S3.remote(
            expand(
                "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.{ext}",
                sample_name=units["sample_name"].values[0],
                ref=config["ref"], breed=breed,
                ext=["gz","gz.tbi"], bucket=config["bucket"]
            )
        ),
       #interval   = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/MoneyCF4/Intervals/{interval}.interval_list"
        interval   = f"{os.path.join(config['hc_intervals'],'{hc_interval}')}.interval_list"
    output:
        vcf = "{bucket}/genotype_gvcfs/{ref}/{hc_interval}/output.vcf.gz",
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"]
    threads: 6
    resources:
         time   = 60,
         mem_mb = 60000
    run:
        # need the index for tool but not included in command
        gvcf = [i for i in input.final_gvcf if ".tbi" not in i][0]

        shell(f'''
            {{params.gatk}} --java-options "-Xmx50g -Xms50g" \
            GenotypeGVCFs \
                -R {{params.ref_fasta}} \
                -O {{output.vcf}} \
                -G StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V {gvcf} \
                -L {{input.interval}}
        ''')

rule fltr_make_sites_only:
    input:
        vcf = "{bucket}/genotype_gvcfs/{ref}/{hc_interval}/output.vcf.gz"
    output:
        var_filtrd_vcf = "{bucket}/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.variant_filtered.vcf.gz",
        sites_only_vcf = "{bucket}/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.sites_only.variant_filtered.vcf.gz"
    params:
        gatk       = config["gatk"],
        excess_het = config["excess_het_threshold"],
        ref_fasta  = config["ref_fasta"]
    resources:
         time   = 30,
         mem_mb = 6000
    shell:
        '''
            {params.gatk} --java-options "-Xmx3g -Xms3g" \
            VariantFiltration \
                --filter-expression "ExcessHet > {params.excess_het}" \
                --filter-name ExcessHet \
                -O {output.var_filtrd_vcf} \
                -V {input.vcf}

            {params.gatk} --java-options "-Xmx3g -Xms3g" \
            MakeSitesOnlyVcf \
                --INPUT {output.var_filtrd_vcf} \
                --OUTPUT {output.sites_only_vcf}
        '''

rule sites_only_gather_vcf:
    input:
        sites_only_vcf = expand(
            "{bucket}/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.sites_only.variant_filtered.vcf.gz",
            ref=config["ref"], bucket=config["bucket"],
            hc_interval=hc_intervals
        )
    output:
        unsrtd_sites_only_vcf = "{bucket}/sites_only_gather_vcf/{ref}/unsrtd.sites_only.vcf.gz",
        gather_sites_only_vcf = "{bucket}/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz.tbi"
    params:
        gatk   = config["gatk"],
        picard = config["picard"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    run:
        vcfs = " --input ".join(map(str,input.sites_only_vcf))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {vcfs} \
                --output {{output.unsrtd_sites_only_vcf}}

            java -jar {{params.picard}} \
                SortVcf \
                I={{output.unsrtd_sites_only_vcf}} \
                O={{output.gather_sites_only_vcf}}

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {{output.gather_sites_only_vcf}}

            rm -f results/sites_only_gather_vcf/tmp.vcf.gz
        ''')

rule indels_var_recal:
    input:
        gather_sites_only_vcf = "{bucket}/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
    output:
        indels_recal    = "{bucket}/recal/{ref}/indels/indels.recal",
        indels_tranches = "{bucket}/recal/{ref}/indels/indels.tranches"
    params:
        gatk                    = config["gatk"],
        recal_tranche_values    = config["indel_recalibration_tranche_values"],
        recal_annotation_values = config["indel_recalibration_annotation_values"],
        dbsnp_indels_vcf        = config["dbsnp_indels_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    run:
        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
        an_values = " -an ".join(map(str,params.recal_annotation_values))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx24g -Xms24g" \
                VariantRecalibrator \
                -V {{input.gather_sites_only_vcf}} \
                -O {{output.indels_recal}} \
                --tranches-file {{output.indels_tranches}} \
                --trust-all-polymorphic \
                -tranche {tranche_values} \
                -an {an_values} \
                -mode INDEL \
                --max-gaussians 4 \
                --resource:dbsnp,known=false,training=true,truth=true,prior=10 {params.dbsnp_indels_vcf}
        ''')


rule snps_var_recal:
    input:
        gather_sites_only_vcf = "{bucket}/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
    output:
        snps_recal    = "{bucket}/recal/{ref}/snps/snps.recal",
        snps_tranches = "{bucket}/recal/{ref}/snps/snps.tranches"
    params:
        gatk                    = config["gatk"],
        recal_tranche_values    = config["snp_recalibration_tranche_values"],
        recal_annotation_values = config["snp_recalibration_annotation_values"],
        dbsnp_snp_vcf           = config["dbsnp_snp_vcf"],
        broad_snp_vcf           = config["broad_snp_vcf"],
        axelsson_snp_vcf        = config["axelsson_snp_vcf"],
        illumina_snp_vcf        = config["illumina_snp_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 16000
    run:
        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
        an_values = " -an ".join(map(str,params.recal_annotation_values))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx12g -Xms3g" \
                VariantRecalibrator \
                -V {{input.gather_sites_only_vcf}} \
                -O {{output.snps_recal}} \
                --tranches-file {{output.snps_tranches}} \
                --trust-all-polymorphic \
                -tranche {tranche_values} \
                -an {an_values} \
                -mode SNP \
                --max-gaussians 6 \
                --resource:illumina,known=true,training=true,truth=true,prior=15.0 {params.illumina_snp_vcf} \
                --resource:broad,known=true,training=true,truth=true,prior=10.0 {params.broad_snp_vcf} \
                --resource:axelsson,known=false,training=true,truth=true,prior=8.0 {params.axelsson_snp_vcf} \
                --resource:dbSNP146,known=true,training=true,truth=true,prior=12.0 {params.dbsnp_snp_vcf}
        ''')

rule apply_recal:
    input:
        input_vcf       = "{bucket}/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.variant_filtered.vcf.gz",
        indels_recal    = "{bucket}/recal/{ref}/indels/indels.recal",
        indels_tranches = "{bucket}/recal/{ref}/indels/indels.tranches",
        snps_recal      = "{bucket}/recal/{ref}/snps/snps.recal",
        snps_tranches   = "{bucket}/recal/{ref}/snps/snps.tranches"
    output:
        recal_vcf       = "{bucket}/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz",
        recal_vcf_index = "{bucket}/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz.tbi"
    params:
        gatk               = config["gatk"],
        indel_filter_level = config["indel_filter_level"],
        snp_filter_level   = config["snp_filter_level"]
    resources:
         time   = 30,
         mem_mb = 16000
    shell:
        '''
            mkdir -p results/apply_recal/{wildcards.hc_interval}/

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -O results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf \
                -V {input.input_vcf} \
                --recal-file {input.indels_recal} \
                --tranches-file {input.indels_tranches} \
                --truth-sensitivity-filter-level {params.indel_filter_level} \
                --create-output-variant-index true \
                -mode INDEL

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -O {output.recal_vcf} \
                -V results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf \
                --recal-file {input.snps_recal} \
                --tranches-file {input.snps_tranches} \
                --truth-sensitivity-filter-level {params.snp_filter_level} \
                --create-output-variant-index true \
                -mode SNP

            rm -f results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf
        '''

rule final_gather_vcfs:
    input:
        recal_vcfs = sorted(
            expand(
                "{bucket}/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz", 
                ref=config["ref"], bucket=config["bucket"],
                hc_interval=hc_intervals
            )
        )
    output:
        final_vcf       = "{bucket}/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
        final_vcf_index = "{bucket}/final_gather_vcfs/joint_genotype.{ref}.vcf.gz.tbi"
    params:
        gatk = config["gatk"],
    threads: 4
    resources:
         time   = 240,
         mem_mb = 22000
    run:
        vcfs = " --input ".join(map(str,input.recal_vcfs))

        shell(f'''
            set -e
            set -o pipefail

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {vcfs} \
                --output {{output.final_vcf}}

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {{output.final_vcf}}
        ''')

rule collect_metrics_on_vcf:
    input:
        final_vcf       = "{bucket}/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
        final_vcf_index = "{bucket}/final_gather_vcfs/joint_genotype.{ref}.vcf.gz.tbi"
    output:
        detail_metrics  = "{bucket}/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_detail_metrics",
        summary_metrics = "{bucket}/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_summary_metrics"
    params:
        gatk               = config["gatk"],
        dbsnp_snp_vcf      = config["dbsnp_snp_vcf"],
        ref_dict           = config["ref_dict"],
        eval_interval_list = config["eval_interval_list"],
        metrics_prefix     = f"{config['bucket']}/collect_metrics_on_vcf/joint_genotype.{config['ref']}"
    threads: 8
    resources:
         time   = 360,
         mem_mb = 22000
    shell:
        '''
            mkdir -p results/collect_metrics_on_vcf/

            {params.gatk} --java-options "-Xmx18g -Xms6g" \
                CollectVariantCallingMetrics \
                --INPUT {input.final_vcf} \
                --DBSNP {params.dbsnp_snp_vcf} \
                --SEQUENCE_DICTIONARY {params.ref_dict} \
                --OUTPUT {params.metrics_prefix} \
                --THREAD_COUNT 8 \
                --TARGET_INTERVALS {params.eval_interval_list}
        '''

rule vep_final_vcf:
    input:
        final_vcf = "{bucket}/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
    output:
        vep_vcf     = "{bucket}/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
        vep_vcf_tbi = "{bucket}/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz.tbi"
    params:
        conda_vep = config["conda_vep"],
        ref_fasta = config["ref_fasta"],
        ref_gtf   = config["ref_gtf"]
    threads: 12
    resources:
         time   = 4320,
         mem_mb = 60000
    run:
        import os
        out_name = os.path.splitext(output.vep_vcf)[0] 
        
        # if canfam4 VEP with the gtf
        if "canfam4" in wildcards.ref:
            shell(f'''
                set +eu

                eval "$(conda shell.bash hook)"
                conda activate {params.conda_vep}

                vep \
                    -i {{input.final_vcf}} \
                    -o {out_name} \
                    --gtf {{params.ref_gtf}} \
                    --fasta {{params.ref_fasta}} \
                    --everything \
                    --force_overwrite \
                    --vcf \
                    --dont_skip

                bgzip --threads 12 -c {out_name} > {{output.vep_vcf}} &&
                tabix -p vcf {{output.vep_vcf}}
            ''')
        # else VEP with the canfam3 annotation from ensembl
        else:
            shell(f'''
                set +eu

                eval "$(conda shell.bash hook)"
                conda activate {params.conda_vep}

                vep \
                  -i {{input.final_vcf}} \
                  --cache \
                  --everything \
                  -o {out_name} \
                  --vcf \
                  --species=canis_familiaris \
                  --offline \
                  --dont_skip

                bgzip --threads 12 -c {out_name} > {{output.vep_vcf}} &&
                tabix -p vcf {{output.vep_vcf}}
            ''')

rule select_common_vars:
    input:
        pop_vcf     = config["pop_vcf"],
       #pop_vcf_tbi = config["pop_vcf_tbi"]
    output:
        common_vars     = "{bucket}/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz",
        common_vars_tbi = "{bucket}/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz.tbi"
    params:
        conda_vep = config["conda_vep"],
        ref_fasta = config["ref_fasta"],
    threads: 8
    resources:
         time   = 2160,
         mem_mb = 32000
    run:
        # check if common_vars already exists for population vcf - generate
        # if not and copy if does
        l = [config["common_vars"],config["common_vars_tbi"]]
        
        if all([os.path.isfile(f) for f in l]):
            shell(f'''
                cp {config["common_vars"]} {{output.common_vars}}
                cp {config["common_vars_tbi"]} {{output.common_vars_tbi}}
            ''')
        else:
            shell('''
                set +eu

                eval "$(conda shell.bash hook)"
                conda activate {params.conda_vep}
                
                module load bcftools

                bcftools view -Oz \
                    -i 'AF[*]>0.005' \
                    {input.pop_vcf} \
                    -o {output.common_vars}

                tabix -p vcf {output.common_vars}
            ''')


rule select_variants_to_table:
    input:            
       #final_vcf      = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
       #common_vars    = f"results/vep_final_vcf/joint_genotype.{config['ref']}.vep.af_nonmajor.vcf.gz",
       #detail_metrics = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
        final_vcf      = "{bucket}/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
        common_vars    = "{bucket}/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz",
        detail_metrics = "{bucket}/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_detail_metrics",
    output:
       #unique_vars           = "results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.unique_vars.vcf",
       #rare_and_common_vars  = "results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.rare_and_common_vars.vcf",
       #rare_vars             = "results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.rare_vars.vcf",
       #unique_vars_vep_split = "results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.unique_vars.vep_split.txt",
       #rare_vars_vep_split   = "results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.rare_vars.vep_split.txt"
        unique_vars           = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf",
        rare_and_common_vars  = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_and_common_vars.vcf",
        rare_vars             = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf",
        unique_vars_vep_split = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vep_split.txt",
        rare_vars_vep_split   = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vep_split.txt"
    params:
        gatk = config["gatk"],
        tmp_dir = "/dev/shm/money_vars_{ref}",
        ref_fasta = config["ref_fasta"],
        pop_vcf = config["pop_vcf"]
    threads: 4
    resources:
         time   = 2880,
         mem_mb = 16000
    shell:
        '''
            module load bcftools    
            mkdir -p {params.tmp_dir} results/select_vars_to_table

            # all unique variants
            {params.gatk} SelectVariants \
                -R {params.ref_fasta} \
                -V {input.final_vcf} \
                -disc {params.pop_vcf} \
                --tmp-dir {params.tmp_dir} \
                -O {output.unique_vars}

            # all common and  rare variants
            {params.gatk} SelectVariants \
                -R {params.ref_fasta} \
                -V {input.final_vcf} \
                -disc {input.common_vars} \
                -O {output.rare_and_common_vars}

            # rare only variants
            {params.gatk} SelectVariants \
                -R {params.ref_fasta} \
                -V {output.rare_and_common_vars} \
                -disc {output.unique_vars} \
                -O {output.rare_vars}

            # vars to table - all unique
            bcftools +split-vep \
                {output.unique_vars} \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
                -d -A tab \
                -o {output.unique_vars_vep_split}

            # vars to table - rare only
            bcftools +split-vep \
                {output.rare_vars} \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
                -d -A tab \
                -o {output.rare_vars_vep_split}
        '''

rule bgzip_and_tabix:
    input:            
        unique_vars = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf",
        rare_vars   = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf",
    output:
        unique_vars_gz  = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz",
        unique_vars_tbi = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz.tbi",
        rare_vars_gz    = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz",
        rare_vars_tbi   = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz.tbi",
    params:
        conda_vep = config["conda_vep"],
    threads: 8
    resources:
         time   = 60,
         mem_mb = 16000
    shell:
        '''
            set +eu

            eval "$(conda shell.bash hook)"
            conda activate {params.conda_vep}
            
            module load bcftools

            # bgzip and tabix the unique vcf
            bgzip -c {input.unique_vars} > {output.unique_vars_gz}
            tabix -p vcf {output.unique_vars_gz}
        
            # and the rare vcf
            bgzip -c {input.rare_vars} > {output.rare_vars_gz}
            tabix -p vcf {output.rare_vars_gz}
        '''

rule final_output:
    input:
       #unique_vars_vep_split = f"results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.unique_vars.vep_split.txt",
       #rare_vars_vep_split   = f"results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.rare_vars.vep_split.txt"
        unique_vars_vep_split = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vep_split.txt",
        rare_vars_vep_split   = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vep_split.txt"
    output:
       #excel_sheet = f"results/final_output/{units['sample_name'].values[0]}_{config['ref']}_tables.xlsx",
        excel_sheet = "{bucket}/final_output/{ref}/{sample_out}_{ref}_tables.xlsx",
       #excel_sheet = "{sample_out}_{ref}_tables.xlsx",
    params:
        conda_vep = config["conda_vep"],
        base_name = "{bucket}/final_output/{ref}"
    threads: 4
    resources:
         time   = 30,
         mem_mb = 16000
    run:
       #import os
       #import pandas as pd
        
        # first add header to the rare and unique tables
        if "canfam3" in wildcards.ref:
            header = [
              'chrom', 'pos', 'ref', 'alt', 'ac', 
              'allele', 'consequence', 'impact', 'symbol', 
              'gene', 'feature_type', 'feature', 'biotype', 
              'exon', 'intron', 'hgvsc', 'hgvsp', 
              'cdna_position', 'cds_position', 'protein_position', 'amino_acids', 
              'codons', 'existing_variation', 'distance', 'strand', 
              'flags', 'variant_class', 'symbol_source', 'hgnc_id', 
              'canonical','mane', 'tsl', 'appris', 
              'ccds', 'ensp', 'swissprot', 'trembl',
              'uniparc', 'gene_pheno', 'sift', 'domains', 
              'mirna', 'af', 'afr_af', 'amr_af', 
              'eas_af', 'eur_af', 'sas_af', 'aa_af', 
              'ea_af', 'gnomad_af', 'gnomad_afr_af', 'gnomad_amr_af', 
              'gnomad_asj_af', 'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_nfe_af', 
              'gnomad_oth_af', 'gnomad_sas_af', 'max_af', 'max_af_pops', 
              'clin_sig', 'somatic', 'pheno', 'pubmed', 
              'check_ref', 'motif_name', 'motif_pos', 'high_inf_pos',
              'motif_score_change'
            ]
        else:
            header = [
                'chrom', 'pos', 'ref', 'alt', 'ac',
                'allele', 'consequence', 'impact', 'symbol', 'gene', 
                'feature_type', 'feature', 'biotype', 'exon', 'intron', 
                'hgvsc', 'hgvsp', 'cdna_position', 'cds_position', 
                'protein_position', 'amino_acids', 'codons', 'existing_variation',
                'distance', 'strand', 'flags', 'variant_class', 'symbol_source', 
                'hgnc_id', 'canonical', 'mane_select', 'mane_plus_clinical', 
                'tsl', 'appris', 'ccds', 'ensp', 'swissprot', 'trembl', 'uniparc', 
                'uniprot_isoform', 'source', 'gene_pheno', 'sift', 'polyphen', 
                'domains', 'mirna', 'hgvs_offset', 'af', 'afr_af', 'amr_af', 
                'eas_af', 'eur_af', 'sas_af', 'aa_af', 'ea_af', 'gnomad_af', 
                'gnomad_afr_af', 'gnomad_amr_af', 'gnomad_asj_af', 
                'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_nfe_af', 
                'gnomad_oth_af', 'gnomad_sas_af', 'max_af', 'max_af_pops', 
                'clin_sig', 'somatic', 'pheno', 'pubmed', 'check_ref', 
                'motif_name', 'motif_pos', 'high_inf_pos', 
                'motif_score_change', 'transcription_factors', 
                'canFam4_genes_ncbi.gtf.gz'
            ]
     
        # for each table, add header and print the rest of table   
        dfs = {}
        for f in [input.unique_vars_vep_split, input.rare_vars_vep_split]:
            tmp = os.path.basename(f).replace(".txt",".reform.txt")
            reform = os.path.join(params.base_name,tmp)
            # get name of each table as key for dict
            name = tmp.split(".vep_split")[0]
            # open table, add header, and write to reform 
            with open(f, "r") as infile, open(reform, "w") as outfile:
                print("\t".join(header), file=outfile)
                for line in infile:
                    print(line.strip(), file=outfile)
            # read fixed table to dataframe
            dfs[name] = pd.read_csv(reform, sep="\t")

        # write each df to the same excel sheet
        with pd.ExcelWriter(
                output.excel_sheet,
                engine="xlsxwriter",
                options={"strings_to_formulas": False}
            ) as writer:
            for k,v in dfs.items():
                v.to_excel(writer,
                           sheet_name=k,
                           index=False)
            writer.save()

rule manifest_and_upload:
    input:
        vep_vcf         = "{bucket}/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
        vep_vcf_tbi     = "{bucket}/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz.tbi",
        unique_vars_gz  = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz",
        unique_vars_tbi = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz.tbi",
        rare_vars_gz    = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz",
        rare_vars_tbi   = "{bucket}/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz.tbi",
        excel_sheet     = "{bucket}/final_output/{ref}/{sample_out}_{ref}_tables.xlsx",
    output:
       #manifest      = "{bucket}/{sample_out}_{ref}_manifest.txt",
       #money_tar_gz  = "{bucket}/final_output/{ref}/{sample_out}_{ref}_K9MM.tar.gz",
        manifest     = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_manifest.txt"),
        money_tar_gz = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_K9MM.tar.gz"),
    params:
       #base_dir = "results/final_output/{ref}/upload/"
    run:
       #os.makedirs(os.path.dirname(output.manifest),exist_ok=True)
        
        with open(output.manifest, "w") as outfile:
            s = f'''
                Included files:
                1) results/vep_final_vcf/joint_genotype.{wildcards.ref}.vep.vcf.gz - annotated VCF and index (.tbi) for all {wildcards.sample_out} variants
                2) results/select_vars_to_table/{wildcards.sample_out}.{wildcards.ref}.unique_vars.vcf.gz - annotated VCF and index (.tbi) for all unique {wildcards.sample_out} variants
                3) results/select_vars_to_table/{wildcards.sample_out}.{wildcards.ref}.rare_vars.vcf.gz - annotated VCF and index (.tbi) for all rare {wildcards.sample_out} variants
                4) results/final_output/{wildcards.sample_out}_{wildcards.ref}_tables.xlsx - contains two sheets
                \t i) {wildcards.sample_out}_{wildcards.ref}.unique_vars - table version of item 2
                \t ii) {wildcards.sample_out}_{wildcards.ref}.rare_vars - table version of item 3
            '''
            print(s, file=outfile)

        shell('''
            tar -czvf {output.money_tar_gz} \
                {input.vep_vcf} \
                {input.vep_vcf_tbi} \
                {input.unique_vars_gz} \
                {input.unique_vars_tbi} \
                {input.rare_vars_gz} \
                {input.rare_vars_tbi} \
                {input.excel_sheet} \
                {output.manifest}
        ''')

rule save_jobs:
    input:
        manifest     = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_manifest.txt"),
        money_tar_gz = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_K9MM.tar.gz"),
    output:
        touch("{bucket}/save_jobs/{breed}_{sample_out}_{ref}_jobs.done")
    params:
        fastqc = config["fastqc"],
        alias  = config["alias"]
   #resources:
   #     time   = 360,
   #     mem_mb = 6000, 
   #     cpus   = 1
    shell:
        '''
            mcli cp --recursive src/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_out}/{wildcards.ref}/jobs/src
            
            mcli cp --recursive slurm.go_money/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_out}/{wildcards.ref}/jobs/slurm.go_money
            
            mcli cp --recursive Jobs/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_out}/{wildcards.ref}/jobs/slurm_logs
            
            mcli cp {wildcards.ref}_money.yaml {wildcards.breed}_{wildcards.sample_out}.go_make_money.slurm go_make_money.smk input.tsv \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_out}/{wildcards.ref}/jobs/
        '''
