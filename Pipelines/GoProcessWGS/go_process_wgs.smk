import pandas as pd
import os

localrules: save_jobs

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
       #S3.remote(
       #    expand(
       #        "{bucket}/wgs/{breed}/{u.sample_name}/{ref}/gvcf/{u.sample_name}.{ref}.g.vcf.{ext}",
       #        u=units.itertuples(), ref=config["ref"],
       #        bucket=config["bucket"], breed=breed,
       #        ext=["gz","gz.tbi"]
       #    )
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
       #expand(
       #    "{bucket}/save_jobs/{breed}_{sample_out}_{ref}_jobs.done",
       #    sample_out=units['sample_name'].values[0], 
       #    bucket=config['bucket'], ref=config["ref"], 
       #    breed=breed, 
       #)
        expand(
            "{bucket}/save_jobs/{breed}/{u.sample_name}_{ref}_jobs.done",
            u=units.itertuples(), bucket=config['bucket'], 
            ref=config["ref"], breed=breed, 
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

rule save_jobs:
    input:
        final_gvcf     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz"),
        final_gvcf_tbi = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz.tbi"),
       #manifest     = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_manifest.txt"),
       #money_tar_gz = S3.remote("{bucket}/wgs/{breed}/{sample_out}/{ref}/money/{sample_out}_{ref}_K9MM.tar.gz"),
    output:
        touch("{bucket}/save_jobs/{breed}/{sample_name}_{ref}_jobs.done")
    params:
        alias  = config["alias"]
    shell:
        '''
            mcli cp --recursive ./src/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/jobs/src/
            
            mcli cp --recursive ./slurm.go_wgs/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/jobs/slurm.go_money/
            
            mcli cp --recursive ./Jobs/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/jobs/slurm_logs/
            
            mcli cp {wildcards.ref}_config.yaml {wildcards.breed}_{wildcards.sample_name}.go_process_wgs.slurm go_process_wgs.smk input.tsv \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/jobs/
        '''
