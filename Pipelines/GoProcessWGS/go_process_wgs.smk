import pandas as pd
import os

include: "src/utils.py"

units = pd.read_table(config["units"],dtype=str).set_index("readgroup_name",drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

# get breed from units
BREED = units["breed"].unique()[0]

sequence_grouping(config["ref_dict"])
# get sequence group intervals without unmapped, with unmapped, and hc caller intervals
intervals, = glob_wildcards(os.path.join("results/seq_group/no_unmap","{interval}.tsv"))
unmap_intervals, = glob_wildcards(os.path.join("results/seq_group/with_unmap","{interval}.tsv"))
hc_intervals, = glob_wildcards(os.path.join(config["hc_intervals"],"{hc_interval}.interval_list"))

rule all:
    input:
       ## fastqc 
        expand("results/fastqc/{u.sample_name}/{u.readgroup_name}/qc.done",
            u=units.itertuples()
        ),
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
        expand("results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
            sample_name=units["sample_name"].values[0],ref=config["ref"],
            interval=unmap_intervals
        ),
       ## gather bams
       #expand("results/gather_bam_files/{u.sample_name}.{ref}.bam",
       #    u=units.itertuples(),ref=config["ref"]
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
        expand("results/merge_gvcfs/{u.sample_name}.{ref}.g.vcf.gz",
            u=units.itertuples(),ref=config["ref"]
        )


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        "results/fastqc/{sample_name}/{readgroup_name}/qc.done"
    params:
        fastqc  = config["fastqc"],
    resources:
         time   = 120,
         mem_mb = 6000, 
         cpus   = 4
    run:
        # get flowcell
        flowcell = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"flowcell"].values[0]

        # output path to fastqc by flowcell
        qc_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "fastqc",
                flowcell
        )
        
        # output path to fastq by flowcell
        fq_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "fastq",
                flowcell
        )

        shell(f'''
            set -e

            mkdir -p {qc_outdir} {fq_outdir}

            {{params.fastqc}} --outdir {qc_outdir} \
                {{input.r1}} {{input.r2}}

            cp -t {fq_outdir} \
                {{input.r1}} {{input.r2}}

            touch {{output}}
        ''')

rule fastqs_to_ubam:
    input:
        unpack(get_fastq),
       #"results/fastqc/{sample_name}/{readgroup_name}/qc.done"
    output:
        ubam = "results/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    params:
        java_opt  = "-Xms6000m",
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
        tmp_dir   = f"/scratch.global/friedlab_{os.environ['USER']}/{{readgroup_name}}.ubam.tmp"
    threads: 6
    resources:
         time   = 540,
         mem_mb = 24000
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

            {{params.gatk}} --java-options {params.java_opt} \
                FastqToSam \
                --TMP_DIR {params.tmp_dir} \
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
        ubam = "results/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    output:
        bwa_log = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.aligned.unmerged.bwa.stderr.log",
        bam     = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.aligned.unmerged.bam"
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
        ubam = "results/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam",
        bam  = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{readgroup_name}.aligned.unmerged.bam"
    output:
        merged_bam = "results/merge_bam_alignment/{sample_name}/{readgroup_name}.{ref}.merged.unsorted.bam"
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
                        "results/merge_bam_alignment/{sample_name}/{readgroup_name}.{ref}.merged.unsorted.bam",
                        sample_name=units["sample_name"].values[0],
                        readgroup_name=list(units["readgroup_name"]),
                        ref=config["ref"]
                      )
    output:
        dedup_bam = "results/mark_duplicates/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
        metrics   = "results/mark_duplicates/{sample_name}.{ref}.duplicate_metrics"
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m -Xmx16g",
        gatk       = config["gatk"],
        tmp_dir    = "/dev/shm/{sample_name}.md.tmp"
    threads: 4
    resources:
         time   = 720,
         mem_mb = 60000
    run:
        # separate bams by --INPUT
        bams = " --INPUT ".join(map(str,input.merged_bams))

        # primary logs dir
        log_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "logs"
        )

        shell(f'''
            set -e

            mkdir -p {{params.tmp_dir}} {log_outdir}

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

            cp {{output.metrics}} {log_outdir}
        ''')

rule sort_and_fix_tags:
    input:
        dedup_bam = "results/mark_duplicates/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam"
    output:
        sorted_bam = "results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        sorted_bai = "results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai"
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m",
        gatk       = config["gatk"],
        ref_fasta  = config["ref_fasta"],
        tmp_dir    = f"/scratch.global/friedlab_{os.environ['USER']}/{{sample_name}}.sort.tmp"
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
        sorted_bam = "results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        interval   = ancient("results/seq_group/no_unmap/{interval}.tsv")
    output:
        recal_csv = "results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv"
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
         mem_mb = 10000
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
                    expand("results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
                        sample_name=units["sample_name"].values[0],
                        ref=config["ref"],
                        interval=intervals
                        ), key=lambda item: int(os.path.basename(item).split(".")[-3].split("_")[1])
        )
    output:
        report = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv"
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
        log_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "logs"
        )

        shell(f'''
            set -e

            {{params.gatk}} --java-options {{params.java_opt}} \
                GatherBQSRReports \
                -I {reports} \
                -O {{output.report}}

            cp {{output.report}} {log_outdir}
        ''')
        
rule apply_bqsr:
    input:
        sorted_bam = "results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        report     = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv",
        interval   = ancient("results/seq_group/with_unmap/{interval}.tsv")
    output:
        recal_bam = "results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
        recal_bai = "results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bai"
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
                expand("results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
                    sample_name=units["sample_name"].values[0],
                    ref=config["ref"],
                    interval=intervals
                    ), key=lambda item: int(os.path.basename(item).split(".")[2].split("_")[1])
        )
    output:
        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam",
        final_bai = "results/gather_bam_files/{sample_name}.{ref}.bai",
        final_md5 = "results/gather_bam_files/{sample_name}.{ref}.bam.md5"
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
        bam_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "bam"
        )

        shell(f'''
            set -e

            mkdir -p {bam_outdir}

            {{params.gatk}} --java-options "-Dsamjdk.compression_level={{params.comp_level}} {{params.java_opt}}" \
                GatherBamFiles \
                --INPUT {bams} \
                --OUTPUT {{output.final_bam}} \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true
            
            cp -t {bam_outdir} {{output}}

        ''')

rule coverage_depth_and_flagstat:
    input:
        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam"
    output:
        doc_smry = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_summary",
        doc_stat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_statistics",
        flagstat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat"
    params:
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
        bam_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "bam"
        )
        
        shell(f'''
            set -e

            java -jar {{params.java_opt}} {{params.gatk3}} \
                -T DepthOfCoverage \
                -R {{params.ref_fasta}} \
                -omitBaseOutput \
                -omitLocusTable \
                -omitIntervals \
                -I {{input.final_bam}} \
                -o results/coverage_depth_and_flagstat/{{wildcards.sample_name}}.{{wildcards.ref}}.depthofcoverage \
                -ct 5 \
                -ct 15 \
                -ct 30 \
                -nt 8
        
            java -jar {{params.java_opt}} {{params.gatk3}} \
                -T FlagStat \
                -R {{params.ref_fasta}} \
                -I {{input.final_bam}} \
                -o results/coverage_depth_and_flagstat/{{wildcards.sample_name}}.{{wildcards.ref}}.flagstat \
                -nct 8

            cp -t {bam_outdir} {{output}}
        ''')

rule haplotype_caller:
    input:
        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam",
        interval = f"{os.path.join(config['hc_intervals'],'{hc_interval}')}.interval_list"
    output:
        hc_gvcf = "results/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz"
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
        flagstat  = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat",
        hc_gvcfs = sorted(
                expand("results/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz",
                    sample_name=units["sample_name"].values[0],
                    ref=config["ref"],
                    hc_interval=hc_intervals
                    ), key=lambda item: int(item.split("/")[-2].split("-")[0])
        )
    output:
        final_gvcf     = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz",
        final_gvcf_tbi = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz.tbi"
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
        gvcf_outdir = os.path.join(
                config["primary"],
                config["ref"],
                BREED,
                wildcards.sample_name,
                "gvcf"
        )

        shell(f'''
            set -e

            mkdir -p {gvcf_outdir}

            {params.gatk} --java-options {params.java_opt}  \
                MergeVcfs \
                --INPUT {gvcfs} \
                --OUTPUT {{output.final_gvcf}}

            cp -t {gvcf_outdir} {{output}}
        ''')

