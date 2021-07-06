import pandas as pd
import os

localrules: manifest_and_upload

include: "src/utils.py"

units = pd.read_table(config["units"],dtype=str).set_index("readgroup_name",drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

# get breed from units
BREED = units["breed"].unique()[0]

# refs
refs = ["canfam3","canfam4"]

for i in refs:
    # gen sequence groupings
    sequence_grouping(config["ref_dict"][i],i)
    # get sequence group intervals without unmapped, with unmapped, and hc caller intervals
    intervals, = glob_wildcards(os.path.join(f"results/seq_group/{i}/no_unmap","{interval}.tsv"))
    unmap_intervals, = glob_wildcards(os.path.join(f"results/seq_group/{i}/with_unmap","{interval}.tsv"))
    hc_intervals, = glob_wildcards(os.path.join(config["hc_intervals"][i],"{hc_interval}.interval_list"))

#sequence_grouping(config["ref_dict"])
## get sequence group intervals without unmapped, with unmapped, and hc caller intervals
#intervals, = glob_wildcards(os.path.join("results/seq_group/no_unmap","{interval}.tsv"))
#unmap_intervals, = glob_wildcards(os.path.join("results/seq_group/with_unmap","{interval}.tsv"))
#hc_intervals, = glob_wildcards(os.path.join(config["hc_intervals"],"{hc_interval}.interval_list"))

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
       #expand("results/merge_bam_alignment/{u.sample_name}/{ref}/{u.readgroup_name}.{ref}.merged.unsorted.bam",
       #    u=units.itertuples(),ref=refs
       #),
       ## mark duplicates
       #expand("results/mark_duplicates/{ref}/{u.sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
       #    u=units.itertuples(),ref=refs
       #),
       ## sort and fix tags
       #expand("results/sort_and_fix_tags/{ref}/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
       #    sample_name=units["sample_name"].values[0],ref=refs
       #),
       ## base recalibrator
        expand("results/base_recal/{ref}/{sample_name}.{ref}.{interval}.recal_data.csv",
            sample_name=units["sample_name"].values[0],ref=refs,
            interval=intervals
        ),
       ## gather bqsr reports
       #expand("results/gather_bqsr_reports/{u.sample_name}.{ref}.recal_data.csv",
       #    u=units.itertuples(),ref=config["ref"]
       #),
       ## apply bqsr - DUE TO THE DIFFERENCE IN INTERVALS WITH OR WITHOUT UNMAPPED,
       ## THIS IS NEEDED IN ADDITION TO THE FINAL MERGEGVCFS...FOR NOW...
       #expand("results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
       #    sample_name=units["sample_name"].values[0],ref=config["ref"],
       #    interval=unmap_intervals
       #),
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
       #expand("results/merge_gvcfs/{u.sample_name}.{ref}.g.vcf.gz",
       #    u=units.itertuples(),ref=config["ref"]
       #)
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


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        touch("results/fastqc/{sample_name}/{readgroup_name}/qc.done")
    params:
        fastqc  = config["fastqc"],
    resources:
         time   = 360,
         mem_mb = 6000, 
         cpus   = 4
    run:
        # get flowcell
        flowcell = units.loc[units["readgroup_name"] == wildcards.readgroup_name,"flowcell"].values[0]

        # only need fastqc and fastqs in one of the working/ref directories to
        # be backed up to secondary - canfam3 chosen arbitrarily
        # output path to fastqc by flowcell
        qc_outdir = os.path.join(
                config["primary"],
               #config["ref"],
                refs[0],
                BREED,
                wildcards.sample_name,
                "fastqc",
                flowcell
        )
        
        # output path to fastq by flowcell
        fq_outdir = os.path.join(
                config["primary"],
               #config["ref"],
                refs[0],
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
        tmp_dir   = f"/scratch.global/friedlab_{os.environ['USER']}/{{readgroup_name}}.ubam.tmp"
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
        ubam = "results/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    output:
        bwa_log = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{ref}/{readgroup_name}.{ref}_aligned.unmerged.bwa.stderr.log",
        bam     = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{ref}/{readgroup_name}.{ref}_aligned.unmerged.bam"
    params:
        comp_level = 5,
        java_opt   = "-Xms3000m",
        gatk       = config["gatk"],
        picard     = config["picard"],
        samtools   = config["samtools"],
        bwa        = config["bwa"],
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        ref_fasta  = lambda wildcards: config["ref_fasta"][wildcards.ref],
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
        bam  = "results/sam_to_fastq_and_bwa_mem/{sample_name}/{ref}/{readgroup_name}.{ref}_aligned.unmerged.bam"
    output:
        merged_bam = "results/merge_bam_alignment/{sample_name}/{ref}/{readgroup_name}.{ref}.merged.unsorted.bam"
    params:
        comp_level = 5,
        java_opt   = "-Xms3000m",
        gatk       = config["gatk"],
        bwa        = config["bwa"],
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        bwa_ver    = "0.7.17-r1188",
        ref_fasta  = lambda wildcards: config["ref_fasta"][wildcards.ref],
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
                        "results/merge_bam_alignment/{sample_name}/{ref}/{readgroup_name}.{ref}.merged.unsorted.bam",
                        sample_name=units["sample_name"].values[0],
                        readgroup_name=list(units["readgroup_name"]),
                        ref=refs
                      )
    output:
        dedup_bam = "results/mark_duplicates/{ref}/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
        metrics   = "results/mark_duplicates/{ref}/{sample_name}.{ref}.duplicate_metrics"
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m -Xmx16g",
        gatk       = config["gatk"],
       #tmp_dir    = "/dev/shm/{sample_name}.md.tmp"
        tmp_dir    = lambda wildcards: f"/dev/shm/{wildcards.ref}/{wildcards.sample_name}.md.tmp"
    threads: 4
    resources:
         time   = 720,
         mem_mb = 240000
    run:
        # keep only bams for the reference to be marked
        ref_bams = [i for i in input.merged_bams if wildcards.ref in i]

        # separate bams by --INPUT
        bams = " --INPUT ".join(map(str,input.ref_bams))

        # primary logs dir
        log_outdir = os.path.join(
                config["primary"],
                wildcards.ref,
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
        dedup_bam = "results/mark_duplicates/{ref}/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam"
    output:
        sorted_bam = "results/sort_and_fix_tags/{ref}/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        sorted_bai = "results/sort_and_fix_tags/{ref}/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai"
    params:
        comp_level = 5,
        java_opt   = "-Xms4000m",
        gatk       = config["gatk"],
       #ref_fasta  = config["ref_fasta"],
        ref_fasta  = lambda wildcards: config["ref_fasta"][wildcards.ref],
       #tmp_dir    = f"/scratch.global/friedlab_{os.environ['USER']}/{{sample_name}}.sort.tmp"
        tmp_dir    = lambda wildcards: f"/scratch.global/friedlab_{os.environ['USER']}/{wildcards.ref}/{wildcards.sample_name}.sort.tmp"
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
        sorted_bam = "results/sort_and_fix_tags/{ref}/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
        interval   = ancient("results/seq_group/{ref}/no_unmap/{interval}.tsv")
    output:
        recal_csv = "results/base_recal/{ref}/{sample_name}.{ref}.{interval}.recal_data.csv"
    params:
        java_opt         = "-Xms4000m",
        gatk             = config["gatk"],
        ref_fasta        = lambda wildcards: config["ref_fasta"][wildcards.ref],
        dbsnp_snp_vcf    = lambda wildcards: config["dbsnp_snp_vcf"][wildcards.ref],
        broad_snp_vcf    = lambda wildcards: config["broad_snp_vcf"][wildcards.ref],
        axelsson_snp_vcf = lambda wildcards: config["axelsson_snp_vcf"][wildcards.ref],
        dbsnp_indels_vcf = lambda wildcards: config["dbsnp_indels_vcf"][wildcards.ref],
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

#rule gather_bqsr_reports:
#    input:
#       #bqsr_reports = sorted(
#       #            expand("results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
#       #                sample_name=units["sample_name"].values[0],
#       #                ref=config["ref"],
#       #                interval=intervals
#       #            )
#       #)
#        bqsr_reports = sorted(
#                    expand("results/base_recal/{sample_name}.{ref}.{interval}.recal_data.csv",
#                        sample_name=units["sample_name"].values[0],
#                        ref=config["ref"],
#                        interval=intervals
#                        ), key=lambda item: int(os.path.basename(item).split(".")[-3].split("_")[1])
#        )
#    output:
#        report = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv"
#    params:
#        java_opt = "-Xms3000m",
#        gatk     = config["gatk"]
#    resources:
#         time   = 10,
#         mem_mb = 4000
#    run:
#        # separate reports by -I
#        reports = " -I ".join(map(str,input.bqsr_reports))
#
#        # primary logs dir
#        log_outdir = os.path.join(
#                config["primary"],
#                config["ref"],
#                BREED,
#                wildcards.sample_name,
#                "logs"
#        )
#
#        shell(f'''
#            set -e
#
#            {{params.gatk}} --java-options {{params.java_opt}} \
#                GatherBQSRReports \
#                -I {reports} \
#                -O {{output.report}}
#
#            cp {{output.report}} {log_outdir}
#        ''')
#        
#rule apply_bqsr:
#    input:
#        sorted_bam = "results/sort_and_fix_tags/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
#        report     = "results/gather_bqsr_reports/{sample_name}.{ref}.recal_data.csv",
#        interval   = ancient("results/seq_group/with_unmap/{interval}.tsv")
#    output:
#        recal_bam = "results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
#        recal_bai = "results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bai"
#    params:
#        java_opt  = "-Xms3000m",
#        gatk      = config["gatk"],
#        ref_fasta = config["ref_fasta"]
#    resources:
#         time   = 180,
#         mem_mb = 10000
#    run:
#        # separate content of interval by -L
#        with open(input.interval,"r") as f:
#            ival = f.read().strip().replace("\t", " -L ")
#
#        shell(f'''
#            {{params.gatk}} --java-options {{params.java_opt}} \
#                ApplyBQSR \
#                -R {{params.ref_fasta}} \
#                -I {{input.sorted_bam}} \
#                -O {{output.recal_bam}} \
#                -L {ival} \
#                -bqsr {{input.report}} \
#                --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
#                --add-output-sam-program-record \
#                --create-output-bam-md5 \
#                --use-original-qualities
#        ''')
#
#rule gather_bam_files:
#    input:
#        recal_bams = sorted(
#                expand("results/apply_bqsr/{sample_name}.{ref}.{interval}.aligned.duplicates_marked.recalibrated.bam",
#                    sample_name=units["sample_name"].values[0],
#                    ref=config["ref"],
#                    interval=intervals
#                    ), key=lambda item: int(os.path.basename(item).split(".")[2].split("_")[1])
#        )
#    output:
#        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam",
#        final_bai = "results/gather_bam_files/{sample_name}.{ref}.bai",
#        final_md5 = "results/gather_bam_files/{sample_name}.{ref}.bam.md5"
#    params:
#        comp_level = 5,
#        java_opt   = "-Xms2000m",
#        gatk       = config["gatk"]
#    threads: 4
#    resources:
#         time   = 180,
#         mem_mb = 16000
#    run:
#        # separate bams by -I
#        bams = " -I ".join(map(str,input.recal_bams))
#        
#        # primary bam dir
#        bam_outdir = os.path.join(
#                config["primary"],
#                config["ref"],
#                BREED,
#                wildcards.sample_name,
#                "bam"
#        )
#
#        shell(f'''
#            set -e
#
#            mkdir -p {bam_outdir}
#
#            {{params.gatk}} --java-options "-Dsamjdk.compression_level={{params.comp_level}} {{params.java_opt}}" \
#                GatherBamFiles \
#                --INPUT {bams} \
#                --OUTPUT {{output.final_bam}} \
#                --CREATE_INDEX true \
#                --CREATE_MD5_FILE true
#            
#            cp -t {bam_outdir} {{output}}
#
#        ''')
#
#rule coverage_depth_and_flagstat:
#    input:
#        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam"
#    output:
#        doc_smry = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_summary",
#        doc_stat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.depthofcoverage.sample_statistics",
#        flagstat = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat"
#    params:
#        comp_level = 5,
#        java_opt   = "-Xmx32000m",
#        gatk3      = config["gatk3"],
#        ref_fasta  = config["ref_fasta"],
#    threads: 8
#    resources:
#         time   = 360,
#         mem_mb = 32000
#    run:
#        # primary bam dir
#        bam_outdir = os.path.join(
#                config["primary"],
#                config["ref"],
#                BREED,
#                wildcards.sample_name,
#                "bam"
#        )
#        
#        shell(f'''
#            set -e
#
#            java -jar {{params.java_opt}} {{params.gatk3}} \
#                -T DepthOfCoverage \
#                -R {{params.ref_fasta}} \
#                -omitBaseOutput \
#                -omitLocusTable \
#                -omitIntervals \
#                -I {{input.final_bam}} \
#                -o results/coverage_depth_and_flagstat/{{wildcards.sample_name}}.{{wildcards.ref}}.depthofcoverage \
#                -ct 5 \
#                -ct 15 \
#                -ct 30 \
#                -nt 8
#        
#            java -jar {{params.java_opt}} {{params.gatk3}} \
#                -T FlagStat \
#                -R {{params.ref_fasta}} \
#                -I {{input.final_bam}} \
#                -o results/coverage_depth_and_flagstat/{{wildcards.sample_name}}.{{wildcards.ref}}.flagstat \
#                -nct 8
#
#            cp -t {bam_outdir} {{output}}
#        ''')
#
#rule haplotype_caller:
#    input:
#        final_bam = "results/gather_bam_files/{sample_name}.{ref}.bam",
#        interval = f"{os.path.join(config['hc_intervals'],'{hc_interval}')}.interval_list"
#    output:
#        hc_gvcf = "results/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz"
#    params:
#        java_opt  = "-Xmx10G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
#        gatk      = config["gatk"],
#        ref_fasta = config["ref_fasta"],
#    threads: 4
#    resources:
#         time   = 720,
#         mem_mb = 8000
#    shell:
#        '''
#            {params.gatk} --java-options "{params.java_opt}" \
#                HaplotypeCaller \
#                -R {params.ref_fasta} \
#                -I {input.final_bam} \
#                -L {input.interval} \
#                -O {output.hc_gvcf} \
#                -contamination 0 -ERC GVCF
#        '''
#
#rule merge_gvcfs:
#    input:
#        flagstat  = "results/coverage_depth_and_flagstat/{sample_name}.{ref}.flagstat",
#        hc_gvcfs = sorted(
#            expand(
#                "results/haplotype_caller/{hc_interval}/{sample_name}.{ref}.g.vcf.gz",
#                sample_name=units["sample_name"].values[0],
#                ref=config["ref"],
#                hc_interval=hc_intervals
#                ), key=lambda item: int(item.split("/")[-2].split("-")[0])
#        )
#    output:
#        final_gvcf     = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz",
#        final_gvcf_tbi = "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz.tbi"
#    params:
#        java_opt  = "-Xmx2000m",
#        gatk      = config["gatk"],
#        ref_fasta = config["ref_fasta"],
#    threads: 4
#    resources:
#         time   = 120,
#         mem_mb = 4000
#    run:
#        # separate bams by -I
#        gvcfs = " -I ".join(map(str,input.hc_gvcfs))
#        
#        # primary gvcf dir
#        gvcf_outdir = os.path.join(
#                config["primary"],
#                config["ref"],
#                BREED,
#                wildcards.sample_name,
#                "gvcf"
#        )
#
#        shell(f'''
#            set -e
#
#            mkdir -p {gvcf_outdir}
#
#            {params.gatk} --java-options {params.java_opt}  \
#                MergeVcfs \
#                --INPUT {gvcfs} \
#                --OUTPUT {{output.final_gvcf}}
#
#            cp -t {gvcf_outdir} {{output}}
#        ''')
#
#rule genotype_gvcfs:
#    input:
#       #ival_db  = "results/import_gvcfs/{interval}",
#       #interval = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/CanFam4/Intervals/{interval}.interval_list"
#        final_gvcf = expand(
#            "results/merge_gvcfs/{sample_name}.{ref}.g.vcf.gz",
#            sample_name=units["sample_name"].values[0],
#            ref=config["ref"]
#        ),
#       #interval   = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/MoneyCF4/Intervals/{interval}.interval_list"
#        interval   = f"{os.path.join(config['hc_intervals'],'{hc_interval}')}.interval_list"
#    output:
#        vcf = "results/genotype_gvcfs/{ref}/{hc_interval}/output.vcf.gz",
#    params:
#        gatk      = config["gatk"],
#        ref_fasta = config["ref_fasta"]
#    threads: 6
#    resources:
#         time   = 60,
#         mem_mb = 120000
#    shell:
#        '''
#            {params.gatk} --java-options "-Xmx50g -Xms50g" \
#            GenotypeGVCFs \
#                -R {params.ref_fasta} \
#                -O {output.vcf} \
#                -G StandardAnnotation \
#                --only-output-calls-starting-in-intervals \
#                -V {input.final_gvcf} \
#                -L {input.interval}
#        '''
#
#rule fltr_make_sites_only:
#    input:
#        vcf = "results/genotype_gvcfs/{ref}/{hc_interval}/output.vcf.gz"
#    output:
#        var_filtrd_vcf = "results/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.variant_filtered.vcf.gz",
#        sites_only_vcf = "results/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.sites_only.variant_filtered.vcf.gz"
#    params:
#        gatk       = config["gatk"],
#        excess_het = config["excess_het_threshold"],
#        ref_fasta  = config["ref_fasta"]
#    resources:
#         time   = 30,
#         mem_mb = 6000
#    shell:
#        '''
#            {params.gatk} --java-options "-Xmx3g -Xms3g" \
#            VariantFiltration \
#                --filter-expression "ExcessHet > {params.excess_het}" \
#                --filter-name ExcessHet \
#                -O {output.var_filtrd_vcf} \
#                -V {input.vcf}
#
#            {params.gatk} --java-options "-Xmx3g -Xms3g" \
#            MakeSitesOnlyVcf \
#                --INPUT {output.var_filtrd_vcf} \
#                --OUTPUT {output.sites_only_vcf}
#        '''
#
#rule sites_only_gather_vcf:
#    input:
#        sites_only_vcf = expand(
#            "results/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.sites_only.variant_filtered.vcf.gz",
#            ref=config["ref"],
#            hc_interval=hc_intervals
#        )
#    output:
#        gather_sites_only_vcf = "results/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
#        gather_sites_only_tbi = "results/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz.tbi"
#    params:
#        gatk   = config["gatk"],
#        picard = config["picard"]
#    threads: 4
#    resources:
#         time   = 240,
#         mem_mb = 24000
#    run:
#        vcfs = " --input ".join(map(str,input.sites_only_vcf))
#
#        shell(f'''
#            set -e
#
#            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
#                GatherVcfsCloud \
#                --ignore-safety-checks \
#                --gather-type BLOCK \
#                --input {vcfs} \
#                --output results/sites_only_gather_vcf/tmp.vcf.gz
#
#            java -jar {{params.picard}} \
#                SortVcf \
#                I=results/sites_only_gather_vcf/tmp.vcf.gz \
#                O={{output.gather_sites_only_vcf}}
#
#            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
#                IndexFeatureFile \
#                --input {{output.gather_sites_only_vcf}}
#
#            rm -f results/sites_only_gather_vcf/tmp.vcf.gz
#        ''')
#
#rule indels_var_recal:
#    input:
#        gather_sites_only_vcf = "results/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
#    output:
#        indels_recal    = "results/recal/{ref}/indels/indels.recal",
#        indels_tranches = "results/recal/{ref}/indels/indels.tranches"
#    params:
#        gatk                    = config["gatk"],
#        recal_tranche_values    = config["indel_recalibration_tranche_values"],
#        recal_annotation_values = config["indel_recalibration_annotation_values"],
#        dbsnp_indels_vcf        = config["dbsnp_indels_vcf"]
#    threads: 4
#    resources:
#         time   = 240,
#         mem_mb = 24000
#    run:
#        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
#        an_values = " -an ".join(map(str,params.recal_annotation_values))
#
#        shell(f'''
#            {{params.gatk}} --java-options "-Xmx24g -Xms24g" \
#                VariantRecalibrator \
#                -V {{input.gather_sites_only_vcf}} \
#                -O {{output.indels_recal}} \
#                --tranches-file {{output.indels_tranches}} \
#                --trust-all-polymorphic \
#                -tranche {tranche_values} \
#                -an {an_values} \
#                -mode INDEL \
#                --max-gaussians 4 \
#                --resource:dbsnp,known=false,training=true,truth=true,prior=10 {params.dbsnp_indels_vcf}
#        ''')
#
#
#rule snps_var_recal:
#    input:
#        gather_sites_only_vcf = "results/sites_only_gather_vcf/{ref}/gather.sites_only.vcf.gz",
#    output:
#        snps_recal    = "results/recal/{ref}/snps/snps.recal",
#        snps_tranches = "results/recal/{ref}/snps/snps.tranches"
#    params:
#        gatk                    = config["gatk"],
#        recal_tranche_values    = config["snp_recalibration_tranche_values"],
#        recal_annotation_values = config["snp_recalibration_annotation_values"],
#        dbsnp_snp_vcf           = config["dbsnp_snp_vcf"],
#        broad_snp_vcf           = config["broad_snp_vcf"],
#        axelsson_snp_vcf        = config["axelsson_snp_vcf"],
#        illumina_snp_vcf        = config["illumina_snp_vcf"]
#    threads: 4
#    resources:
#         time   = 240,
#         mem_mb = 16000
#    run:
#        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
#        an_values = " -an ".join(map(str,params.recal_annotation_values))
#
#        shell(f'''
#            {{params.gatk}} --java-options "-Xmx12g -Xms3g" \
#                VariantRecalibrator \
#                -V {{input.gather_sites_only_vcf}} \
#                -O {{output.snps_recal}} \
#                --tranches-file {{output.snps_tranches}} \
#                --trust-all-polymorphic \
#                -tranche {tranche_values} \
#                -an {an_values} \
#                -mode SNP \
#                --max-gaussians 6 \
#                --resource:illumina,known=true,training=true,truth=true,prior=15.0 {params.illumina_snp_vcf} \
#                --resource:broad,known=true,training=true,truth=true,prior=10.0 {params.broad_snp_vcf} \
#                --resource:axelsson,known=false,training=true,truth=true,prior=8.0 {params.axelsson_snp_vcf} \
#                --resource:dbSNP146,known=true,training=true,truth=true,prior=12.0 {params.dbsnp_snp_vcf}
#        ''')
#
#rule apply_recal:
#    input:
#        input_vcf       = "results/fltr_make_sites_only/{ref}/{hc_interval}/filtr.{hc_interval}.variant_filtered.vcf.gz",
#        indels_recal    = "results/recal/{ref}/indels/indels.recal",
#        indels_tranches = "results/recal/{ref}/indels/indels.tranches",
#        snps_recal      = "results/recal/{ref}/snps/snps.recal",
#        snps_tranches   = "results/recal/{ref}/snps/snps.tranches"
#    output:
#        recal_vcf       = "results/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz",
#        recal_vcf_index = "results/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz.tbi"
#    params:
#        gatk               = config["gatk"],
#        indel_filter_level = config["indel_filter_level"],
#        snp_filter_level   = config["snp_filter_level"]
#    resources:
#         time   = 30,
#         mem_mb = 16000
#    shell:
#        '''
#            set -e
#
#            mkdir -p results/apply_recal/{wildcards.hc_interval}/
#
#            {params.gatk} --java-options "-Xmx15g -Xms5g" \
#                ApplyVQSR \
#                -O results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf \
#                -V {input.input_vcf} \
#                --recal-file {input.indels_recal} \
#                --tranches-file {input.indels_tranches} \
#                --truth-sensitivity-filter-level {params.indel_filter_level} \
#                --create-output-variant-index true \
#                -mode INDEL
#
#            {params.gatk} --java-options "-Xmx15g -Xms5g" \
#                ApplyVQSR \
#                -O {output.recal_vcf} \
#                -V results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf \
#                --recal-file {input.snps_recal} \
#                --tranches-file {input.snps_tranches} \
#                --truth-sensitivity-filter-level {params.snp_filter_level} \
#                --create-output-variant-index true \
#                -mode SNP
#
#            rm -f results/apply_recal/{wildcards.hc_interval}/tmp.indel.recalibrated.vcf
#        '''
#
#rule final_gather_vcfs:
#    input:
#        recal_vcfs = sorted(
#            expand(
#                "results/apply_recal/{ref}/{hc_interval}/recal.{hc_interval}.vcf.gz", 
#                ref=config["ref"],
#                hc_interval=hc_intervals
#            )
#        )
#    output:
#        final_vcf       = "results/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
#        final_vcf_index = "results/final_gather_vcfs/joint_genotype.{ref}.vcf.gz.tbi"
#    params:
#        gatk = config["gatk"],
#    threads: 4
#    resources:
#         time   = 240,
#         mem_mb = 22000
#    run:
#        vcfs = " --input ".join(map(str,input.recal_vcfs))
#
#        shell(f'''
#            set -e
#            set -o pipefail
#
#            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
#                GatherVcfsCloud \
#                --ignore-safety-checks \
#                --gather-type BLOCK \
#                --input {vcfs} \
#                --output {{output.final_vcf}}
#
#            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
#                IndexFeatureFile \
#                --input {{output.final_vcf}}
#        ''')
#
#rule collect_metrics_on_vcf:
#    input:
#        final_vcf       = "results/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
#        final_vcf_index = "results/final_gather_vcfs/joint_genotype.{ref}.vcf.gz.tbi"
#    output:
#        detail_metrics  = "results/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_detail_metrics",
#        summary_metrics = "results/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_summary_metrics"
#    params:
#        gatk               = config["gatk"],
#        dbsnp_snp_vcf      = config["dbsnp_snp_vcf"],
#        ref_dict           = config["ref_dict"],
#        eval_interval_list = config["eval_interval_list"],
#        metrics_prefix     = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}"
#    threads: 8
#    resources:
#         time   = 360,
#         mem_mb = 22000
#    shell:
#        '''
#            set -e
#
#            mkdir -p results/collect_metrics_on_vcf/
#
#            {params.gatk} --java-options "-Xmx18g -Xms6g" \
#                CollectVariantCallingMetrics \
#                --INPUT {input.final_vcf} \
#                --DBSNP {params.dbsnp_snp_vcf} \
#                --SEQUENCE_DICTIONARY {params.ref_dict} \
#                --OUTPUT {params.metrics_prefix} \
#                --THREAD_COUNT 8 \
#                --TARGET_INTERVALS {params.eval_interval_list}
#        '''
#
#rule vep_final_vcf:
#    input:
#        final_vcf = "results/final_gather_vcfs/joint_genotype.{ref}.vcf.gz",
#    output:
#        vep_vcf     = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
#        vep_vcf_tbi = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz.tbi"
#    params:
#        conda_vep = config["conda_vep"],
#        ref_fasta = config["ref_fasta"],
#        ref_gtf   = config["ref_gtf"]
#    threads: 12
#    resources:
#         time   = 4320,
#         mem_mb = 60000
#    run:
#        import os
#        out_name = os.path.splitext(output.vep_vcf)[0] 
#        
#        # if canfam4 VEP with the gtf
#        if "canfam4" in wildcards.ref:
#            shell(f'''
#                set +eu
#
#                eval "$(conda shell.bash hook)"
#                conda activate {params.conda_vep}
#
#                vep \
#                    -i {{input.final_vcf}} \
#                    -o {out_name} \
#                    --gtf {{params.ref_gtf}} \
#                    --fasta {{params.ref_fasta}} \
#                    --everything \
#                    --force_overwrite \
#                    --vcf \
#                    --dont_skip
#
#                bgzip --threads 12 -c {out_name} > {{output.vep_vcf}} &&
#                tabix -p vcf {{output.vep_vcf}}
#            ''')
#        # else VEP with the canfam3 annotation from ensembl
#        else:
#            shell(f'''
#                set +eu
#
#                eval "$(conda shell.bash hook)"
#                conda activate {params.conda_vep}
#
#                vep \
#                  -i {{input.final_vcf}} \
#                  --cache \
#                  --everything \
#                  -o {out_name} \
#                  --vcf \
#                  --no_stats \
#                  --species=canis_familiaris \
#                  --offline \
#                  --dont_skip
#
#                bgzip --threads 12 -c {out_name} > {{output.vep_vcf}} &&
#                tabix -p vcf {{output.vep_vcf}}
#            ''')
#
################################################################################
###### ONLY NEEDED TO GENERATE THE COMMON VARIANTS ONE WHEN A NEW POP VCF ######
###### IS CREATED ##############################################################
################################################################################
#
## ADD A CHECK TO ONLY CREATE THE INDEX IF NECESSARY???
#rule select_common_vars:
#    input:
#        pop_vcf = config["pop_vcf"],
#    output:
#       #common_vars     = "results/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz",
#       #common_vars_tbi = "results/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz.tbi"
#        common_vars = config["common_vars"]
#    params:
#        conda_vep = config["conda_vep"],
#        ref_fasta = config["ref_fasta"],
#    threads: 8
#    resources:
#         time   = 2160,
#         mem_mb = 32000
#    shell:
#        '''
#            set +eu
#
#            eval "$(conda shell.bash hook)"
#            conda activate {params.conda_vep}
#            
#            module load bcftools
#
#            bcftools view -Oz \
#                -i 'AF[*]>0.005' \
#                {input.pop_vcf} \
#                -o {output.common_vars}
#
#            tabix -p vcf {output.common_vars}
#        '''
#
#
#rule select_variants_to_table:
#    input:            
#        final_vcf      = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
#        common_vars    = "results/select_common_vars/joint_genotype.{ref}.vep.af_nonmajor.vcf.gz",
#        detail_metrics = "results/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_detail_metrics",
#    output:
#        unique_vars           = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf",
#        rare_and_common_vars  = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_and_common_vars.vcf",
#        rare_vars             = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf",
#        unique_vars_vep_split = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vep_split.txt",
#        rare_vars_vep_split   = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vep_split.txt"
#    params:
#        gatk = config["gatk"],
#        tmp_dir = "/dev/shm/money_vars",
#        ref_fasta = config["ref_fasta"],
#        pop_vcf = config["pop_vcf"]
#    threads: 4
#    resources:
#         time   = 2880,
#         mem_mb = 16000
#    shell:
#        '''
#            set -e
#    
#            module load bcftools    
#            mkdir -p {params.tmp_dir} results/select_vars_to_table
#
#            # all unique variants
#            {params.gatk} SelectVariants \
#                -R {params.ref_fasta} \
#                -V {input.final_vcf} \
#                -disc {params.pop_vcf} \
#                --tmp-dir {params.tmp_dir} \
#                -O {output.unique_vars}
#
#            # all common and  rare variants
#            {params.gatk} SelectVariants \
#                -R {params.ref_fasta} \
#                -V {input.final_vcf} \
#                -disc {input.common_vars} \
#                -O {output.rare_and_common_vars}
#
#            # rare only variants
#            {params.gatk} SelectVariants \
#                -R {params.ref_fasta} \
#                -V {output.rare_and_common_vars} \
#                -disc {output.unique_vars} \
#                -O {output.rare_vars}
#
#            # vars to table - all unique
#            bcftools +split-vep \
#                {output.unique_vars} \
#                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
#                -d -A tab \
#                -o {output.unique_vars_vep_split}
#
#            # vars to table - rare only
#            bcftools +split-vep \
#                {output.rare_vars} \
#                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
#                -d -A tab \
#                -o {output.rare_vars_vep_split}
#        '''
#
#rule bgzip_and_tabix:
#    input:            
#        unique_vars = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf",
#        rare_vars   = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf",
#    output:
#        unique_vars_gz  = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz",
#        unique_vars_tbi = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz.tbi",
#        rare_vars_gz    = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz",
#        rare_vars_tbi   = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz.tbi",
#    params:
#        conda_vep = config["conda_vep"],
#    threads: 8
#    resources:
#         time   = 60,
#         mem_mb = 16000
#    shell:
#        '''
#            set +eu
#
#            eval "$(conda shell.bash hook)"
#            conda activate {params.conda_vep}
#            
#            module load bcftools
#
#            # bgzip and tabix the unique vcf
#            bgzip -c {input.unique_vars} > {output.unique_vars_gz}
#            tabix -p vcf {output.unique_vars_gz}
#        
#            # and the rare vcf
#            bgzip -c {input.rare_vars} > {output.rare_vars_gz}
#            tabix -p vcf {output.rare_vars_gz}
#        '''
#
#rule final_output:
#    input:
#       #unique_vars_vep_split = f"results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.unique_vars.vep_split.txt",
#       #rare_vars_vep_split   = f"results/select_vars_to_table/{units['sample_name'].values[0]}.{config['ref']}.rare_vars.vep_split.txt"
#        unique_vars_vep_split = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vep_split.txt",
#        rare_vars_vep_split   = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vep_split.txt"
#    output:
#       #excel_sheet = f"results/final_output/{units['sample_name'].values[0]}_{config['ref']}_tables.xlsx",
#        excel_sheet = "results/final_output/{ref}/{sample_out}_{ref}_tables.xlsx",
#       #excel_sheet = "{sample_out}_{ref}_tables.xlsx",
#    params:
#        conda_vep = config["conda_vep"],
#        base_name = "results/final_output/{ref}"
#    threads: 4
#    resources:
#         time   = 30,
#         mem_mb = 16000
#    run:
#       #import os
#       #import pandas as pd
#        
#        # first add header to the rare and unique tables
#       #header = [
#       #  'chrom', 'pos', 'ref', 'alt', 'ac', 
#       #  'allele', 'consequence', 'impact', 'symbol', 
#       #  'gene', 'feature_type', 'feature', 'biotype', 
#       #  'exon', 'intron', 'hgvsc', 'hgvsp', 
#       #  'cdna_position', 'cds_position', 'protein_position', 'amino_acids', 
#       #  'codons', 'existing_variation', 'distance', 'strand', 
#       #  'flags', 'variant_class', 'symbol_source', 'hgnc_id', 
#       #  'canonical','mane', 'tsl', 'appris', 
#       #  'ccds', 'ensp', 'swissprot', 'trembl',
#       #  'uniparc', 'gene_pheno', 'sift', 'domains', 
#       #  'mirna', 'af', 'afr_af', 'amr_af', 
#       #  'eas_af', 'eur_af', 'sas_af', 'aa_af', 
#       #  'ea_af', 'gnomad_af', 'gnomad_afr_af', 'gnomad_amr_af', 
#       #  'gnomad_asj_af', 'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_nfe_af', 
#       #  'gnomad_oth_af', 'gnomad_sas_af', 'max_af', 'max_af_pops', 
#       #  'clin_sig', 'somatic', 'pheno', 'pubmed', 
#       #  'check_ref', 'motif_name', 'motif_pos', 'high_inf_pos',
#       #  'motif_score_change'
#       #]
#        header = [
#            'chrom', 'pos', 'ref', 'alt', 'ac',
#            'allele', 'consequence', 'impact', 'symbol', 'gene', 
#            'feature_type', 'feature', 'biotype', 'exon', 'intron', 
#            'hgvsc', 'hgvsp', 'cdna_position', 'cds_position', 
#            'protein_position', 'amino_acids', 'codons', 'existing_variation',
#            'distance', 'strand', 'flags', 'variant_class', 'symbol_source', 
#            'hgnc_id', 'canonical', 'mane_select', 'mane_plus_clinical', 
#            'tsl', 'appris', 'ccds', 'ensp', 'swissprot', 'trembl', 'uniparc', 
#            'uniprot_isoform', 'source', 'gene_pheno', 'sift', 'polyphen', 
#            'domains', 'mirna', 'hgvs_offset', 'af', 'afr_af', 'amr_af', 
#            'eas_af', 'eur_af', 'sas_af', 'aa_af', 'ea_af', 'gnomad_af', 
#            'gnomad_afr_af', 'gnomad_amr_af', 'gnomad_asj_af', 
#            'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_nfe_af', 
#            'gnomad_oth_af', 'gnomad_sas_af', 'max_af', 'max_af_pops', 
#            'clin_sig', 'somatic', 'pheno', 'pubmed', 'check_ref', 
#            'motif_name', 'motif_pos', 'high_inf_pos', 
#            'motif_score_change', 'transcription_factors', 
#            'canFam4_genes_ncbi.gtf.gz'
#        ]
#     
#        # for each table, add header and print the rest of table   
#        dfs = {}
#        for f in [input.unique_vars_vep_split, input.rare_vars_vep_split]:
#            tmp = os.path.basename(f).replace(".txt",".reform.txt")
#            reform = os.path.join(params.base_name,tmp)
#            # get name of each table as key for dict
#            name = tmp.split(".vep_split")[0]
#            # open table, add header, and write to reform 
#            with open(f, "r") as infile, open(reform, "w") as outfile:
#                print("\t".join(header), file=outfile)
#                for line in infile:
#                    print(line.strip(), file=outfile)
#            # read fixed table to dataframe
#            dfs[name] = pd.read_csv(reform, sep="\t")
#
#        # write each df to the same excel sheet
#        with pd.ExcelWriter(
#                output.excel_sheet,
#                engine="xlsxwriter",
#                options={"strings_to_formulas": False}
#            ) as writer:
#            for k,v in dfs.items():
#                v.to_excel(writer,
#                           sheet_name=k,
#                           index=False)
#            writer.save()
#
#rule manifest_and_upload:
#    input:
#        vep_vcf         = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
#        vep_vcf_tbi     = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz.tbi",
#        unique_vars_gz  = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz",
#        unique_vars_tbi = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.unique_vars.vcf.gz.tbi",
#        rare_vars_gz    = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz",
#        rare_vars_tbi   = "results/select_vars_to_table/{ref}/{sample_out}.{ref}.rare_vars.vcf.gz.tbi",
#        excel_sheet     = "results/final_output/{ref}/{sample_out}_{ref}_tables.xlsx",
#    output:
#       #excel_sheet = f"results/final_output/{units['sample_name'].values[0]}_{config['ref']}_tables.xlsx",
#        manifest      = "results/{sample_out}_{ref}_manifest.txt",
#        money_tar_gz  = "results/final_output/{ref}/{sample_out}_{ref}_K9MM.tar.gz",
#       #excel_sheet = "{sample_out}_{ref}_tables.xlsx",
#    params:
#       #base_dir = "results/final_output/{ref}/upload/"
#    run:
#       #os.makedirs(os.path.dirname(output.manifest),exist_ok=True)
#        
#        with open(output.manifest, "w") as outfile:
#            s = f'''
#                Included files:
#                1) results/vep_final_vcf/joint_genotype.{wildcards.ref}.vep.vcf.gz - annotated VCF and index (.tbi) for all {wildcards.sample_out} variants
#                2) results/select_vars_to_table/{wildcards.sample_out}.{wildcards.ref}.unique_vars.vcf.gz - annotated VCF and index (.tbi) for all unique {wildcards.sample_out} variants
#                3) results/select_vars_to_table/{wildcards.sample_out}.{wildcards.ref}.rare_vars.vcf.gz - annotated VCF and index (.tbi) for all rare {wildcards.sample_out} variants
#                4) results/final_output/{wildcards.sample_out}_{wildcards.ref}_tables.xlsx - contains two sheets
#                \t i) {wildcards.sample_out}_{wildcards.ref}.unique_vars - table version of item 2
#                \t ii) {wildcards.sample_out}_{wildcards.ref}.rare_vars - table version of item 3
#            '''
#            print(s, file=outfile)
#
#        shell('''
#            tar -czvf {output.money_tar_gz} \
#                {input.vep_vcf} \
#                {input.vep_vcf_tbi} \
#                {input.unique_vars_gz} \
#                {input.unique_vars_tbi} \
#                {input.rare_vars_gz} \
#                {input.rare_vars_tbi} \
#                {input.excel_sheet} \
#                {output.manifest}
#        ''')
#
