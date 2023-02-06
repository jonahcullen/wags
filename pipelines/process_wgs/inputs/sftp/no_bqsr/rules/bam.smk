
rule fastqs_to_ubam:
    input:
        unpack(get_fastq),
    output:
        ubam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.unmapped.bam"
    params:
        java_opt          = "-Xms6000m",
        ref_fasta         = config['ref_fasta'],
       #tmp_dir           = f"/dev/shm/{os.environ['USER']}/{{readgroup_name}}_{config['ref']}/",
        tmp_dir           = config['tmp_dir']['fastq_tmp'],
        library_name      = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"library_name"].values[0],
        platform_unit     = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"platform_unit"].values[0],
        run_date          = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"run_date"].values[0],
        platform_name     = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"platform_name"].values[0],
        sequencing_center = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"sequencing_center"].values[0],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.fastqs_to_ubam.benchmark.txt"
    threads: 4
    resources:
         time   = 400,
         mem_mb = lambda wildcards, attempt: 2**(attempt-1)*60000,
    shell:
        '''
            mkdir -p {params.tmp_dir}

            gatk --java-options {params.java_opt} \
                FastqToSam \
                --TMP_DIR {params.tmp_dir} \
                --FASTQ {input.r1} \
                --FASTQ2 {input.r2} \
                --OUTPUT {output.ubam} \
                --READ_GROUP_NAME {wildcards.readgroup_name} \
                --SAMPLE_NAME {wildcards.sample_name} \
                --LIBRARY_NAME {params.library_name} \
                --PLATFORM_UNIT {params.platform_unit} \
                --RUN_DATE {params.run_date} \
                --PLATFORM {params.platform_name} \
                --SEQUENCING_CENTER {params.sequencing_center}
        '''

rule mark_adapters:
    input:
        ubam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.unmapped.bam"
    output:
        mark_adapt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.mark_adapt.unmapped.bam",
        metrics    = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.mark_adapt.metrics.txt"
    params:
        tmp_dir = f"/dev/shm/{os.environ['USER']}/{{readgroup_name}}_mark_adapt/"
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.mark_adapters.benchmark.txt"
    threads: 4
    resources:
         time   = 800,
         mem_mb = lambda wildcards, attempt: 2**(attempt-1)*20000,
    shell:
        '''
            mkdir -p {params.tmp_dir}
            
            java -Xmx8G -jar /opt/wags/src/picard.jar  \
                MarkIlluminaAdapters \
                --INPUT {input.ubam} \
                --OUTPUT {output.mark_adapt} \
                --METRICS {output.metrics} \
                --TMP_DIR {params.tmp_dir}
        '''

rule sam_to_fastq_and_bwa_mem:
    input:
        mark_adapt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.mark_adapt.unmapped.bam",
    output:
        bwa_log = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.{ref}_aligned.unmerged.bwa.stderr.log",
        bam     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.{ref}_aligned.unmerged.bam"
    params:
        tmp_dir   = f"/dev/shm/{os.environ['USER']}/{{readgroup_name}}_sam_fastq/",
        java_opt  = "-Xms3000m",
        bwa_cl    = "mem -K 100000000 -p -v 3 -t 16 -Y",
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.sam_to_fq_bwa.benchmark.txt"
    threads: 16
    resources:
         time   = 1440,
         mem_mb = lambda wildcards, attempt: 2**(attempt-1)*40000,
    shell:
        '''
            set -o pipefail
            set -e

            mkdir -p {params.tmp_dir}

            # set the bash variable needed for the command-line
            bash_ref_fasta={params.ref_fasta}

		    java -Dsamjdk.compression_level=5 {params.java_opt} -jar /opt/wags/src/picard.jar \
                SamToFastq \
                I={input.mark_adapt} \
                FASTQ=/dev/stdout \
                CLIPPING_ATTRIBUTE=XT \
                CLIPPING_ACTION=2 \
                INTERLEAVE=true \
                NON_PF=true \
                TMP_DIR={params.tmp_dir} \
            | \
            bwa {params.bwa_cl} $bash_ref_fasta /dev/stdin -  2> >(tee {output.bwa_log} >&2) \
            | \
            samtools view -1 - > {output.bam}
        '''

# get bwa version
#BWA_VER=$(bwa 2>&1 | rg -e '^Version' | sed 's/Version: //' 2>&1)

rule merge_bam_alignment:
    input:
        ubam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.unmapped.bam",
        bam  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.{ref}_aligned.unmerged.bam"
    output:
        merged_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.{ref}.merged.unsorted.bam"
    params:
        java_opt   = "-Xms3000m",
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        bwa_ver    = "0.7.17-r1188",
        ref_fasta  = config['ref_fasta']
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.merge_bams.benchmark.txt"
    threads: 4
    resources:
        time   = 720,
        mem_mb = lambda wildcards, attempt: 2**(attempt-1)*12000,
    shell:
        '''
            gatk --java-options "-Dsamjdk.compression_level=5 {params.java_opt}" \
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
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{readgroup_name}.{ref}.merged.unsorted.bam",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'], 
            readgroup_name=list(units['readgroup_name']),
        )
    output:
        dedup_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
        metrics   = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.duplicate_metrics",
    params:
        bams       = lambda wildcards, input: " --INPUT ".join(map(str,input.merged_bams)),
        java_opt   = "-Xms4000m -Xmx16g",
        tmp_dir    = "/dev/shm/{sample_name}_{ref}.md.tmp"
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.mark_duplicates.benchmark.txt"
    threads: 4
    resources:
         time   = 720,
         mem_mb = lambda wildcards, attempt: 2**(attempt-1)*60000,
    shell:
        '''
            mkdir -p {params.tmp_dir}

            gatk --java-options "-Dsamjdk.compression_level=5 {params.java_opt}" \
                MarkDuplicates \
                --TMP_DIR {params.tmp_dir} \
                --INPUT {params.bams} \
                --OUTPUT {output.dedup_bam} \
                --METRICS_FILE {output.metrics} \
                --VALIDATION_STRINGENCY SILENT \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --ASSUME_SORT_ORDER "queryname" \
                --CREATE_MD5_FILE true
        '''

rule sort_and_fix_tags:
    input:
        dedup_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.unsorted.duplicates_marked.bam",
    output:
        sorted_bam = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam"),
        sorted_bai = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai"),
        sorted_md5 = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam.md5"),
    params:
        java_opt  = "-Xms4000m",
        tmp_dir   = config['tmp_dir']['sort_tmp'],
        ref_fasta = config['ref_fasta']
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.sort_and_fix.benchmark.txt"
    threads: 12
    resources:
         time   = 720,
         mem_mb = lambda wildcards, attempt: 2**(attempt-1)*24000,
    shell:
        '''
            set -o pipefail

            gatk --java-options "-Dsamjdk.compression_level=5 {params.java_opt}" \
                SortSam \
                --TMP_DIR {params.tmp_dir} \
                --INPUT {input.dedup_bam} \
                --OUTPUT /dev/stdout \
                --SORT_ORDER "coordinate" \
                --CREATE_INDEX false \
                --CREATE_MD5_FILE false \
            | \
            gatk --java-options "-Dsamjdk.compression_level=5 {params.java_opt}" \
                SetNmMdAndUqTags \
                --INPUT /dev/stdin \
                --OUTPUT {output.sorted_bam} \
                --CREATE_INDEX true \
                --CREATE_MD5_FILE true \
                --REFERENCE_SEQUENCE {params.ref_fasta}
       '''

# optional argument to left align bam
if config['left_align']:
    rule left_align_bam:
        input:
            sorted_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam",
            sorted_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai"
        output:
            left_bam = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam"),
            left_bai = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai"),
            left_md5 = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam.md5"),
        params:
            java_opt  = "-Xms4000m",
            ref_fasta = config['ref_fasta']
        benchmark:
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.sort_fix_left_align.benchmark.txt"
        threads: 4
        resources:
             time   = 480,
             mem_mb = lambda wildcards, attempt: 2**(attempt-1)*24000,
        shell:
            '''
                gatk --java-options "-Dsamjdk.compression_level=5 {params.java_opt}" \
                    LeftAlignIndels \
                    -I {input.sorted_bam} \
                    -O {output.left_bam} \
                    --create-output-bam-index true \
                    --create-output-bam-md5 true \
                    -R {params.ref_fasta}
           '''

