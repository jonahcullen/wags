import pandas as pd

def get_fastq(wildcards):
    '''
    Get fastq files of given sample-unit.
    '''
    fastqs = units.loc[(wildcards.readgroup_name), ["fastq_1", "fastq_2"]].dropna()
    
    return {"r1": fastqs.fastq_1, "r2": fastqs.fastq_2}


#units = pd.read_table(config["units"], dtype=str).set_index(["sample","time"], drop=False)
units = pd.read_table(config["units"],dtype=str).set_index("readgroup_name",drop=False)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
print(units)

rule all:
    input:
       ## fastqc 
       #expand("results/{u.sample_name}/{u.readgroup_name}/{u.platform_unit}/hello",
       #    u=units.itertuples()
       #),
       ## fastq to ubam
       #expand("results/fastqs_to_ubam/{u.sample_name}/{u.readgroup_name}.unmapped.bam",
       #    u=units.itertuples()
       #)
       ## sam_to_fastq_and_bwa_mem
       #expand("results/sam_to_fastq_and_bwa_mem/{u.sample_name}/{u.readgroup_name}.aligned.unmerged.bam",
       #    u=units.itertuples()
       #),
       ## merge bams
        expand("results/merge_bam_alignment/{u.sample_name}/{u.readgroup_name}.merged.unsorted.bam",
            u=units.itertuples()
        )


#expand("results/step01/{breed}/{breed}_snp_infor.Rdat",
#       #        breed=units.breed.unique()

#rule test:
#    output:
#        hello = "results/{sample_name}/{platform_unit}/hello"
#    shell:
#        '''
#            touch {output.hello}        
#        '''

rule fastqc:
    input:
        unpack(get_fastq)
    params:
        fastqc = config["fastqc"],
        outdir = "results/{sample_name}/{readgroup_name}/{platform_unit}/"
    resources:
         time   = 120,
         mem_mb = 6000, 
         cpus   = 1
    output:
        "results/{sample_name}/{readgroup_name}/{platform_unit}/hello"
    shell:
        '''
            {params.fastqc} --outdir {params.outdir} \
                {input}
        '''

rule fastqs_to_ubam:
    input:
        unpack(get_fastq)
    output:
        ubam = "results/fastqs_to_ubam/{sample_name}/{readgroup_name}.unmapped.bam"
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"]
    resources:
         time   = 360,
         mem_mb = 12000, 
         cpus   = 6
    run:

        # get fastq meta data for logging in the bam
        library_name = units.loc[units['readgroup_name'] == wildcards.readgroup_name,'library_name'].values[0]
        platform_unit = units.loc[units['readgroup_name'] == wildcards.readgroup_name,'platform_unit'].values[0]
        run_date = units.loc[units['readgroup_name'] == wildcards.readgroup_name,'run_date'].values[0]
        platform_name = units.loc[units['readgroup_name'] == wildcards.readgroup_name,'platform_name'].values[0]
        sequencing_center = units.loc[units['readgroup_name'] == wildcards.readgroup_name,'sequencing_center'].values[0]

        f'''
            {{params.gatk}} --java-options "-Xmx6000m" \
                FastqToSam \
                --FASTQ {{input.r1}} \
                --FASTQ2 {{input.r2}} \
                --OUTPUT {output.ubam} \
                --READ_GROUP_NAME {{wildcards.readgroup_name}} \
                --SAMPLE_NAME {{wildcards.sample_name}} \
                --LIBRARY_NAME {library_name} \
                --PLATFORM_UNIT {platform_unit} \
                --RUN_DATE {run_date} \
                --PLATFORM {platform_name} \
                --SEQUENCING_CENTER {sequencing_center}
        '''

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
        ref_fasta  = config["ref_fasta"]
    resources:
         time   = 1440,
         mem_mb = 32000, 
         cpus   = 16
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
        merged_bam = "results/merge_bam_alignment/{sample_name}/{readgroup_name}.merged.unsorted.bam"
    params:
        comp_level = 5,
        java_opt   = "-Xms3000m",
        gatk       = config["gatk"],
        bwa        = config["bwa"],
        bwa_cl     = "mem -K 100000000 -p -v 3 -t 16 -Y",
        ref_fasta  = config["ref_fasta"]
    resources:
         time   = 1440,
         mem_mb = 32000, 
         cpus   = 16
    shell:
        '''
            # set the bash variable needed for the command-line
            bash_ref_fasta={params.ref_fasta}

            # get bwa version
            BWA_VER=$(/panfs/roc/groups/0/fried255/shared/gatk4_workflow/tools/bwa-0.7.17/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //' 2>&1)

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
                --PROGRAM_GROUP_VERSION "${BWA_VER}" \
                --PROGRAM_GROUP_COMMAND_LINE "{params.bwa_cl}" \
                --PROGRAM_GROUP_NAME "bwamem" \
                --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
                --ALIGNER_PROPER_PAIR_FLAGS true \
                --UNMAP_CONTAMINANT_READS true
        '''

#  call GetBwaVersion {
#    input: 
#      bwa = bwa
#  }
#
#  # Align flowcell-level unmapped input bams in parallel
#  scatter (unmapped_bam in read_lines(CreateFoFN.fofn_list)) {
#
#    # Get the basename, i.e. strip the filepath and the extension
#    String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)
#
#    # Map reads to reference
#    call SamToFastqAndBwaMem {
#      input:
#        bwa = bwa,
#        picard = picard,
#        samtools = samtools,
#        input_bam = unmapped_bam,
#        bwa_commandline = bwa_commandline,
#        output_bam_basename = bam_basename + ".unmerged",
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
#        ref_dict = ref_dict,
#        compression_level = compression_level
#     }
#
#    # Merge original uBAM and BWA-aligned BAM 
#    call MergeBamAlignment {
#      input:
#        gatk = gatk,
#        unmapped_bam = unmapped_bam,
#        bwa_commandline = bwa_commandline,
#        bwa_version = GetBwaVersion.version,
#        aligned_bam = SamToFastqAndBwaMem.output_bam,
#        output_bam_basename = bam_basename + ".aligned.unsorted",
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
#        ref_dict = ref_dict,
#        compression_level = compression_level
#    }
#  }
#task SamToFastqAndBwaMem {
#  String input_bam
#  String bwa_commandline
#  String output_bam_basename
#  String ref_fasta
#  String ref_fasta_index
#  String ref_dict
#
#  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
#  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
#  # references such as b37 and hg19.
#  String ref_amb
#  String ref_ann
#  String ref_bwt
#  String ref_pac
#  String ref_sa
#
#  Int compression_level
#
#  String bwa
#  String picard
#  String samtools
#  String java_opt
#
#  Int mem_gb
#  Int cpu
#  String walltime
#
#  command <<<
#    set -o pipefail
#    set -e
#
#    # set the bash variable needed for the command-line
#    bash_ref_fasta=${ref_fasta}
#
#		java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${picard} \
#      SamToFastq \
#			INPUT=${input_bam} \
#			FASTQ=/dev/stdout \
#			INTERLEAVE=true \
#			NON_PF=true \
#    | \
#		${bwa}${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
#    | \
#		${samtools} view -1 - > ${output_bam_basename}.bam
#
#  >>>
#  runtime {
#    backend: "slurm"
#    mem_gb: mem_gb
#    cpu: cpu
#    walltime: walltime
#  }
#  output {
#    File output_bam = "${output_bam_basename}.bam"
#    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
#  }
#}
#
## Merge original input uBAM file with BWA-aligned BAM file
#task MergeBamAlignment {
#  String unmapped_bam
#  String bwa_commandline
#  String bwa_version
#  String aligned_bam
#  String output_bam_basename
#  String ref_fasta
#  String ref_fasta_index
#  String ref_dict
#
#  Int compression_level
#
#  String gatk
#  String java_opt
#
#  Int mem_gb
#  Int cpu
#  String walltime
#
#  command {
#    # set the bash variable needed for the command-line
#    bash_ref_fasta=${ref_fasta}
#    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
#      MergeBamAlignment \
#      --VALIDATION_STRINGENCY SILENT \
#      --EXPECTED_ORIENTATIONS FR \
#      --ATTRIBUTES_TO_RETAIN X0 \
#      --ALIGNED_BAM ${aligned_bam} \
#      --UNMAPPED_BAM ${unmapped_bam} \
#      --OUTPUT ${output_bam_basename}.bam \
#      --REFERENCE_SEQUENCE ${ref_fasta} \
#      --PAIRED_RUN true \
#      --SORT_ORDER "unsorted" \
#      --IS_BISULFITE_SEQUENCE false \
#      --ALIGNED_READS_ONLY false \
#      --CLIP_ADAPTERS false \
#      --MAX_RECORDS_IN_RAM 2000000 \
#      --ADD_MATE_CIGAR true \
#      --MAX_INSERTIONS_OR_DELETIONS -1 \
#      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
#      --PROGRAM_RECORD_ID "bwamem" \
#      --PROGRAM_GROUP_VERSION "${bwa_version}" \
#      --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
#      --PROGRAM_GROUP_NAME "bwamem" \
#      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
#      --ALIGNER_PROPER_PAIR_FLAGS true \
#      --UNMAP_CONTAMINANT_READS true
#  }
#  runtime {
#    backend: "slurm"
#    mem_gb: mem_gb
#    cpu: cpu
#    walltime: walltime
#  }
#  output {
#    File output_bam = "${output_bam_basename}.bam"
#  }
#}
