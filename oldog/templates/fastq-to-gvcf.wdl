##Copyright Broad Institute, 2018
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
## ```readgroup   fastq_pair1_file_path   fastq_pair2_file_path   sample_name   library_name   platform_unit   run_date   platform_name   sequecing_center``` 
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
## Cromwell version support 
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.
##
##
##
## Copyright Broad Institute, 2019
## 
## This WDL pipeline implements data pre-processing according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements 
## - GATK 4 or later
## - BWA 0.7.15-r1140
## - Picard 2.16.0-SNAPSHOT
## - Samtools 1.3.1 (using htslib 1.3.1)
## - Python 2.7
##
## Cromwell version support 
## - Successfully tested on v37
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.
##
##
##
## Copyright Broad Institute, 2019
## 
## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
## When executed the workflow scatters the HaplotypeCaller tool over a sample
## using an intervals list file. The output file produced will be a
## single gvcf file which can be used by the joint-discovery workflow.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support 
## - Successfully tested on v37
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToGvcf_GATK4 {

  # fastq specific input
  Array[String] sample_name 
  Array[String] fastq_1 
  Array[String] fastq_2 
  Array[String] readgroup_name 
  Array[String] library_name 
  Array[String] platform_unit 
  Array[String] run_date 
  Array[String] platform_name 
  Array[String] sequencing_center 
  Array[String] flowcell

  String primary_path
  String ubam_list_name
  String out_path  

  String sample
  String ref_name

  String unmapped_bam_suffix
  String base_file_name = sample + "." + ref_name
 
  # reference input 
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbsnp_snp_vcf
  File dbsnp_snp_vcf_index
  File broad_snp_vcf
  File broad_snp_vcf_index
  File axelsson_snp_vcf
  File axelsson_snp_vcf_index   
  File dbsnp_indels_vcf
  File dbsnp_indels_vcf_index

  String bwa_commandline = " mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
  Int compression_level 

  Int flowcell_small_disk
  Int flowcell_medium_disk
  Int agg_small_disk
  Int agg_medium_disk
  Int agg_large_disk
  
  File input_bam
  File input_bam_index
  File scattered_calling_intervals_list
  
  Boolean? make_gvcf
  Boolean making_gvcf = select_first([make_gvcf,true])
 
  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  #is the input a cram file?
  Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"

 #String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
 #String vcf_basename = sample_basename
  String output_suffix = if making_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_filename = base_file_name + output_suffix

  # Tools
  String samtools
  String gatk
  String gatk3
  String bwa
  String picard

###############################################################################
#
# Task calls for converting paired fastq files to uBAMs
#
###############################################################################

  # Convert multiple pairs of input fastqs in parallel
  scatter (i in range(length(readgroup_name))) {

    # Get fastqc output
    call RunFastQC {
      input:
        target = primary_path + "/fastqc/" + flowcell[i],
        r1 = fastq_1[i],
        r2 = fastq_2[i]
    }

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        gatk = gatk,
        sample_name = sample_name[i],
        fastq_1 = fastq_1[i],
        fastq_2 = fastq_2[i],
        readgroup_name = readgroup_name[i],
        library_name = library_name[i],
        platform_unit = platform_unit[i],
        run_date = run_date[i],
        platform_name = platform_name[i],
        sequencing_center = sequencing_center[i],
        target = primary_path + "/fastq/" + flowcell[i]
    }
  }

  # Create a file with a list of the generated ubams
  call CreateFoFN {
    input:
      array_of_files = PairedFastQsToUnmappedBAM.output_bam,
      fofn_name = ubam_list_name,
  }

  # Copy file containing list of generated ubams to appropriate directory
  call CopyFoFN {
    input:
      list = CreateFoFN.fofn_list,
      target = out_path
  }

###############################################################################
#
# Task calls for generating variant call-ready BAM from uBAMs
#
###############################################################################

  # Get the version of BWA to include in the PG record in the header of the BAM produced by MergeBamAlignment. 
  call GetBwaVersion {
    input: 
      bwa = bwa
  }

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in read_lines(CreateFoFN.fofn_list)) {

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)

    # Map reads to reference
    call SamToFastqAndBwaMem {
      input:
        bwa = bwa,
        picard = picard,
        samtools = samtools,
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = bam_basename + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        compression_level = compression_level
     }

    # Merge original uBAM and BWA-aligned BAM 
    call MergeBamAlignment {
      input:
        gatk = gatk,
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = bam_basename + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        compression_level = compression_level
    }
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.

  call MarkDuplicates {
    input:
      gatk = gatk,
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      compression_level = compression_level,
      target = primary_path + "/logs/"
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      gatk = gatk,
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      compression_level = compression_level
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        gatk = gatk,
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_snp_vcf = dbsnp_snp_vcf,
        dbsnp_snp_vcf_index = dbsnp_snp_vcf_index,
        broad_snp_vcf = broad_snp_vcf,
        broad_snp_vcf_index = broad_snp_vcf_index,
        axelsson_snp_vcf = axelsson_snp_vcf,
        axelsson_snp_vcf_index = axelsson_snp_vcf_index,
        dbsnp_indels_vcf = dbsnp_indels_vcf,
        dbsnp_indels_vcf_index = dbsnp_indels_vcf_index,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }  
  }  
  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      gatk = gatk,
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      target = primary_path + "/logs/"
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        gatk = gatk,
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      gatk = gatk,
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      compression_level = compression_level,
      target = primary_path + "/bam/"
  }

  call DepthOfCoverageAndFlagStat {
    input:
      gatk3 = gatk3,
      input_bam = GatherBamFiles.output_bam,
      doc_output_bam_basename = base_file_name + ".depthofcoverage",
      fs_output_bam_basename = base_file_name + ".flagstat",
      ref_fasta = ref_fasta,
      target = primary_path + "/bam/"
  }

###############################################################################
#
# Task calls for generating gVCF from BAM
#
###############################################################################

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        gatk = gatk,
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bam_index,
        interval_list = interval_file,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        make_gvcf = making_gvcf,
    }
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      gatk = gatk,
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      target = primary_path + "/gvcf/"
  }

###############################################################################
#
# Outputs that will be retained when execution is complete
#
###############################################################################

  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
    File unmapped_bam_list = CreateFoFN.fofn_list
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File analysis_ready_bam = GatherBamFiles.output_bam
    File analysis_ready_bam_index = GatherBamFiles.output_bam_index
    File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
    File output_vcf = MergeGVCFs.output_vcf
    File output_vcf_index = MergeGVCFs.output_vcf_index
  }
}




###############################################################################
#
# Task definitions for converting paired fastq files to uBAMs
#
###############################################################################

task RunFastQC {
    String fastqc
    String target
    String r1
    String r2

    command {
      dir=${target} && \
      mkdir -p "$dir" && \
      ${fastqc} --outdir "$dir" ${r1} ${r2}
    }
}

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  # Command parameters
  String sample_name
  String fastq_1
  String fastq_2
  String readgroup_name
  String library_name
  String platform_unit
  String run_date
  String platform_name
  String sequencing_center

  String target
  String gatk

  Int mem_gb
  Int cpu
  String walltime

  command {
    ${gatk} --java-options "-Xmx6000m" \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ${readgroup_name} \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME ${library_name} \
    --PLATFORM_UNIT ${platform_unit} \
    --RUN_DATE ${run_date} \
    --PLATFORM ${platform_name} \
    --SEQUENCING_CENTER ${sequencing_center}

    dir=${target}
    mkdir -p "$dir"
    cp -t "$dir" ${fastq_1} ${fastq_2}
  }
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}

task CreateFoFN {
  # Command parameters
  Array[String] array_of_files
  String fofn_name
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
}

task CopyFoFN {

  String list
  String target 

  command {
    dir=${target} && \
    mkdir -p "$dir" && \
    cp ${list} "$dir"
  }
}

###############################################################################
#
# Task definitions for generating variant call-ready BAM from uBAMs
#
###############################################################################

# Get version of BWA
task GetBwaVersion {

  String bwa  

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ${bwa} 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  output {
    String version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  String input_bam
  String bwa_commandline
  String output_bam_basename
  String ref_fasta
  String ref_fasta_index
  String ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  String ref_amb
  String ref_ann
  String ref_bwt
  String ref_pac
  String ref_sa

  Int compression_level

  String bwa
  String picard
  String samtools
  String java_opt

  Int mem_gb
  Int cpu
  String walltime

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

		java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${picard} \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true \
    | \
		${bwa}${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
    | \
		${samtools} view -1 - > ${output_bam_basename}.bam

  >>>
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  String unmapped_bam
  String bwa_commandline
  String bwa_version
  String aligned_bam
  String output_bam_basename
  String ref_fasta
  String ref_fasta_index
  String ref_dict

  Int compression_level

  String gatk
  String java_opt

  Int mem_gb
  Int cpu
  String walltime

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${output_bam_basename}.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
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
      --PROGRAM_GROUP_VERSION "${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true
  }
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  
  Int compression_level
  
  String gatk
  String java_opt
  String tmp_dir
  
  String target

  Int mem_gb
  Int cpu
  String walltime

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    mkdir -p ${tmp_dir} &&
    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt} -Xmx16g" \
      MarkDuplicates \
      --TMP_DIR ${tmp_dir} \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true && \
      dir=${target} && \
      mkdir -p "$dir" && \
      cp ${metrics_filename} ${target}
  }
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
 #runtime {
 #  #continueOnReturnCode: [0, 250]
 #   continueOnReturnCode: [0]
 #}
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  String input_bam
  String output_bam_basename
  String ref_dict
  String ref_fasta
  String ref_fasta_index
  
  Int compression_level

  String gatk
  String java_opt_sort
  String java_opt_fix

  command {
    set -o pipefail

    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_sort}" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_fix}" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  }
  runtime {
    continueOnReturnCode: [0, 250]
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict  

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  String input_bam
  String input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  String dbsnp_snp_vcf
  String dbsnp_snp_vcf_index
  String broad_snp_vcf
  String broad_snp_vcf_index
  String axelsson_snp_vcf
  String axelsson_snp_vcf_index
  String dbsnp_indels_vcf
  String dbsnp_indels_vcf_index
  
  String ref_dict
  String ref_fasta
  String ref_fasta_index
  
  String gatk
  String java_opt

  Int mem_gb
  Int cpu
  String walltime

  command { 
    ${gatk} --java-options "${java_opt}" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbsnp_snp_vcf} \
      --known-sites ${broad_snp_vcf} \
      --known-sites ${axelsson_snp_vcf} \
      --known-sites ${dbsnp_indels_vcf} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {
  Array[String] input_bqsr_reports
  String output_report_filename

  String gatk
  String java_opt

  String target

  command {
    ${gatk} --java-options "${java_opt}" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename} && \
      cp ${output_report_filename} ${target}
  }
  runtime {
   #continueOnReturnCode: [0, 250]
    continueOnReturnCode: [0]
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  String input_bam
  String input_bam_index
  String output_bam_basename
  String recalibration_report
  Array[String] sequence_group_interval
  String ref_dict
  String ref_fasta
  String ref_fasta_index

  String gatk
  String java_opt

  Int mem_gb
  Int cpu
  String walltime

  command {  
    ${gatk} --java-options "${java_opt}" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[String] input_bams
  String output_bam_basename

  Int compression_level
  
  String gatk
  String java_opt

  String target

  command {
    ${gatk} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
      GatherBamFiles \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true && \
      dir=${target} && \
      mkdir -p "$dir" && \
      cp -t ${target} ${output_bam_basename}.bam \
                      ${output_bam_basename}.bai \
                      ${output_bam_basename}.bam.md5
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Get depth of coverage from variant call ready bam
task DepthOfCoverageAndFlagStat {
    String input_bam
    String doc_output_bam_basename
    String fs_output_bam_basename
    String ref_fasta

    String gatk3
    String target

    command {
      java -jar -Xmx32g ${gatk3} -T DepthOfCoverage \
        -R ${ref_fasta} \
        -omitBaseOutput \
        -omitLocusTable \
        -omitIntervals \
        -I ${input_bam} \
        -o ${doc_output_bam_basename} \
        -ct 5 \
        -ct 15 \
        -ct 30 \
        -nt 8 && \
        cp -t ${target} ${doc_output_bam_basename}.sample_summary \
                ${doc_output_bam_basename}.sample_statistics && \
      java -jar -Xmx32g ${gatk3} -T FlagStat \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -o ${fs_output_bam_basename} \
        -nct 8 && \
        cp ${fs_output_bam_basename} ${target}
    }
}

###############################################################################
#
# Task definitions for generating gVCF from BAM
#
###############################################################################

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Boolean make_gvcf

  String gatk
  String? java_options
  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

  # Runtime parameters
  Int heap_max
#  String pbs
  Int mem_gb
  Int cpu
  String walltime

  command <<<
  set -e
  
    ${gatk} --java-options "-Xmx${heap_max}G ${java_opt}" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -L ${interval_list} \
      -O ${output_filename} \
      -contamination ${default=0 contamination} ${true="-ERC GVCF" false="" make_gvcf}
  >>>
  runtime {
    backend: "slurm"
    mem_gb: mem_gb
    cpu: cpu
    walltime: walltime
  }
  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_filename

  String gatk

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1

  String target

  command <<<
  set -e

    ${gatk} --java-options "-Xmx${command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ${sep=' --INPUT ' input_vcfs} \
      --OUTPUT ${output_filename} && \
      dir=${target} && \
      mkdir -p "$dir" && \
      cp -t ${target} ${output_filename} ${output_filename}.tbi
  >>>  
  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}
