
rule fastqc:
    input:
        unpack(get_fastq)
    output:
        r1_zip  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R1_fastqc.zip"),
        r1_html = S3.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R1_fastqc.html"),
        r2_zip  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R2_fastqc.zip"),
        r2_html = S3.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R2_fastqc.html"),
    params:
        outdir = "{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}",
        r1_html = lambda wildcards, input: os.path.basename(input.r1.replace(".fastq.gz","_fastqc.html")),
        r1_zip  = lambda wildcards, input: os.path.basename(input.r1.replace(".fastq.gz","_fastqc.zip")),
        r2_html = lambda wildcards, input: os.path.basename(input.r2.replace(".fastq.gz","_fastqc.html")),
        r2_zip  = lambda wildcards, input: os.path.basename(input.r2.replace(".fastq.gz","_fastqc.zip")),
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}.fastqc.benchmark.txt"
    threads: 6
    resources:
         time   = 360,
         mem_mb = 6000, 
    shell:
        '''
            set -e

            fastqc -t {threads} \
                --outdir {params.outdir} \
                {input.r1} {input.r2}

            mv {params.outdir}/{params.r1_html} {output.r1_html}
            mv {params.outdir}/{params.r2_html} {output.r2_html}
            mv {params.outdir}/{params.r1_zip} {output.r1_zip}
            mv {params.outdir}/{params.r2_zip} {output.r2_zip}
        '''

rule flagstat:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
        flagstat = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.flagstat.txt"),
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.flagstat.benchmark.txt"
    threads: 8
    resources:
         time   = 60,
         mem_mb = 12000
    shell:
        '''
            samtools flagstat \
                {input.final_bam} \
                --thread 8 \
                > {output.flagstat}
        '''

rule bqsr_analyze_covariates:
    input:
        bqsr_before = "{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.recal_data.txt",
        bqsr_after  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.second_recal_data.txt",
    output:
        ac_pdf = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.analyze_cov.pdf"),
        ac_csv = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.analyze_cov.csv"),
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.bqsr_analyze_cov.benchmark.txt"
    resources:
         time   = 10,
         mem_mb = 4000
    shell:
        '''
            gatk AnalyzeCovariates \
                -before {input.bqsr_before} \
                -after {input.bqsr_after} \
                -plots {output.ac_pdf} \
                -csv {output.ac_csv}
        '''

rule qualimap_bamqc:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
        report_pdf  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/report.pdf",
        report_html = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/qualimapReport.html",
        report_txt  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/genome_results.txt"),
        raw_dat     = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/raw_data_qualimapReport/"),
        figures     = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/images_qualimapReport/"),
        css         = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/css/"),
    params:
        out_dir = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/",
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/{sample_name}.multimap_bamqc.benchmark.txt"
    threads: 12
    resources:
         time   = 120,
         mem_mb = 24000
    shell:
        '''
            unset DISPLAY

            qualimap bamqc \
                -bam {input.final_bam} \
                -outdir {params.out_dir} \
                -outformat PDF:HTML \
                -nt {threads} \
                --java-mem-size=24G
        '''

# NOTE THERE IS AN ISSUE WTIH CALLING IT DBSNP IF A USER DOES NOT HAVE THAT RESOURCE...
rule collect_metrics_on_vcf:
    input:
        final_vcf       = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz"),
        final_vcf_index = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz.tbi"),
        acgt_ivals      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.N50.interval_list",
    output:
        detail_metrics  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.variant_calling_detail_metrics"),
        summary_metrics = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.variant_calling_summary_metrics"),
    params:
        dbsnp_snp_vcf  = config['known_sites']['dbsnp_snp_vcf'],
        ref_dict       = config['ref_dict'],
        metrics_prefix = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}"
    threads: 8
    resources:
         time   = 360,
         mem_mb = 22000
    shell:
        '''
            gatk --java-options "-Xmx18g -Xms6g" \
                CollectVariantCallingMetrics \
                --INPUT {input.final_vcf} \
                --DBSNP {params.dbsnp_snp_vcf} \
                --SEQUENCE_DICTIONARY {params.ref_dict} \
                --OUTPUT {params.metrics_prefix} \
                --THREAD_COUNT 8 \
                --TARGET_INTERVALS {input.acgt_ivals}
        '''

rule bcftools_stats:
    input:
        final_vcf       = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz"),
        final_vcf_index = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz.tbi")
    output:
        all_stats = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}.{ref}.vchk"),
    params:
        ref_fasta = config['ref_fasta'],
    threads: 1
    resources:
         time   = 720,
         mem_mb = 6000
    shell:
        '''
            bcftools stats \
                -F {params.ref_fasta} \
                -s - {input.final_vcf} \
                > {output.all_stats}
        '''

rule bcftools_plot:
    input:
        all_stats = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}.{ref}.vchk"),
    output:
        summary = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}/summary.pdf")
    params:
        prefix = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}",
    threads: 1
    resources:
         time   = 60,
         mem_mb = 6000
    shell:
        '''
            plot-vcfstats \
                -p {params.prefix} \
                {input.all_stats}
        '''

#rule qc_cohort:
#    input:
#        html_veps = sorted(
#            expand(
#                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/recal.{interval}.vep.vcf_summary.html", 
#                bucket=config['bucket'],
#                ref=config['ref'],
#                date=config['date'],
#                interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
#            )
#        ),
#        len_barplt      = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths_mqc.tiff"),
#        all_stats       = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.{ref}.vchk"),
#        summary         = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort/summary.pdf"),
#        detail_metrics  = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_detail_metrics"),
#        summary_metrics = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_summary_metrics"),
#    output: 
#        S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html")
#    params:
#        outdir = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/"
#    threads: 4
#    resources:
#        time   = 360,
#        mem_mb = 12000
#    shell:
#        '''
#            multiqc {wildcards.bucket} \
#                --interactive \
#                --force \
#                -o {params.outdir}
#        '''

rule multiqc:
    input:
        r1_zip          = S3.remote(expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/fastqc/{{sample_name}}/{u.flowcell}/{u.readgroup_name}_R1_fastqc.zip",
            u=units.itertuples()
        )),
        r2_zip          = S3.remote(expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/fastqc/{{sample_name}}/{u.flowcell}/{u.readgroup_name}_R2_fastqc.zip",
            u=units.itertuples()
        )),
        flagstat        = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.flagstat.txt"),
        bqsr_before     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.recal_data.txt",
        bqsr_after      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/logs/{sample_name}.{ref}.second_recal_data.txt",
        ac_csv          = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.analyze_cov.csv"),
        metrics         = expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/{{ref}}/bam/{u.readgroup_name}.mark_adapt.metrics.txt",
            u=units.itertuples()
        ),
        report_txt      = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/genome_results.txt"),
        raw_dat         = rules.qualimap_bamqc.output.raw_dat,
        all_stats       = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}.{ref}.vchk"),
        summary         = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/stats/{breed}_{sample_name}/summary.pdf"),
        detail_metrics  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.variant_calling_detail_metrics"),
        summary_metrics = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.variant_calling_summary_metrics"),
    output: 
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc.log"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_picard_dups.txt"),
       #S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_picard_mark_illumina_adapters.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_samtools_flagstat.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_fastqc.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_general_stats.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_sources.txt"),
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_data.json"),
    params:
        outdir = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc"
    shell:
        '''
            multiqc {wildcards.bucket} \
                -x '*duplicates_marked*' -x '*group_*' -e snippy \
                --interactive \
                --force \
                -o {params.outdir}
        '''

