
rule plot_interval_lengths:
    input:
        lengths = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/collapsed_lengths.csv"
    output:
        len_barplt = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths_mqc.tiff")
    threads: 1
    resources:
         time   = 30,
         mem_mb = 6000
    script:
        '../src/interval_plot.R'

rule collect_metrics_on_vcf:
    input:
        final_vcf       = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz"),
        final_vcf_index = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz.tbi"),
        eval_ival_list  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list"
    output:
        detail_metrics  = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_detail_metrics"),
        summary_metrics = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_summary_metrics"),
    params:
        dbsnp_snp_vcf  = config['dbsnp_snp_vcf'],
        ref_dict       = config['ref_dict'],
        metrics_prefix = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort"
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
                --TARGET_INTERVALS {input.eval_ival_list}
        '''

rule bcftools_stats:
    input:
        final_vcf = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz"),
    output:
        all_stats = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.{ref}.vchk"),
    params:
        ref_fasta = config['ref_fasta'],
        conda_env = config['conda_envs']['qc']
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
        all_stats = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.{ref}.vchk"),
    output:
        summary = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort/summary.pdf")
    params:
        prefix    = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort",
        conda_env = config['conda_envs']['qc']
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

rule qc_cohort:
    input:
        html_veps = sorted(
            expand(
                "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/recal.{interval}.vep.vcf_summary.html", 
                bucket=config['bucket'],
                ref=config['ref'],
                date=config['date'],
                interval=[str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
            )
        ),
        len_barplt      = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths_mqc.tiff"),
        all_stats       = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.{ref}.vchk"),
        summary         = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort/summary.pdf"),
        detail_metrics  = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_detail_metrics"),
        summary_metrics = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort.variant_calling_summary_metrics"),
    output: 
        S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html")
    params:
        outdir = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/"
    threads: 4
    resources:
        time   = 360,
        mem_mb = 12000
    shell:
        '''
            multiqc {wildcards.bucket} \
                --interactive \
                --force \
                -o {params.outdir}
        '''

