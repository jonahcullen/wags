
rule plot_interval_lengths:
    input:
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz"
    output:
        len_barplt = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths_mqc.tiff"
    params:
        lengths = "{}/wgs/pipeline/{}/{}/intervals/collapsed_lengths.csv".format(
            config['bucket'], config['ref'], config['date']
        )
    threads: 1
    resources:
         time   = 30,
         mem_mb = 6000
    script:
        '../src/interval_plot.R'

rule bcftools_stats:
    input:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz",
    output:
        all_stats = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/join_call.{ref}.{date}.vchk",
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
        all_stats = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/join_call.{ref}.{date}.vchk",
    output:
        summary = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort/summary.pdf"
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

def get_vep_htmls(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{vep_interval}-scattered.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf_summary.html", 
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        vep_interval=INTERVALS
    ))

rule qc_cohort:
    input:
        get_vep_htmls,
        len_barplt      = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/interval_lengths_mqc.tiff",
        all_stats       = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/join_call.{ref}.{date}.vchk",
        summary         = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/{ref}_{date}_cohort/summary.pdf",
    output: 
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/multiqc_report.html"
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

