
rule vars_to_table:
    input:
        ival_vep      = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz", 
        ival_vep_tbi  = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz.tbi", 
    output:
        ival_vars = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.table"
    params:
        ref_fasta = config['ref_fasta'],
    threads: 12
    resources:
         time   = 840,
         mem_mb = 12000
    shell:
        '''
            gatk \
                VariantsToTable \
                -R {params.ref_fasta} \
                -V {input.ival_vep} \
                -F CHROM -F POS -F REF -F ALT -F FILTER -F AF -F HOM-REF -F HET -F HOM-VAR -F NO-CALL -F CSQ \
                --show-filtered \
                -O {output.ival_vars}
        '''

def get_var_tables(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{vep_interval}-scattered.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.table",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        vep_interval=INTERVALS
    ))

rule all_var_to_table:
    input:
        get_var_tables
    output:
        all_table = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{ref}.{date}.table")
    threads: 4
    resources:
         time   = 120,
         mem_mb = 12000
    shell:
        '''
            awk 'FNR == 1 && NR != 1 {{ next; }} {{ print; }}' {input} > {output.all_table}
        '''


#rule all_var_to_table:
#    input:
#        vep_vcf = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz"),
#        vep_tbi = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi"),
#    output:
#        all_vars_table = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{ref}.{date}.table")
#    params:
#        ref_fasta = config['ref_fasta'],
#    threads: 12
#    resources:
#         time   = 840,
#         mem_mb = 12000
#    shell:
#        '''
#            gatk \
#                VariantsToTable \
#                -R {params.ref_fasta} \
#                -V {input.vep_vcf} \
#                -F CHROM -F POS -F REF -F ALT -F FILTER -F AF -F HOM-REF -F HET -F HOM-VAR -F NO-CALL -F CSQ \
#                --show-filtered \
#                -O {output.all_vars_table}
#        '''
