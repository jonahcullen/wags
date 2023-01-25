
rule all_var_to_table:
    input:
        vep_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
        vep_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi",
    output:
        all_vars_table = "{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{ref}.{date}.table"
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
                -V {input.vep_vcf} \
                -F CHROM -F POS -F REF -F ALT -F FILTER -F AF -F HOM-REF -F HET -F HOM-VAR -F NO-CALL -F CSQ \
                --show-filtered \
                -O {output.all_vars_table}
        '''
