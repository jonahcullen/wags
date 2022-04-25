
rule all_var_to_table:
    input:
        vep_vcf       = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vep.vcf.gz",keep_local=True),
        vep_vcf_index = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vep.vcf.gz.tbi",keep_local=True),
    output:
        all_vars_table = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/var_to_table/all.{date}.{ref}.table")
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
