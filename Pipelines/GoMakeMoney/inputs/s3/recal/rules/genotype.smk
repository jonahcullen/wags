
rule genotype_gvcfs:
    input:
        final_gvcf = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz"),
        interval   = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/00{split}-scattered.interval_list"
    output:
        vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_{split}/output.vcf.gz",
    params:
        ref_fasta = config['ref_fasta']
    threads: 6
    resources:
         time   = 60,
         mem_mb = 60000
    shell:
        '''
            gatk --java-options "-Xmx50g -Xms50g" \
                GenotypeGVCFs \
                -R {params.ref_fasta} \
                -O {output.vcf} \
                -G StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V {input.final_gvcf} \
                -L {input.interval}
        '''
