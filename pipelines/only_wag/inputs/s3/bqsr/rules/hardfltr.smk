
rule nonsnps_hard_fltr:
    input:
       #nonsnp_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz",
       #nonsnp_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz.tbi"
        nonsnp_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/nonsnps.vcf.gz",
        nonsnp_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/nonsnps.vcf.gz.tbi"
    output:
        nonsnp_filtered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/nonsnp_fltr.vcf.gz",
        nonsnp_filtered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/nonsnp_fltr.vcf.gz.tbi"
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx24g -Xms3g" \
                VariantFiltration \
                -V {input.nonsnp_unfiltered_vcf} \
                -O {output.nonsnp_filtered_vcf} \
                --filter-name "QD2" --filter-expression "QD < 2.0" \
                --filter-name "FS200" --filter-expression "FS > 200.0" \
                --filter-name "ReadPosRankSum-2" --filter-expression "ReadPosRankSum < -2.0" \
                --filter-name "SOR-10" --filter-expression "SOR > 10.0" \
                --verbosity ERROR
        '''

rule snps_hard_fltr:
    input:
       #snp_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz",
       #snp_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz.tbi"
        snp_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz",
        snp_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz.tbi"
    output:
       #snp_filtered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/snp_fltr.vcf.gz",
       #snp_filtered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/snp_fltr.vcf.gz.tbi"
        snp_filtered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/snp_fltr.vcf.gz",
        snp_filtered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/snp_fltr.vcf.gz.tbi"
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx24g -Xms3g" \
                VariantFiltration \
                -V {input.snp_unfiltered_vcf} \
                -O {output.snp_filtered_vcf} \
                --filter-name "QD2" -filter "QD < 2.0" \
                --filter-name "QUAL30" -filter "QUAL < 30.0" \
                --filter-name "SOR3" -filter "SOR > 3.0" \
                --filter-name "FS60" -filter "FS > 60.0" \
                --filter-name "MQ40" -filter "MQ < 40.0" \
                --filter-name "MQRankSum-12.5" -filter "MQRankSum < -12.5" \
                --filter-name "ReadPosRankSum-8" -filter "ReadPosRankSum < -8.0" \
                --verbosity ERROR
        '''

