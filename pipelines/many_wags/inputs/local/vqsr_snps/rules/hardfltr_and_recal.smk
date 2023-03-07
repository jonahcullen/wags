
rule nonsnps_hard_fltr:
    input:
        nonsnp_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz",
        nonsnp_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz.tbi"
    output:
        nonsnp_filtered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz",
        nonsnp_filtered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz.tbi"
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


rule snps_var_recal:
    input:
        snp_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz",
        snp_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz.tbi"
    output:
        snps_recal    = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.recal",
        snps_tranches = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches",
        snps_plot     = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.plot.R",
        snps_plot_pdf = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.plot.R.pdf",
        tranches_pdf  = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches.pdf",
    params:
        tranche_values    = " -tranche ".join(map(str,config["snp_recalibration_tranche_values"])),
        annotation_values = " -an ".join(map(str,config["snp_recalibration_annotation_values"])),
        k9hd_axiom_vcf    = config["k9hd_axiom_vcf"],
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx24g -Xms3g" \
                VariantRecalibrator \
                -V {input.snp_sites_only_vcf} \
                --trust-all-polymorphic \
                -tranche {params.tranche_values} \
                -an {params.annotation_values} \
                -mode SNP \
                --max-gaussians 6 \
                --resource:array,known=false,training=true,truth=true,prior=12 {params.k9hd_axiom_vcf} \
                -O {output.snps_recal} \
                --tranches-file {output.snps_tranches} \
                --rscript-file {output.snps_plot}
        '''

rule snps_apply_recal:
    input:
        snps_var_filtrd_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.snps_only.variant_filtered.vcf.gz",
        snps_recal          = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.recal",
        snps_tranches       = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches",
    output:
        recal_snps_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal_snps.{interval}.vcf.gz",
        recal_vcf_tbi  = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal_snps.{interval}.vcf.gz.tbi"
    params:
        snp_filter_level   = config["snp_filter_level"]
    resources:
         time   = 30,
         mem_mb = 16000
    shell:
        '''
            gatk --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -V {input.snps_var_filtrd_vcf} \
                --recal-file {input.snps_recal} \
                --tranches-file {input.snps_tranches} \
                --truth-sensitivity-filter-level {params.snp_filter_level} \
                --create-output-variant-index true \
                -mode SNP \
                -O {output.recal_snps_vcf} \
        '''

