
rule indels_var_recal:
    input:
        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    output:
        indels_recal    = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.recal",
        indels_tranches = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.tranches",
        indels_plot     = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.plot.R",
        indels_plot_pdf = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.plot.R.pdf"
    params:
        tranche_values    = " -tranche ".join(map(str,config["indel_recalibration_tranche_values"])),
        annotation_values = " -an ".join(map(str,config["indel_recalibration_annotation_values"])),
        dbsnp_indels_vcf  = config["dbsnp_indels_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx24g -Xms24g" \
                VariantRecalibrator \
                -V {input.gather_sites_only_vcf} \
                --trust-all-polymorphic \
                -tranche {params.tranche_values} \
                -an {params.annotation_values} \
                -mode INDEL \
                --max-gaussians 4 \
                --resource:dbsnp,known=false,training=true,truth=true,prior=10 {params.dbsnp_indels_vcf} \
                -O {output.indels_recal} \
                --tranches-file {output.indels_tranches} \
                --rscript-file {output.indels_plot}
        '''

rule snps_var_recal:
    input:
        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    output:
        snps_recal    = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.recal",
        snps_tranches = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches",
        snps_plot     = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.plot.R",
        snps_plot_pdf = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.plot.R.pdf",
        tranches_pdf  = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches.pdf",
    params:
        tranche_values    = " -tranche ".join(map(str,config["snp_recalibration_tranche_values"])),
        annotation_values = " -an ".join(map(str,config["snp_recalibration_annotation_values"])),
        dbsnp_snp_vcf     = config["dbsnp_snp_vcf"],
        broad_snp_vcf     = config["broad_snp_vcf"],
        axelsson_snp_vcf  = config["axelsson_snp_vcf"],
        illumina_snp_vcf  = config["illumina_snp_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx24g -Xms3g" \
                VariantRecalibrator \
                -V {input.gather_sites_only_vcf} \
                --trust-all-polymorphic \
                -tranche {params.tranche_values} \
                -an {params.annotation_values} \
                -mode SNP \
                --max-gaussians 6 \
                --resource:illumina,known=false,training=true,truth=true,prior=15 {params.illumina_snp_vcf} \
                --resource:axelsson,known=false,training=true,truth=true,prior=12 {params.axelsson_snp_vcf} \
                --resource:broad,known=false,training=true,truth=false,prior=10 {params.broad_snp_vcf} \
                --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp_snp_vcf} \
                -O {output.snps_recal} \
                --tranches-file {output.snps_tranches} \
                --rscript-file {output.snps_plot}
        '''

rule apply_recal:
    input:
        var_filtrd_vcf  = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        indels_recal    = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.recal",
        indels_tranches = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/indels.tranches",
        snps_recal      = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.recal",
        snps_tranches   = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/snps.tranches",
    output:
        recal_vcf       = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal.{interval}.vcf.gz",
        recal_vcf_index = "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal.{interval}.vcf.gz.tbi"
    params:
        apply_recal_dir    = lambda wildcards: f"{wildcards.bucket}/wgs/pipeline/{wildcards.ref}/{wildcards.date}/var_recal/apply/wags_{wildcards.interval}",
        indel_filter_level = config["indel_filter_level"],
        snp_filter_level   = config["snp_filter_level"]
    resources:
         time   = 30,
         mem_mb = 16000
    shell:
        '''
            set -e

            mkdir -p {params.apply_recal_dir}

            gatk --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -V {input.var_filtrd_vcf} \
                --recal-file {input.indels_recal} \
                --tranches-file {input.indels_tranches} \
                --truth-sensitivity-filter-level {params.indel_filter_level} \
                --create-output-variant-index true \
                -mode INDEL \
                -O {params.apply_recal_dir}/tmp.indel.recalibrated.vcf \

            gatk --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -V {params.apply_recal_dir}/tmp.indel.recalibrated.vcf \
                --recal-file {input.snps_recal} \
                --tranches-file {input.snps_tranches} \
                --truth-sensitivity-filter-level {params.snp_filter_level} \
                --create-output-variant-index true \
                -mode SNP \
                -O {output.recal_vcf} \

            rm -f {params.apply_recal_dir}/tmp.indel.recalibrated.vcf
        '''
