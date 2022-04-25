
rule filter_x_chrom:
    input:
        final_vcf       = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz"),
        final_vcf_index = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.vcf.gz.tbi")
    output:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/joint_genotype.{ref}.{chrom}.vcf.gz",
    params:
       #gatk = config["gatk"],
       #vcfs = lambda wildcards, input: " --input ".join(map(str,input.recal_vcfs)),
    threads: 4
    resources:
         time   = 60,
         mem_mb = 16000
    shell:
        '''
            # change chr39 to X to extract
            if [ ${wildcards.chrom} = "chr39" ]; then
                chrom=chrX
            else
                chrom=${wildcards.chrom}
            fi

            bcftools view \
                --min-alleles 2 \
                --max-alleles 2 \
                --types snps \
                --regions $chrom \
                --exclude ' GT="." '
                {input.final_vcf} \
            | \
            bcftools filter \
                --exclude "F_MISSING > 0.01" \
                -Oz \
                -o {output.chrom_vcf}
        '''

rule filter_x_chrom:
    input:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/joint_genotype.{ref}.{chrom}.vcf.gz",
    output:
        phase_vcf       = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_genotype.{ref}.{chrom}.phased.vcf.gz",
        phase_vcf_index = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_genotype.{ref}.{chrom}.phased.vcf.gz.tbi",
    params:
        eff_pop_size = '200',
        window       = '120',
        overlap      = '10'
    threads: 24
    resources:
         time   = 60,
         mem_mb = 248000
    shell:
        '''
            java -jar -Xmx246g /opt/wags/src/beagle.18May20.d20.jar \
                gt={input.chrom_vcf} \
                ne={params.eff_pop_size} \
                nthreads={threads} \
                map=canFam3.linkage.map.wgs \
                window={params.window} \
                overlap={params.overlap} \
                out=joint_genotype.canfam3.snps.chrxxx.phased

            tabix -p vcf joint_genotype.canfam3.snps.chrxxx.phased.vcf.gz
        '''
