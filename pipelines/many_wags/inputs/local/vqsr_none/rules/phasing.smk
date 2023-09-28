
rule filter_x_chrom:
    input:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.vcf.gz",
        final_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.vcf.gz.tbi"
    output:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/joint_call.{ref}.{chrom}.vcf.gz",
    threads: 4
    resources:
         time   = 60,
         mem_mb = 16000
    shell:
        '''
            # change chr39 to X to extract
            chrom={wildcards.chrom}
            if [ $chrom = "chr39" ]; then
                chrom=chrX
            fi

            bcftools view \
                --min-alleles 2 \
                --max-alleles 2 \
                --types snps \
                --regions $chrom \
                --exclude ' GT="." ' \
                {input.final_vcf} \
            | \
            bcftools filter \
                --exclude "F_MISSING > 0.01" \
                -Oz \
                -o {output.chrom_vcf}
        '''

rule phase_x_chrom:
    input:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/joint_call.{ref}.{chrom}.vcf.gz",
    output:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_call.{ref}.{chrom}.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_call.{ref}.{chrom}.phased.vcf.gz.tbi",
    params:
        link_map     = "/home/refgen/dog/canfam3/canFam3.linkage.map.wgs",
        eff_pop_size = 200,
        window       = 120,
        overlap      = 10,
        out_prefix   = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_call.{ref}.{chrom}.phased",
    threads: 24
    resources:
        partition = config['partitions']['phasing'],
        time      = 480,
        mem_mb    = 248000
    shell:
        '''
            java -jar -Xmx246g /opt/wags/src/beagle.18May20.d20.jar \
                gt={input.chrom_vcf} \
                ne={params.eff_pop_size} \
                nthreads={threads} \
                map={params.link_map} \
                window={params.window} \
                overlap={params.overlap} \
                out={params.out_prefix}

            tabix -p vcf {output.phase_vcf}
        '''

rule concat_phased:
    input:
        phase_vcf = expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/joint_call.{ref}.{chrom}.phased.vcf.gz",
            bucket=config['bucket'],
            ref=config['ref'],
            date=config['date'],
            chrom=[f"chr{i}" for i in range(1,39+1)]
        ),
    output:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 720,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.phase_vcf} \
                {input.phase_vcf}

            tabix -p vcf {output.phase_vcf}
        '''

rule genomebrowse_index:
    input:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz.tbi",
    output:
        "{bucket}/wgs/pipeline/{ref}/{date}/phasing/joint_call.{ref}.{date}.snps.phased.vcf.gz.covtsf",
    params:
        ref_dir = os.path.dirname(config['ref_fasta']),
    threads: 4
    resources:
        time   = 2880,
        mem_mb = 12000
    shell:
        '''
            gautil coverage \
                {input.phase_vcf} \
                --refFolder={params.ref_dir}
        '''
