
def get_unfiltered_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.generate_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"wags_{interval}.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        interval=INTERVALS
    ))

rule gather_unfltr_vcf:
    input:
        get_unfiltered_vcfs
    output:
        var_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/all_vars.vcf.gz",
        var_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/all_vars.vcf.gz.tbi"
    params:
        vcfs    = lambda wildcards, input: " --input ".join(map(str,input)),
        tmp_vcf = lambda wildcards: f"{wildcards.bucket}/wgs/pipeline/{wildcards.ref}/{wildcards.date}/unfltr_vcf/tmp.vcf.gz",
        tmp_dir = config['tmp_dir']['unfilt_gather_vcf']
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 20000
    shell:
        '''
            set -e

            gatk --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {params.vcfs} \
                --output {params.tmp_vcf}
           
            java -jar -Xmx190g -Xms6g /opt/wags/src/picard.jar \
                SortVcf \
                TMP_DIR={params.tmp_dir} \
                I={params.tmp_vcf} \
                O={output.var_unfiltered_vcf}
            
            rm -f {params.tmp_vcf}
            rm -rf {params.tmp_dir}
        '''

rule select_nonsnps:
    input:
        var_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/all_vars.vcf.gz",
        var_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/all_vars.vcf.gz.tbi"
    output:
        nonsnp_unfiltered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz",
        nonsnp_unfiltered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/unfltr_vcf/nonsnps.vcf.gz.tbi"
    threads: 4
    resources:
         time   = 240,
         mem_mb = 14000
    shell:
        '''
            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.var_unfiltered_vcf} \
                -O {output.nonsnp_unfiltered_vcf} \
                -xl-select-type SNP
        '''

rule fltr_make_sites_only:
    input:
        vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz" 
    output:
        var_filtrd_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.variant_filtered.vcf.gz",
        sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.sites_only.variant_filtered.vcf.gz"
    params:
        excess_het = "54.69",
        ref_fasta  = config["ref_fasta"]
    resources:
         time   = 30,
         mem_mb = 12000
    shell:
        '''
            gatk --java-options "-Xmx3g -Xms3g" \
            VariantFiltration \
                --filter-expression "ExcessHet > {params.excess_het}" \
                --filter-name ExcessHet \
                -O {output.var_filtrd_vcf} \
                -V {input.vcf}

            gatk --java-options "-Xmx3g -Xms3g" \
            MakeSitesOnlyVcf \
                --INPUT {output.var_filtrd_vcf} \
                --OUTPUT {output.sites_only_vcf}
        '''

def get_filtered_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.generate_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"wags_{interval}.interval_list"))
    # return list of recal vcfs
    return expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.sites_only.variant_filtered.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        interval=INTERVALS
    )

rule sites_only_gather_vcf:
    input:
        get_filtered_vcfs
    output:
        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    params:
        vcfs    = lambda wildcards, input: " --input ".join(map(str,input)),
        tmp_vcf = lambda wildcards: f"{wildcards.bucket}/wgs/pipeline/{wildcards.ref}/{wildcards.date}/sites_only_gather_vcf/tmp.vcf.gz",
        tmp_dir = config['tmp_dir']['sites_only_gather_vcf']
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 20000
    shell:
        '''
            set -e

            gatk --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {params.vcfs} \
                --output {params.tmp_vcf}
           
            java -jar /opt/wags/src/picard.jar \
                SortVcf \
                TMP_DIR={params.tmp_dir} \
                I={params.tmp_vcf} \
                O={output.gather_sites_only_vcf}
            
            rm -f {params.tmp_vcf}
            rm -rf {params.tmp_dir}
        '''

rule select_snps:
    input:
        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    output:
        snp_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz",
        snp_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.snp_sites_only.vcf.gz.tbi"
    threads: 4
    resources:
         time   = 240,
         mem_mb = 14000
    shell:
        '''
            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.gather_sites_only_vcf} \
                -O {output.snp_sites_only_vcf} \
                -select-type SNP
        '''

rule fltr_snps_only_var_filtrd:
    input:
        var_filtrd_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.variant_filtered.vcf.gz",
    output:
        snps_var_filtrd_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/fltr.{interval}.snps_only.variant_filtered.vcf.gz",
    params:
        excess_het = "54.69",
        ref_fasta  = config["ref_fasta"]
    resources:
         time   = 30,
         mem_mb = 12000
    shell:
        '''
            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.var_filtrd_vcf} \
                -O {output.snps_var_filtrd_vcf} \
                -select-type SNP
        '''

