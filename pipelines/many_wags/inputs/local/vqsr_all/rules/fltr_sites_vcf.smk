
rule fltr_make_sites_only:
    input:
        vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz" 
    output:
        var_filtrd_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz"
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
        "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz",
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
        tmp_vcf = lambda wildcards: "{}/wgs/pipeline/{}/{}/sites_only_gather_vcf/tmp.vcf.gz".format(
            wildcards.bucket, wildcards.ref, wildcards.date
        ),
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

