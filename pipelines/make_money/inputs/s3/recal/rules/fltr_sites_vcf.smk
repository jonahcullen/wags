
rule fltr_make_sites_only:
    input:
        vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_{split}/output.vcf.gz",
    output:
        var_filtrd_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_00{split}/filtr.00{split}.variant_filtered.vcf.gz",
        sites_only_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_00{split}/filtr.00{split}.sites_only.variant_filtered.vcf.gz"
    params:
        excess_het = config['excess_het_threshold'],
        ref_fasta  = config['ref_fasta']
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

def get_sites_only_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals up to scatter_size set in config (default: 50)
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"00{split}-scattered.interval_list"))
    # return list of split intervals var_filtered.vcf.gz
    return expand(
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_00{split}/filtr.00{split}.sites_only.variant_filtered.vcf.gz",
        bucket = config['bucket'],
        breed=breed,
        sample_name = sample_name,
        ref = config['ref'],
        split = SPLIT
    )

rule sites_only_gather_vcf:
    input:
        get_sites_only_vcfs
    output:
        gather_sites_only_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    params:
        vcfs    = lambda wildcards, input: " --input ".join(map(str,input)),
        tmp_vcf = lambda wildcards: f"{config['bucket']}/wgs/{breed}/{sample_name}/{config['ref']}/money/sites_only_gather_vcf/tmp.vcf.gz",
        tmp_dir = config['tmp_dir']['sites_only_gather_vcf']
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 200000
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

