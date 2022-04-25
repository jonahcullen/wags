
#def get_vcfs(wildcards):
#    interval_dir = checkpoints.generate_intervals.get(**wildcards).output.ivals
#    intervals = glob_wildcards(f"{interval_dir}/{{interval}}.interval_list").interval
#   #vcfs = expand(rules.fltr_make_sites_only.output.var_filtrd_vcf, **wildcards, interval = intervals)
#    vcfs = expand(
#        "{bucket}/wgs/friedlab/pipeline/{ref}/{date}/genotype_gvcfs/{interval}/filtr.{interval}.variant_filtered.vcf.gz", 
#        bucket = config['bucket'],
#        ref = config['ref'],
#        date = config['date'],
#        interval = intervals
#    )
#    return vcfs

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

rule sites_only_gather_vcf:
    input:
        sites_only_vcf = expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.variant_filtered.vcf.gz",
            bucket = config['bucket'],
            ref = config['ref'],
            date = config['date'],
            interval = [str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
        )
    output:
        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    params:
        vcfs    = lambda wildcards, input: " --input ".join(map(str,input.sites_only_vcf)),
        tmp_vcf = lambda wildcards: f"{wildcards.bucket}/wgs/pipeline/{wildcards.ref}/{wildcards.date}/sites_only_gather_vcf/tmp.vcf.gz",
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

#rule sites_only_gather_vcf:
#    input:
#        sites_only_vcf = expand(
#            "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/filtr.{interval}.variant_filtered.vcf.gz",
#            bucket = config['bucket'],
#            ref = config['ref'],
#            date = config['date'],
#            interval = [str(i).zfill(4) for i in range(0,config['num_intervals']+1)]
#        )
#    output:
#        gather_sites_only_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz",
#        gather_sites_only_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
#    params:
#        vcfs    = lambda wildcards, input: " --input ".join(map(str,input.sites_only_vcf)),
#        tmp_vcf = lambda wildcards: f"{wildcards.bucket}/wgs/pipeline/{wildcards.ref}/{wildcards.date}/sites_only_gather_vcf/tmp.vcf.gz",
#        tmp_dir = config['tmp_dir']['sites_only_gather_vcf']
#    threads: 4
#   #resources:
#   #     time   = 240,
#   #     mem_mb = 24000
#    shell:
#        '''
#            set -e
#
#            mkdir -p {params.tmp_dir}
#
#            gatk --java-options "-Xmx18g -Xms6g" \
#                GatherVcfsCloud \
#                --tmp-dir {params.tmp_dir} \
#                --ignore-safety-checks \
#                --gather-type BLOCK \
#                --input {params.vcfs} \
#                --output {params.tmp_vcf}
#            
#            java -jar /opt/wags/src/picard.jar \
#                SortVcf \
#                TMP_DIR={params.tmp_dir} \
#                I={params.tmp_vcf} \
#                O={output.gather_sites_only_vcf}
#            
#            gatk --java-options "-Xmx18g -Xms6g" \
#                IndexFeatureFile \
#                --input {output.gather_sites_only_vcf}
#
#            rm -f {params.tmp_vcf}
#        '''
