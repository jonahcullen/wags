
def get_unfiltered_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{interval}-scattered.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/genotype_gvcfs/money_{interval}/output.vcf.gz",
       #"{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz",
        bucket=config['bucket'],
        breed=breed,
        sample_name=sample_name,
        ref=config['ref'],
        interval=INTERVALS
    ))

rule gather_unfltr_vcf:
    input:
        get_unfiltered_vcfs
    output:
        var_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz",
        var_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz.tbi"
    params:
        vcfs    = lambda wildcards, input: " --input ".join(map(str,input)),
       #tmp_vcf = lambda wildcards: f"{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/money/unfltr_vcf/tmp.vcf.gz",
        tmp_vcf = lambda wildcards: "{}/wgs/{}/{}/{}/money/unfltr_vcf/tmp.vcf.gz".format(
            wildcards.bucket, wildcards.breed, wildcards.sample_name, wildcards.ref
        ),
        tmp_dir = config['tmp_dir']['unfilt_gather_vcf']
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
        var_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz",
        var_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz.tbi"
    output:
        nonsnp_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/nonsnps.vcf.gz",
        nonsnp_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/nonsnps.vcf.gz.tbi"
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

rule select_snps:
    input:
        var_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz",
        var_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/all_vars.vcf.gz.tbi"
    output:
        snp_unfiltered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/snp_only.vcf.gz",
        snp_unfiltered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/unfltr_vcf/snp_only.vcf.gz.tbi"
    threads: 4
    resources:
         time   = 240,
         mem_mb = 14000
    shell:
        '''
            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.var_unfiltered_vcf} \
                -O {output.snp_unfiltered_vcf} \
                -select-type SNP
        '''

