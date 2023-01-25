
def get_recal_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.generate_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"wags_{interval}.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal_snps.{interval}.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        interval=INTERVALS
    ))

rule gather_snp_recal_vcfs:
    input:
        get_recal_vcfs
    output:
        final_snp_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/snps.{ref}.vcf.gz",
        final_snp_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/snps.{ref}.vcf.gz.tbi"
    params:
        vcfs = lambda wildcards, input: " --input ".join(map(str,input)),
    threads: 4
    resources:
         time   = 240,
         mem_mb = 22000
    shell:
        '''
            set -e

            gatk --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {params.vcfs} \
                --output {output.final_snp_vcf}

            gatk --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {output.final_snp_vcf}
        '''

rule combine_snps_nonsnps:
    input:
        final_snp_vcf       = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/snps.{ref}.vcf.gz",
        final_snp_tbi       = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/snps.{ref}.vcf.gz.tbi",
        nonsnp_filtered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz",
        nonsnp_filtered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz.tbi"
    output:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz",
        final_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi",
    threads: 12
    resources:
         time   = 600,
         mem_mb = 24000
    shell:
        '''
            bcftools concat \
                --threads 12 \
                --allow-overlaps \
                {input.final_snp_vcf} {input.nonsnp_filtered_vcf} \
                -Oz \
                -o {output.final_vcf}

            tabix -p vcf {output.final_vcf}
        '''

rule vep_by_interval:
    input:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz",
        final_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi",
        interval  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/import/wags_{interval}.interval_list",
    output:
        final_interval    = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/wags_{interval}/joint_call.{interval}.vcf.gz",
        interval_vep      = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/joint_call.{interval}.vep.vcf.gz", 
        interval_vep_tbi  = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/joint_call.{interval}.vep.vcf.gz.tbi", 
        interval_vep_html = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/joint_call.{interval}.vep.vcf_summary.html", 
    params:
        out_name = lambda wildcards, output: os.path.splitext(output.interval_vep)[0],
        ref_fasta = config["ref_fasta"],
        ref_gtf   = config["ref_gtf"]
    threads: 6
    resources:
         time   = 720,
         mem_mb = 60000
    shell:
        '''
            set -e

            source activate ensembl-vep

            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.final_vcf} \
                -O {output.final_interval} \
                -L {input.interval}

            vep \
                -i {output.final_interval} \
                -o {params.out_name} \
                --gtf {params.ref_gtf} \
                --fasta {params.ref_fasta} \
                --fork 4 \
                --everything \
                --force_overwrite \
                --vcf \
                --dont_skip

            bgzip --threads 6 -c {params.out_name} > {output.interval_vep}
            tabix -p vcf {output.interval_vep}
        '''

def get_vep_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.generate_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"wags_{interval}.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/wags_{interval}/joint_call.{interval}.vep.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        interval=INTERVALS
    ))

rule final_gather_veps:
    input:
        get_vep_vcfs
    output:
        vep_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
        vep_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi",
    params:
        vcf_tmp = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_genotype.{ref}.TMP.gz",
        veps    = lambda wildcards, input: " --input ".join(map(str,input)),
    threads: 24
    resources:
         time   = 1440,
         mem_mb = 22000
    shell:
        '''
            set -e

            gatk --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {params.veps} \
                --output {params.vcf_tmp}

            zcat {params.vcf_tmp} | bgzip --threads 24 -c > {output.vep_vcf} &&
            tabix -p vcf {output.vep_vcf}

            rm -f {params.vcf_tmp}
        '''

