
rule combine_snps_nonsnps:
    input:
        snp_filtered_vcf    = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/snp_fltr.vcf.gz",
        snp_filtered_tbi    = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/snp_fltr.vcf.gz.tbi",
        nonsnp_filtered_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/nonsnp_fltr.vcf.gz",
        nonsnp_filtered_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/hardflt_vcf/nonsnp_fltr.vcf.gz.tbi"
    output:
        final_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz",
        final_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz.tbi",
    threads: 12
    resources:
         time   = 600,
         mem_mb = 24000
    shell:
        '''
            bcftools concat \
                --threads {threads} \
                --allow-overlaps \
                {input.snp_filtered_vcf} {input.nonsnp_filtered_vcf} \
                -Oz \
                -o {output.final_vcf}

            tabix -p vcf {output.final_vcf}
        '''

rule vep_by_interval:
    input:
        final_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz",
        final_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vcf.gz.tbi",
        interval  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{vep_interval}-scattered.interval_list"
    output:
        final_interval    = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/split/wags_{vep_interval}/{sample_name}.{vep_interval}.vcf.gz",
        interval_vep      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/vep/wags_{vep_interval}/{sample_name}.{vep_interval}.vep.vcf.gz",
        interval_vep_tbi  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/vep/wags_{vep_interval}/{sample_name}.{vep_interval}.vep.vcf.gz.tbi",
        interval_vep_html = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/vep/wags_{vep_interval}/{sample_name}.{vep_interval}.vep.vcf_summary.html",
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
            set +eu
            source activate ensembl-vep
            set -e

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
                --fork {threads} \
                --everything \
                --force_overwrite \
                --vcf \
                --dont_skip

            bgzip --threads {threads} -c {params.out_name} > {output.interval_vep}
            tabix -p vcf {output.interval_vep}
        '''
           #vep \
           #    -i {output.final_interval} \
           #    -o {params.out_name} \
           #    --gtf {params.ref_gtf} \
           #    --fasta {params.ref_fasta} \
           #    --fork {threads} \
           #    --force_overwrite \
           #    --vcf \
           #    --custom file={params.ref_gtf},short_name=UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat,format=gtf \
           #    --custom file={params.phylop},short_name=PhyloP_score,format=bed,type=exact,coords=0 \
           #    --dont_skip \
           #    --protein \
           #    --variant_class \
           #    --biotype

def get_vep_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals up to scatter_size set in config (default: 50)
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{interval}-scattered.interval_list"))
    # return list of split intervals recal.vcf.gz
    return sorted(expand(
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/vep/wags_{vep_interval}/{sample_name}.{vep_interval}.vep.vcf.gz",
        bucket = config['bucket'],
        breed = breed,
        sample_name = sample_name,
        ref = config['ref'],
        vep_interval = INTERVALS
    ))

rule final_gather_veps:
    input:
        get_vep_vcfs
    output:
        vep_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz",
        vep_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz.tbi",
    params:
        vcf_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/joint_genotype.{ref}.TMP.gz",
        veps    = lambda wildcards, input: " --input ".join(map(str,input)),
    threads: 12
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

            zcat {params.vcf_tmp} | bgzip --threads {threads} -c > {output.vep_vcf} &&
            tabix -p vcf {output.vep_vcf}

            rm -f {params.vcf_tmp}
        '''

