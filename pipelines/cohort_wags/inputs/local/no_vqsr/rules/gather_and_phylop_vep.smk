# generate vep intervals by runs of missing bases
rule vep_scatter_intervals:
    output:
        acgt_ivals = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/vep/acgt.interval_list"
    params:
        ref_fasta = config['ref_fasta'],
        contig_ns = config['nrun_length'],
    threads: 1
    resources:
         time   = 20,
         mem_mb = 8000
    shell:
        '''
            java -jar /opt/wags/src/picard.jar \
                ScatterIntervalsByNs \
                R={params.ref_fasta} \
                OT=ACGT \
                N={params.contig_ns} \
                O={output.acgt_ivals}
            
            sed -i '/^chrUn/d' {output.acgt_ivals}
        '''

checkpoint split_intervals:
    input:
        acgt_ivals = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/vep/acgt.interval_list"
    output:
        directory("{bucket}/wgs/pipeline/{ref}/{date}/intervals/vep/scattered")
    params:
        ref_fasta    = config['ref_fasta'],
        scatter_size = config['scatter_size'],
    threads: 1
    resources:
         time   = 20,
         mem_mb = 8000
    shell:
        '''
            gatk SplitIntervals \
                -R {params.ref_fasta} \
                -L {input.acgt_ivals} \
                --scatter-count {params.scatter_size} \
                --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
                -O {output}
        '''

rule combine_snps_nonsnps:
    input:
        snp_filtered_vcf    = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/snp_fltr.vcf.gz",
        snp_filtered_tbi    = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/snp_fltr.vcf.gz.tbi",
        nonsnp_filtered_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz",
        nonsnp_filtered_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/hardflt_vcf/nonsnp_fltr.vcf.gz.tbi"
    output:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz",
        final_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi",
    threads: 12
    resources:
         time   = 2880,
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

rule pre_vep_select_variants:
    input:
        final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz",
        final_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi",
        interval  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/vep/scattered/{vep_interval}-scattered.interval_list"
    output:
        final_ival = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.vcf.gz",
    threads: 6
    resources:
         time   = 1440,
         mem_mb = 24000
    shell:
        '''
            gatk --java-options "-Xmx12g -Xms3g" \
                SelectVariants \
                -V {input.final_vcf} \
                -O {output.final_ival} \
                -L {input.interval}
        '''

rule SIFT_vep_by_interval:
    input:
        final_ival = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.vcf.gz",
    output:
        ival_vep      = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz", 
        ival_vep_tbi  = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz.tbi", 
        ival_vep_html = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf_summary.html", 
    params:
        out_name = lambda wildcards, output: os.path.splitext(output.ival_vep)[0],
        ref_fasta = config["ref_fasta"],
        ref_gtf   = config["ref_gtf"],
        phylop    = config["phylop"],
       #sift      = config["sift"]
    threads: 6
    resources:
         time   = 2880,
         mem_mb = 60000
    shell:
        '''
            set +eu
            source activate ensembl-vep
            set -eu

            vep \
                -i {input.final_ival} \
                -o {params.out_name} \
                --fasta {params.ref_fasta} \
                --fork {threads} \
                --force_overwrite \
                --vcf \
                --custom file={params.ref_gtf},short_name=UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat,format=gtf \
                --custom file={params.phylop},short_name=PhyloP_score,format=bed,type=exact,coords=0 \
                --dont_skip \
                --protein \
                --variant_class \
                --biotype

            bgzip --threads {threads} -c {params.out_name} > {output.ival_vep}
            tabix -p vcf {output.ival_vep}
        '''
               #--plugin TranscriptAnnotator,file={params.sift},prefix=0 \
               #--dir_plugin /opt/wags/src/VepPlugins \

def get_vep_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{vep_interval}-scattered.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        vep_interval=INTERVALS
    ))

rule final_gather_veps:
    input:
        get_vep_vcfs
    output:
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi",
        vep_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
    params:
        veps = lambda wildcards, input: " --input ".join(map(str,input)),
    threads: 1
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
                --output {output.vep_vcf}

            tabix -p vcf {output.vep_vcf}
        '''

rule gbindex_vep_vcf:
    input:
        vep_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep_exact.vcf.gz",
        vep_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep_exact.vcf.gz.tbi",
    output:
        "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep_exact.vcf.gz.covtsf",
    params:
        ref_dir = os.path.dirname(config['ref_fasta']),
    threads: 4
    resources:
        time   = 2880,
        mem_mb = 12000
    shell:
        '''
            gautil coverage \
                {input.vep_vcf} \
                --refFolder={params.ref_dir}
        '''

