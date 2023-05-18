
def get_recal_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.generate_intervals.get(**wildcards).output[0]
    # variable number of intervals 
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"wags_{interval}.interval_list"))
    # return list of recal vcfs
    return sorted(expand(
        "{bucket}/wgs/pipeline/{ref}/{date}/var_recal/apply/wags_{interval}/recal.{interval}.vcf.gz",
        bucket=config['bucket'],
        ref=config['ref'],
        date=config['date'],
        interval=INTERVALS
    ))

rule final_gather_vcfs:
    input:
        get_recal_vcfs
    output:
        final_vcf = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz"),
        final_tbi = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi")
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
                --output {output.final_vcf}

            gatk --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {output.final_vcf}
        '''

# generate vep intervals by runs of missing bases
rule scatter_intervals:
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

rule vep_by_interval:
    input:
        final_vcf = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz"),
        final_tbi = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vcf.gz.tbi"),
        interval  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/vep/scattered/{vep_interval}-scattered.interval_list"
    output:
        final_interval    = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/split/vep_{vep_interval}/joint_call.{vep_interval}.vcf.gz",
        interval_vep      = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz", 
        interval_vep_tbi  = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf.gz.tbi", 
        interval_vep_html = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/vep/vep_{vep_interval}/joint_call.{vep_interval}.vep.vcf_summary.html", 
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
                --fork {threads} \
                --everything \
                --force_overwrite \
                --vcf \
                --dont_skip

            bgzip --threads {threads} -c {params.out_name} > {output.interval_vep}
            tabix -p vcf {output.interval_vep}
        '''

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
        vep_vcf = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz"),
        vep_tbi = S3.remote("{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi"),
    params:
        vcf_tmp = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.TMP.gz",
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

            zcat {params.vcf_tmp} | bgzip --threads {threads} -c > {output.vep_vcf} &&
            tabix -p vcf {output.vep_vcf}

            rm -f {params.vcf_tmp}
        '''

