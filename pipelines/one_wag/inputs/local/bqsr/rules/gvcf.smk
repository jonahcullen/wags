
rule scatter_intervals:
    output:
        acgt_ivals = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.interval_list",
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
        '''

checkpoint split_intervals:
    input:
        acgt_ivals = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.interval_list",
    output:
        directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered")
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

rule haplotype_caller:
    input:
        final_cram = "{bucket}/wgs/{breed}/{sample_name}/{ref}/cram/{sample_name}.{ref}.cram"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/cram/{sample_name}.{ref}.left_aligned.cram",
        final_crai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/cram/{sample_name}.{ref}.cram.crai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/cram/{sample_name}.{ref}.left_aligned.cram.crai",
        interval   = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{split}-scattered.interval_list"
    output:
        hc_gvcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{sample_name}.{split}.g.vcf.gz",
        hc_tbi  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{sample_name}.{split}.g.vcf.gz.tbi"
    params:
        java_opt  = "-Xmx24G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/benchmarks/{sample_name}.{split}.hc.benchmark.txt"
    threads: 4
    resources:
         time   = 1080,
         mem_mb = 26000
    shell:
        '''
            gatk --java-options "{params.java_opt}" \
                HaplotypeCaller \
                -R {params.ref_fasta} \
                -I {input.final_cram} \
                -L {input.interval} \
                -O {output.hc_gvcf} \
                -contamination 0 -ERC GVCF
        '''

def get_gvcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    # variable number of intervals up to scatter_size set in config (default: 50)
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"{split}-scattered.interval_list"))
    # return list of split intervals
    return expand(os.path.join(ivals_dir,"{sample_name}.{split}.g.vcf.gz"),sample_name=sample_name,split=SPLIT)

rule merge_gvcfs:
    input:
        get_gvcfs
    output:
        final_gvcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz",
        final_tbi  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz.tbi",
    params:
        gvcfs     = lambda wildcards, input: " -INPUT ".join(map(str,input)),
        java_opt  = "-Xmx2000m",
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.merge_hc.benchmark.txt"
    threads: 4
    resources:
         time   = 120,
         mem_mb = 4000
    shell:
        '''
            gatk --java-options {params.java_opt} \
                MergeVcfs \
                --INPUT {params.gvcfs} \
                --OUTPUT {output.final_gvcf}
        '''

