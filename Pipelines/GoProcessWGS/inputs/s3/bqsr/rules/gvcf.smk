
rule hc_intervals:
    output:
        acgt_ivals  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.N50.interval_list",
        split_ivals = expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/{{ref}}/gvcf/hc_intervals/scattered/00{split}-scattered.interval_list",
            split=list(map("{:02d}".format, list(range(0,config['scatter_size']))))
        )
    params:
        contig_ns    = config['contig_n_size'],
        scatter_size = config['scatter_size'],
        ref_fasta    = config['ref_fasta'],
        split_dir    = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered"
    threads: 1
    resources:
         time   = 20,
         mem_mb = 8000
    shell:
        '''
            set -e

            java -jar /opt/wags/src/picard.jar \
                ScatterIntervalsByNs \
                R={params.ref_fasta} \
                OT=ACGT \
                N={params.contig_ns} \
                O={output.acgt_ivals}

            gatk SplitIntervals \
                -R {params.ref_fasta} \
                -L {output.acgt_ivals} \
                --scatter-count {params.scatter_size} \
                --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
                -O {params.split_dir}
        '''

rule haplotype_caller:
    input:
        final_bam = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam")
            if not config['left_align'] else S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam"),
        final_bai = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai")
            if not config['left_align'] else S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai"),
        interval  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/00{split}-scattered.interval_list"
    output:
        hc_gvcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{sample_name}.00{split}.g.vcf.gz"
    params:
        java_opt  = "-Xmx10G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10",
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/benchmarks/{sample_name}.00{split}.hc.benchmark.txt"
    threads: 4
    resources:
         time   = 720,
         mem_mb = 8000
    shell:
        '''
            gatk --java-options "{params.java_opt}" \
                HaplotypeCaller \
                -R {params.ref_fasta} \
                -I {input.final_bam} \
                -L {input.interval} \
                -O {output.hc_gvcf} \
                -contamination 0 -ERC GVCF
        '''

rule merge_gvcfs:
    input:
        hc_gvcfs = sorted(
            expand(
                "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{sample_name}.00{split}.g.vcf.gz",
                bucket=config["bucket"],
                ref=config['ref'],
                breed=breed,
                sample_name=sample_name,
                split=list(map("{:02d}".format, list(range(0,config['scatter_size']))))
                ), key=lambda item: int(os.path.basename(item).split(".")[1])
        )
    output:
        final_gvcf     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz"),
        final_gvcf_tbi = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz.tbi"),
    params:
        gvcfs     = lambda wildcards, input: " -INPUT ".join(map(str,input.hc_gvcfs)),
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
            gatk --java-options {params.java_opt}  \
                MergeVcfs \
                --INPUT {params.gvcfs} \
                --OUTPUT {output.final_gvcf}
        '''

