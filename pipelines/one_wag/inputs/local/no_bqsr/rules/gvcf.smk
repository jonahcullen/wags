
if "Thoroughbred" not in config['ref']:
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
else:
    localrules: t2t_scatter_intervals
    rule t2t_scatter_intervals:
        output:
            acgt_ivals = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.interval_list",
        params:
            ref_dict = config['ref_dict'],
            tb_gaps  = "/home/refgen/horse/Thoroughbred/resources/TB-T2T_gaps.csv"
        threads: 1
        run:
            dict_lines, chrom_lengths = read_genome_dict(params.ref_dict)
            new_ranges = process_gaps(params.tb_gaps, chrom_lengths)
            # apply natural sort
            sorted_new_ranges = natural_sort(new_ranges)
            # ugly swap
            chr_m_ivals = [ival for ival in sorted_new_ranges if ival.startswith('chrM')]
            sorted_chr_ivals = [ival for ival in sorted_new_ranges if not ival.startswith('chrM')]
            chr_x_index = next((i for i, v in enumerate(sorted_chr_ivals) if v.startswith('chrX')), -1)
            # NOTE - this does not quite work as chrM only appears after one of the X ivals
            # but technically does not matter as the interval gets split out but should be
            # fixed
            if chr_x_index != -1:
                sorted_chr_ivals[chr_x_index+1:chr_x_index+1] = chr_m_ivals
            else:
                sorted_chr_ivals.extend(chr_m_ivals)
            # write to output
            with open(output.acgt_ivals, 'w') as f_out:
                for line in dict_lines:
                    f_out.write(line)
                for line in sorted_chr_ivals:
                    f_out.write(line + '\n')


    localrules: t2t_split_intervals
    checkpoint t2t_split_intervals:
        input:
            acgt_ivals = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/acgt.interval_list",
        output:
            split = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered")
        params:
            ref_fasta    = config['ref_fasta'],
            scatter_size = config['scatter_size'],
        threads: 1
        resources:
             time   = 20,
             mem_mb = 8000
        run:
            os.makedirs(output.split, exist_ok=True)
            chr_ivals = []
            hap1_ivals = []
            hap2_ivals = []
            unass_ivals = []
            headers = []

            with open(input.acgt_ivals, 'r') as f_in:
                for line in f_in:
                    if line.startswith('@'):
                        headers.append(line)
                    elif line.startswith('chr'):
                        chr_ivals.append(line)
                    elif line.startswith('haplotype1'):
                        hap1_ivals.append(line)
                    elif line.startswith('haplotype2'):
                        hap2_ivals.append(line)
                    elif line.startswith('unassigned'):
                        unass_ivals.append(line)

            def write_ivals(f, ivals):
                with open(f, 'w') as f_out:
                    f_out.writelines(headers)
                    f_out.writelines(ivals)

            for i,ival in enumerate(chr_ivals):
                file_name = os.path.join(output.split, f"{i:04d}-scattered.interval_list")
                write_ivals(file_name, [ival])

            haplotype1_file = os.path.join(output.split, f"{len(chr_ivals):04d}-scattered.interval_list")
            write_ivals(haplotype1_file, hap1_ivals)
            haplotype2_file = os.path.join(output.split, f"{len(chr_ivals)+1:04d}-scattered.interval_list")
            write_ivals(haplotype2_file, hap2_ivals)
            unassigned_file = os.path.join(output.split, f"{len(chr_ivals)+2:04d}-scattered.interval_list")
            write_ivals(unassigned_file, unass_ivals)

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
    # if using a t2t assembly (eg the horse/Thoroughbred/T2T_TB.v5.fa) haplotype
    # calling will occur using known assembly gaps (generated above) where known
    if "Thoroughbred" in config['ref']:
        ivals_dir = checkpoints.t2t_split_intervals.get(**wildcards).output[0]
    else:
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

