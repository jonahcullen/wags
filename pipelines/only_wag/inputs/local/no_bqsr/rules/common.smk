# check if vcf with only common variants (AF > 0.005) has been generated for 
# cohort VCF and generate if not
rule scatter_commons:
    output:
        acgt_ivals = "common/common_intervals/acgt.interval_list",
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

checkpoint split_commons:
    input:
        acgt_ivals = "common/common_intervals/acgt.interval_list",
    output:
        directory("common/common_intervals/scattered")
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

rule select_common_vars:
    input:
        interval = "common/common_intervals/scattered/{common_interval}-scattered.interval_list"
    output:
        interval_bed = "common/common_intervals/scattered/{common_interval}-scattered.bed",
        common_vcf   = "common/common_intervals/scattered/common_{common_interval}.vep.vcf.gz",
    params:
        pop_vcf     = config['pop_vcf'],
        allele_freq = config['allele_freq']
    threads: 4
    resources:
         time   = 2880,
         mem_mb = 16000
    shell:
        '''
            # convert to bed format
            awk 'BEGIN {{FS="\t";OFS="\t"}} !/^@/ {{print $1, $2-1, $3}}' \
                {input.interval} > {output.interval_bed}

            bcftools view -Oz \
                -R {output.interval_bed} \
                -i 'AF[*]>{params.allele_freq}' \
                {params.pop_vcf} \
                -o {output.common_vcf}

            tabix -p vcf {output.common_vcf}
        '''

def get_common_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_commons.get(**wildcards).output[0]
    # variable number of intervals up to scatter_size set in config (default: 50)
    INTERVALS, = glob_wildcards(os.path.join(ivals_dir,"{interval}-scattered.interval_list"))
    # return list of split intervals recal.vcf.gz
    return sorted(expand(
       #"{bucket}/wgs/{breed}/{sample_name}/{ref}/money/common/wags_{vep_interval}/common_{vep_interval}.vep.vcf.gz",
        "common/common_intervals/scattered/common_{common_interval}.vep.vcf.gz",
        common_interval = INTERVALS
    ))



if not os.path.isfile(config['common_vcf']):
    rule gather_common_vcfs:
        input:
            get_common_vcfs
        output:
            common_vcf = config['common_vcf'],
            common_tbi = config['common_vcf']+'.tbi'
        params:
           #vcf_tmp = 'hello',
           #vcf_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/common/common_gather/joint_genotype.{ref}.TMP.gz",
            vcf_tmp = "common/common_vars.TMP.gz",
            commons = lambda wildcards, input: " --input ".join(map(str,input)),
        threads: 12
        resources:
             time   = 2160,
             mem_mb = 32000
        shell:
            '''
                set -e

                gatk --java-options "-Xmx18g -Xms6g" \
                    GatherVcfsCloud \
                    --ignore-safety-checks \
                    --gather-type BLOCK \
                    --input {params.commons} \
                    --output {params.vcf_tmp}

                zcat {params.vcf_tmp} | bgzip --threads {threads} -c > {output.common_vcf} &&
                tabix -p vcf {output.common_vcf}

                rm -f {params.vcf_tmp}
            '''

