
rule sv_delly:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai",
    output:
        delly_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.delly.tmp.bcf",
    params:
        conda_env = config['conda_envs']['delly'],
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.delly.benchmark.txt"
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 60000
    shell:
        '''
            source activate {params.conda_env}

            export OMP_NUM_THREADS={threads}

            # call svs for each sv type
            delly call \
                -t {wildcards.sv_type} \
                -g {params.ref_fasta} \
                -o {output.delly_tmp} \
                {input.final_bam}
        '''

rule sv_delly_filter:
    input:
        delly_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.delly.tmp.bcf",
    output:
        sv_bcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.filter_delly.bcf.gz",
        sv_csi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.filter_delly.bcf.gz.csi"
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.delly.filter.benchmark.txt"
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 60000
    shell:
        '''
            # filter for pass and save as compressed bcf
            bcftools filter \
                -O b \
                -o {output.sv_bcf} \
                -i "FILTER == 'PASS'" \
                {input.delly_tmp}

            # index filtered output
            bcftools index {output.sv_bcf}
        '''

rule sv_delly_concat:
    input:
        sorted(expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.filter_delly.bcf.gz",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
            sv_type=["BND","DEL","DUP","INS","INV"],
        ))
    output:
        sv_gz  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz",
        sv_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz.tbi"
    params:
        vcf_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.tmp.vcf",
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly_concat.benchmark.txt"
    threads: 12
    resources:
         time   = 480,
         mem_mb = 60000
    shell:
        '''
            set -e

            bcftools concat \
                -a \
                -O v \
                -o {params.vcf_tmp} \
                {input}

            # bgzip and index
            bgzip --threads {threads} -c {params.vcf_tmp} > {output.sv_gz}
            tabix -p vcf {output.sv_gz}

            rm -f {params.vcf_tmp}
        '''

rule sv_gridss:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai",
    output:
        gridss_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.bam",
        sv_gz      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz",
        sv_tbi     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz.tbi"
    params:
        gridss_tmp  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.tmp.vcf",
        gridss_filt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.tmp.filt.vcf",
        work_dir    = lambda wildcards, output: os.path.dirname(output.sv_gz),
        conda_env   = config['conda_envs']['gridss'],
        ref_fasta   = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.sv_gridss.benchmark.txt"
    threads: 8
    resources:
         time   = 1440,
         mem_mb = 36000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e
            
            # due to an issue with certain which aliases found on some hpcs
            # need to unset which
            unset -f which

            gridss \
                -t 8 \
                -r {params.ref_fasta} \
                -o {params.gridss_tmp} \
                -a {output.gridss_bam} \
                --jvmheap 32g \
                -w {params.work_dir} \
                {input.final_bam}

            # removed -i "FILTER == '.'" as no records were returned
            # unclear if issue sample or larger...
            # filter for pass and save as uncompressed vcf
            bcftools filter \
                -O v \
                -o {params.gridss_filt} \
                {params.gridss_tmp}
            
            # bgzip and index
            bgzip --threads {threads} -c {params.gridss_filt} > {output.sv_gz}
            tabix -p vcf {output.sv_gz}
        '''

rule sv_smoove:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai",
    output:
        smoove_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.smoove.{ref}.vcf.gz", 
        smoove_csi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.smoove.{ref}.vcf.gz.csi"
    params:
        tmp_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.{ref}-smoove.genotyped.vcf.gz",
        tmp_csi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.{ref}-smoove.genotyped.vcf.gz.csi",
        out_dir   = lambda wildcards, output: os.path.dirname(output.smoove_vcf),
        base_name = (
            lambda wildcards, input: re.split(
                r'\.aligned|\.left_aligned|\.bam', 
                os.path.basename(input.final_bam)
            )[0]
        ),
        ref_fasta = config['ref_fasta'],
        conda_env = config['conda_envs']['smoove'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.sv_smoove.benchmark.txt"
    threads: 4
    resources:
         time   = 1440,
         mem_mb = 60000
    shell:
        '''
            source activate {params.conda_env}

            smoove call \
                --outdir {params.out_dir} \
                --name {params.base_name} \
                --fasta {params.ref_fasta} \
                -p {threads} \
                --genotype {input.final_bam}

            mv {params.tmp_vcf} {output.smoove_vcf}
            mv {params.tmp_csi} {output.smoove_csi}
        ''' 

rule sv_manta:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai" 
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai",
    output:
        config     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/runWorkflow.py",
        pickle     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/runWorkflow.py.config.pickle",
        sv_gz      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz",
        sv_tbi     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz.tbi",
        cand_stat  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svCandidateGenerationStats.tsv",
        graph_stat = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svLocusGraphStats.tsv",
        align_stat = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/alignmentStatsSummary.txt",
    params:
        manta_tmp  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/variants/diploidSV.vcf.gz",
        manta_filt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/variants/diploidSV.filt.vcf",
        work_dir    = lambda wildcards, output: os.path.dirname(output.sv_gz),
        conda_env   = config['conda_envs']['manta'],
        ref_fasta   = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.benchmark.txt"
    threads: 12
    resources:
         time   = 1440,
         mem_mb = 60000,
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e

            configManta.py \
                --bam {input.final_bam} \
                --reference {params.ref_fasta} \
                --runDir {params.work_dir}
           
            # cd to working dir
            cd {params.work_dir}

            ./runWorkflow.py \
                --quiet \
                -m local \
                -j {threads}

            # cd back to working dir
            cd -

            # filter for pass and save as uncompressed vcf
            bcftools filter \
                -O v \
                -o {params.manta_filt} \
                -i "FILTER == 'PASS'" \
                {params.manta_tmp}
            
            # bgzip and index
            bgzip --threads {threads} -c {params.manta_filt} > {output.sv_gz}
            tabix -p vcf {output.sv_gz}
        '''

rule sv_done:
    input:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz.tbi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz.tbi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.smoove.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/smoove/{sample_name}.smoove.{ref}.vcf.gz.csi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz.tbi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svCandidateGenerationStats.tsv",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svLocusGraphStats.tsv",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/alignmentStatsSummary.txt",
    output:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/sv.done"
    shell:
        '''
            touch {output}
        '''

