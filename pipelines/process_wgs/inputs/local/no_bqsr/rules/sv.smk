
# no longer tracking left aligned versus not here - could rethink...
rule sv_delly:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
        sv_bcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.{ref}.bcf.gz",
        sv_csi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.{ref}.bcf.gz.csi"
    params:
        delly_tmp = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/tmp.{sv_type}.bcf",
        conda_env = config['conda_envs']['delly'],
        ref_fasta = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.delly.benchmark.txt"
    threads: 12
    resources:
         time   = 720,
         mem_mb = 60000
    shell:
        '''
            source activate {params.conda_env}

            export OMP_NUM_THREADS={threads}

            # call svs for each sv type
            delly call \
                -t {wildcards.sv_type} \
                -g {params.ref_fasta} \
                -o {params.delly_tmp} \
                {input.final_bam}

            # filter for pass and save as compressed bcf
            bcftools filter \
                -O b \
                -o {output.sv_bcf} \
                -i "FILTER == 'PASS'" \
                {params.delly_tmp}

            # index filtered output
            bcftools index {output.sv_bcf}
        '''

rule merge_delly_calls:
    input:
        expand(
            "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.{sv_type}.{ref}.bcf.gz",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config['ref'],
            sv_type=["BND","DEL","DUP","INS","INV"],
        )
    output:
        sv_gz  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz",
        sv_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf.gz.tbi"
    params:
        sv_vcf = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.delly.{ref}.vcf"
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/delly/{sample_name}.merge_delly_calls.benchmark.txt"
    threads: 4
    resources:
         time   = 480,
         mem_mb = 60000
    shell:
        '''
            bcftools concat \
                -a \
                -O v \
                -o {params.sv_vcf} \
                {input}

            # bgzip and index
            bgzip --threads {threads} -c {params.sv_vcf} > {output.sv_gz}
            tabix -p vcf {output.sv_gz}
        '''

rule sv_gridss:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
        gridss_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.bam",
        sv_gz      = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz",
        sv_tbi     = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.vcf.gz.tbi"
    params:
        gridss_tmp  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.tmp.vcf",
        gridss_bam  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.bam",
        gridss_filt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.gridss.{ref}.tmp.filt.vcf",
        work_dir    = lambda wildcards, output: os.path.dirname(output.sv_gz),
        conda_env   = config['conda_envs']['gridss'],
        ref_fasta   = config['ref_fasta'],
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/gridss/{sample_name}.sv_gridss.benchmark.txt"
    threads: 8
    resources:
         time   = 720,
         mem_mb = 36000
    shell:
        '''
            source activate {params.conda_env}

            gridss \
                -t 8 \
                -r {params.ref_fasta} \
                -o {params.gridss_tmp} \
                -a {output.gridss_bam} \
                --jvmheap 32g \
                -w {params.work_dir} \
                {input.final_bam} \


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

rule sv_lumpy:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
        sv_gz  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.gz",
        sv_tbi = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.gz.tbi"
    params:
        lumpy_tmp  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.tmp",
        lumpy_filt = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.tmp.filt",
        lumpy_sort = "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.tmp.filt.vcf",
        work_dir   = lambda wildcards, output: os.path.basename(output.sv_gz),
        conda_env  = config['conda_envs']['lumpy'],
        ref_fasta  = config['ref_fasta'],
        ref_dict   = config['ref_dict']
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.sv_lumpy.benchmark.txt"
    threads: 8
    resources:
         time   = 2880,
         mem_mb = 36000
    shell:
        '''
            source activate {params.conda_env}

            lumpyexpress \
                -B {input.final_bam} \
                -o {params.lumpy_tmp} \
                -T {params.work_dir}

            # filter for pass and save as uncompressed vcf
            bcftools filter \
                -O v \
                -o {params.lumpy_filt} \
                -i "FILTER == '.'" \
                {params.lumpy_tmp}
           
            # sort using ref dict
            gatk SortVcf \
                -SD {params.ref_dict} \
                -I {params.lumpy_filt} \
                -O {params.lumpy_sort}

            # bgzip and index
            bgzip --threads {threads} -c {params.lumpy_sort} > {output.sv_gz}
            tabix -p vcf {output.sv_gz}
        '''

rule sv_manta:
    input:
        final_bam = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bam"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bam",
        final_bai = "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.bai"
            if not config['left_align'] else "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.bai",
    output:
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
    threads: 24
    resources:
         time   = 480,
         mem_mb = 36000
    shell:
        '''
            source activate {params.conda_env}

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
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/lumpy/{sample_name}.lumpy.{ref}.vcf.gz.tbi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/{sample_name}.manta.diploidSV.{ref}.vcf.gz.tbi",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svCandidateGenerationStats.tsv",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/svLocusGraphStats.tsv",
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/manta/results/stats/alignmentStatsSummary.txt",
    output:
        touch("{bucket}/wgs/{breed}/{sample_name}/{ref}/svar/sv.done")


