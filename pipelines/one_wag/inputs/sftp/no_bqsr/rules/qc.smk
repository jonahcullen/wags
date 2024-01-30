
rule fastqc:
    input:
        unpack(get_fastq)
    output:
        r1_zip  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R1_fastqc.zip"),
        r1_html = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R1_fastqc.html"),
        r2_zip  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R2_fastqc.zip"),
        r2_html = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}_R2_fastqc.html"),
    params:
        outdir = "{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}",
        r1_html = lambda wildcards, input: os.path.basename(re.sub(r'\.(fastq|fq)\.gz$', '_fastqc.html', input.r1)),
        r1_zip  = lambda wildcards, input: os.path.basename(re.sub(r'\.(fastq|fq)\.gz$', '_fastqc.zip', input.r1)),
        r2_html = lambda wildcards, input: os.path.basename(re.sub(r'\.(fastq|fq)\.gz$', '_fastqc.html', input.r2)),
        r2_zip  = lambda wildcards, input: os.path.basename(re.sub(r'\.(fastq|fq)\.gz$', '_fastqc.zip', input.r2)),
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/fastqc/{sample_name}/{flowcell}/{readgroup_name}.fastqc.benchmark.txt"
    threads: 6
    resources:
         time   = 360,
         mem_mb = 6000, 
    shell:
        '''
            set -e

            fastqc -t {threads} \
                --outdir {params.outdir} \
                {input.r1} {input.r2}

            mv {params.outdir}/{params.r1_html} {output.r1_html}
            mv {params.outdir}/{params.r2_html} {output.r2_html}
            mv {params.outdir}/{params.r1_zip} {output.r1_zip}
            mv {params.outdir}/{params.r2_zip} {output.r2_zip}
        '''

rule flagstat:
    input:
        sorted_bam = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam")
            if not config['left_align'] else SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam"),
        sorted_bai = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai")
            if not config['left_align'] else SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai"),
    output:
        flagstat = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.flagstat.txt"),
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.flagstat.benchmark.txt"
    threads: 8
    resources:
         time   = 60,
         mem_mb = 12000
    shell:
        '''
            samtools flagstat \
                {input.sorted_bam} \
                --thread 8 \
                > {output.flagstat}
        '''

rule qualimap_bamqc:
    input:
        sorted_bam = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bam")
            if not config['left_align'] else SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bam"),
        sorted_bai = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.aligned.duplicate_marked.sorted.bai")
            if not config['left_align'] else SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/{sample_name}.{ref}.left_aligned.duplicate_marked.sorted.bai"),
    output:
        report_pdf  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/report.pdf",
        report_html = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/qualimapReport.html",
        report_txt  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/genome_results.txt"),
        raw_dat     = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/raw_data_qualimapReport/"),
        figures     = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/images_qualimapReport/"),
        css         = directory("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/css/"),
    params:
        out_dir = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/",
    benchmark:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/{sample_name}.multimap_bamqc.benchmark.txt"
    threads: 12
    resources:
         time   = 120,
         mem_mb = 60000
    shell:
        '''
            unset DISPLAY

            qualimap bamqc \
                -bam {input.sorted_bam} \
                -outdir {params.out_dir} \
                -outformat PDF:HTML \
                -nt {threads} \
                --java-mem-size=60G
        '''

rule multiqc:
    input:
        r1_zip      = SFTP.remote(expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/fastqc/{{sample_name}}/{u.flowcell}/{u.readgroup_name}_R1_fastqc.zip",
            u=units.itertuples()
        )),
        r2_zip      = SFTP.remote(expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/fastqc/{{sample_name}}/{u.flowcell}/{u.readgroup_name}_R2_fastqc.zip",
            u=units.itertuples()
        )),
        flagstat    = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/{sample_name}.{ref}.flagstat.txt"),
        metrics     = expand(
            "{{bucket}}/wgs/{{breed}}/{{sample_name}}/{{ref}}/bam/{u.readgroup_name}.mark_adapt.metrics.txt",
            u=units.itertuples()
        ),
        report_txt  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/qualimap/genome_results.txt"),
        raw_dat     = rules.qualimap_bamqc.output.raw_dat,
    output: 
        SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html"),
        mqc_log    = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc.log"),
        mqc_bamqc  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt"),
        mqc_dups   = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_picard_dups.txt"),
       #mqc_adapt  = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_picard_mark_illumina_adapters.txt"),
        mqc_flag   = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_samtools_flagstat.txt"),
        mqc_fastqc = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_fastqc.txt"),
        mqc_gen    = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_general_stats.txt"),
        mqc_src    = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_sources.txt"),
        mqc_json   = SFTP.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_data/multiqc_data.json"),
    params:
        outdir = "{bucket}/wgs/{breed}/{sample_name}/{ref}/qc"
    shell:
        '''
            multiqc {wildcards.bucket} \
                -x '*duplicates_marked*' -x '*group_*' \
                --interactive \
                --force \
                -o {params.outdir}
        '''

