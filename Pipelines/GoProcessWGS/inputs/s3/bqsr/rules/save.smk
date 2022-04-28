
rule upload_fastqs:
    input:
        unpack(get_fastq)
    output:
        touch("{bucket}/fastqc/{breed}_{sample_name}/{readgroup_name}.upload")
    params:
        alias    = config["alias"],
        flowcell = lambda wildcards: units.loc[units["readgroup_name"] == wildcards.readgroup_name,"flowcell"].values[0],
    run:
        # Due to the way minio client (mc) works, the container version of mc
        # will update config files locally to match the container version. This
        # can be a problem and break a users local mc setup if the versions
        # are not compatiable. As a work around, the run directive cannot use
        # a container and thus the below shell comand will be executed using
        # the user installed and configured mc

        shell(f'''
            mc cp {input.r1} {input.r2} \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/fastq/{params.flowcell}/
        ''')

rule upload_pipe_and_logs:
    input:
        S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/qc/multiqc_report.html"),
        final_gvcf     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz"),
        final_gvcf_tbi = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/{sample_name}.{ref}.g.vcf.gz.tbi"),
    output:
        touch("{bucket}/{breed}_{sample_name}_{ref}.done"),
    params:
        alias = config["alias"],
    run:
        # Due to the way minio client (mc) works, the container version of mc
        # will update config files locally to match the container version. This
        # can be a problem and break a users local mc setup if the versions
        # are not compatiable. As a work around, the run directive cannot use
        # a container and thus the below shell comand will be executed using
        # the user installed and configured mc

        shell('''
            set -e

            mc cp --recursive ./src/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/pipeline/

            mc cp --recursive ./slurm.go_wgs/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/pipeline/

            mc cp {wildcards.ref}_config.yaml {wildcards.breed}_{wildcards.sample_name}.go_process_wgs.slurm go_process_wgs.smk input.tsv \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/pipeline/

            mc cp --recursive ./.logs/ \
                {params.alias}/{wildcards.bucket}/wgs/{wildcards.breed}/{wildcards.sample_name}/{wildcards.ref}/
        ''')
