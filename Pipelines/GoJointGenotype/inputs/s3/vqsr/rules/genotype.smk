
rule import_gvcfs:
    input:
       #gvcf_list = "gvcfs.list",
        gvcf_list = "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/inputs.list",
        interval  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/import/wags_{interval}.interval_list",
    output:
        directory("{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/wags_{interval}")
    params:
        tmp_dir = lambda wildcards: f"/dev/shm/{wildcards.interval}",
        batch   = config["batch_size"],
    threads: 6
    resources:
         time   = 4800,
         mem_mb = 60000
    shell:
        ''' 
            mkdir -p {params.tmp_dir}

            gatk --java-options "-Xms50g -Xmx50g" \
            GenomicsDBImport \
                --genomicsdb-workspace-path {output} \
                --batch-size {params.batch} \
                -L {input.interval} \
                --sample-name-map {input.gvcf_list} \
                --tmp-dir {params.tmp_dir} \
                --genomicsdb-shared-posixfs-optimizations true \
                --reader-threads 5 \
                -ip 500
        '''

rule genotype_gvcfs:
    input:
        ival_db  = "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/wags_{interval}",
        interval = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/import/wags_{interval}.interval_list"
    output:
        vcf = "{bucket}/wgs/pipeline/{ref}/{date}/genotype_gvcfs/wags_{interval}/output.vcf.gz",
    params:
        ref_fasta = config["ref_fasta"]
    threads: 6
    resources:
         time   = 4800,
         mem_mb = 60000
    shell:
        '''
            gatk --java-options "-Xmx50g -Xms50g" \
            GenotypeGVCFs \
                -R {params.ref_fasta} \
                -O {output.vcf} \
                -G StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V gendb://{input.ival_db} \
                -L {input.interval}
        '''
