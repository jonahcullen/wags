
rule download_gvcfs:
    input:
        S3.remote(
            "{bucket}/wgs/{breed}/{sample}/{ref}/gvcf/{sample}.{ref}.g.vcf.gz",
            keep_local=True
        ),
        S3.remote(
            "{bucket}/wgs/{breed}/{sample}/{ref}/gvcf/{sample}.{ref}.g.vcf.gz.tbi",
            keep_local=True
        )
    output:
        gvcf = "{bucket}/wgs/{breed}/{sample}/{ref}/gvcf/{sample}.{ref}.done"
    threads: 1
    resources:
         time   = 60,
         mem_mb = 6000
    shell: "touch {output.gvcf}"

rule input_list:
    input:
        gvcfs = expand(
            "{bucket}/wgs/{u.breed}/{u.sample}/{ref}/gvcf/{u.sample}.{ref}.done",
            u=units.itertuples(),
            bucket=config['bucket'],
            ref=config['ref'],
        )
    output:
        gvcf_list = "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/inputs.list"
    run:
        with open(output.gvcf_list, "w") as f:
            for i in input.gvcfs:
                dog = os.path.basename(i).split(".")[0]
                gvcf = i.replace("done","g.vcf.gz")
                f.write(dog + "\t" + gvcf + "\n") 
