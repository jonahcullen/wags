
rule download_gvcfs:
    input:
        S3.remote(
            "friedlab/wgs/{breed}/{dogid}/{ref}/gvcf/{dogid}.{ref}.g.vcf.gz",
            keep_local=True
        ),
        S3.remote(
            "friedlab/wgs/{breed}/{dogid}/{ref}/gvcf/{dogid}.{ref}.g.vcf.gz.tbi",
            keep_local=True
        )
    output:
        gvcf = "friedlab/wgs/{breed}/{dogid}/{ref}/gvcf/{dogid}.{ref}.done"
    threads: 1
    resources:
         time   = 60,
         mem_mb = 6000
    shell: "touch {output.gvcf}"

rule input_list:
    input:
        gvcfs = expand(
            "friedlab/wgs/{u.breed}/{u.dogid}/{ref}/gvcf/{u.dogid}.{ref}.done",
            u=units.itertuples(),
            ref=config["ref"],
        )
    output:
        gvcf_list = "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/inputs.list"
    run:
        with open(output.gvcf_list, "w") as f:
            for i in input.gvcfs:
                dog = os.path.basename(i).split(".")[0]
                gvcf = i.replace("done","g.vcf.gz")
                f.write(dog + "\t" + gvcf + "\n") 
