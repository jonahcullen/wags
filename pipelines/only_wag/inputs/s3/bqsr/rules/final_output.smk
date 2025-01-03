
rule reformat_vep_split:
    input:
        vep_vcf   = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz",
        vep_split = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.{var_type}_vars.vep_split.txt",
    output:
        reform = "{bucket}/compare_pop/final_output/{ref}/{breed}_{sample_name}.{ref}.{var_type}_vars.vep_split.reform.txt",
    threads: 1
    resources:
         time   = 30,
         mem_mb = 2000
    run:
        import re
        import gzip

        def extract_csq_header(f):
            with gzip.open(f, "rt") as f_in:
                for line in f_in:
                    if line.startswith("##INFO=<ID=CSQ"):
                        mat = re.search(r"Format: (.+?)\"", line)
                        if mat:
                            return mat.group(1).split("|")
            raise ValueError("CSQ header not found in vcf")

        # get base and consequence columns directly
        base_cols = ["chrom", "pos", "ref", "alt", "ac"]
        csq_cols = extract_csq_header(input.vep_vcf)
        header = base_cols + csq_cols

        with open(input.vep_split, "r") as f_in, open(output.reform, "w") as f_out:
            print("\t".join(header), file=f_out)
            for line in f_in:
                print(line.strip(), file=f_out)

rule convert_to_excel:
    input:
        reforms = expand(
            "{bucket}/compare_pop/final_output/{ref}/{breed}_{sample_name}.{ref}.{var_type}_vars.vep_split.reform.txt",
            bucket=config['bucket'],
            breed=breed,
            sample_name=sample_name,
            ref=config["ref"],
            var_type=["unique", "rare"]
        )
    output:
        excel_sheet = "{bucket}/compare_pop/final_output/{ref}/{breed}_{sample_name}.{ref}.tables.xlsx",
    threads: 1
    resources:
         time   = 30,
         mem_mb = 32000
    run:
        import pandas as pd

        dfs = {}
        # read each input reform vep split into a df
        for i in input.reforms:
            var_type = os.path.basename(i).split(".")[-4]
            print(f"Processing file: {i}, type: {var_type}")
            try:
                df = pd.read_csv(i, sep="\t", dtype=str, engine="python")
                if df.columns[0].startswith("chrom"):
                    print(f"Headers are correctly aligned for {var_type}.")
                else:
                    raise ValueError(f"Headers misaligned in {i}. Check the file format.")
                print(f"df for {var_type} head:\n{df.head()}")
            except Exception as e:
                print(f"error reading file {i}: {e}")
    
            dfs[var_type] = df

        print(f"Loaded DataFrames: {list(dfs.keys())}")
        # write to excel tabs
        with pd.ExcelWriter(output.excel_sheet, engine="xlsxwriter") as writer:
            for var_type, df in dfs.items():
                print(f"Writing tab: {var_type} with {len(df)} rows.")
                try:
                    df.to_excel(writer, sheet_name=var_type, index=False)
                except Exception as e:
                    print(f"Error writing tab {var_type}: {e}")
                    raise

localrules: manifest_and_archive
rule manifest_and_archive:
    input:
        vep_vcf         = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz"),
        vep_vcf_tbi     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz.tbi"),
        unique_vars_gz  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz",
        unique_vars_tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz.tbi",
        rare_vars_gz    = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz",
        rare_vars_tbi   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz.tbi",
        excel_sheet     = "{bucket}/compare_pop/final_output/{ref}/{breed}_{sample_name}.{ref}.tables.xlsx",
    output:
        manifest     = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/{breed}_{sample_name}.{ref}.manifest.txt"),
        money_tar_gz = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/{breed}_{sample_name}.{ref}.K9MM.tar.gz"),
    params:
        staging_dir = "staging_tarball"
    threads: 1
    resources:
         time   = 30,
         mem_mb = 16000
    run:
        import shutil
        import tempfile
        # make the staging directory and prepae files dict
        os.makedirs(params.staging_dir, exist_ok=True)

        with tempfile.TemporaryDirectory() as staging_dir:
            # map files to descriptions
            file_map = {
                input.vep_vcf: "VEP annotated VCF for all variants",
                input.vep_vcf_tbi: "index for the VEP annotated VCF",
                input.unique_vars_gz: "VCF of unique variants",
                input.unique_vars_tbi: "index for the unique variants VCF",
                input.rare_vars_gz: "VCF of rare variants",
                input.rare_vars_tbi: "index for the rare variants VCF",
                input.excel_sheet: "excel file of unique and rare variants"
            }

            # copy files to staging directory
            for src in file_map.keys():
                dest = os.path.join(staging_dir, os.path.basename(src))
                shutil.copy(src, dest)

            # write the manifest
            manifest_path = os.path.join(staging_dir, os.path.basename(output.manifest))
            with open(manifest_path, "w") as manifest:
                manifest.write(f"Manifest for {wildcards.breed}_{wildcards.sample_name}.{wildcards.ref} money (tar)ball\n\n")
                manifest.write("Included files:\n")
                for file, description in file_map.items():
                    manifest.write(f"{os.path.basename(file)}: {description}\n")

            # create money (tar)ball
            shell(f"tar -czvf {output.money_tar_gz} -C {staging_dir} .")

            # copy the manifest to final destination
            shutil.copy(manifest_path, output.manifest)

