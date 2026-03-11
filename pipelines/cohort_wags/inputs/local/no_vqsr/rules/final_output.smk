
rule add_genotype_counts:
    input:
        vep_vcf   = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
        vep_split = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.{var_type}_vars.vep_split.txt"
    output:
        reform = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.{var_type}_vars.vep_split.gt_cts.tsv",
    params:
        gt_counts = Path(workflow.basedir) / 'src' / 'gt_counts_per_csq.py',
        joint_mode = config['joint']
    singularity: '/projects/standard/fried255/cull0084/.sif/dog/UU_Cfam_GSD_1.0_ROSY.NO_SVS.sif'
    threads: 1
    resources:
         time   = 30,
         mem_mb = 32000
    shell:
        '''
            python {params.gt_counts} \
                --vcf {input.vep_vcf} \
                --vep_split {input.vep_split} \
                --output {output.reform}
        '''

vep_merge_rules = {
    "unique": rules.merge_uniq_txts,
    "rare": rules.merge_rare_txts
}

def get_vep_split(wildcards):
    if config['joint']:
        return rules.add_genotype_counts.output.reform
    else:
        return vep_merge_rules[wildcards.var_type].output[0]

##########################################
# NOTE HARDCODED THE OUTPUT NAMES WHICH NEEDS TO BE IN THE CONFIG DURING SETUP
##########################################
rule reformat_vep_split:
    input:
        vep_vcf   = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
        vep_split = get_vep_split
    output:
        reform = "{bucket}/compare_pop/final_output/{ref}/{date}/sadie_negative_joint_money.{ref}.{var_type}_vars.vep_split.reform.tsv",
    params:
        joint_mode = config['joint']
    threads: 1
    resources:
         time   = 30,
         mem_mb = 2000
    run:
        import re
        import csv
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
        base_cols = ["chrom", "pos", "locus", "ref", "alt", "ac"]
        count_cols = ["csq_allele", "csq_count", "hom_ref", "het", "hom_var", "missing"]
        csq_cols = extract_csq_header(input.vep_vcf)

        header = base_cols + (count_cols if params.joint_mode else []) + csq_cols

        with open(input.vep_split, "r") as f_in, open(output.reform, "w") as f_out:
            print("\t".join(header), file=f_out)
            for line in f_in:
                fields = line.strip().split("\t")
                locus = f"{fields[0]}:{fields[1]}"

                base_out = [fields[0], fields[1], locus, fields[2], fields[3], fields[4]]
                csq = fields[5:-6] if params.joint_mode else fields[5:]
                counts = fields[-6:] if params.joint_mode else []

                out_fields = base_out + counts + csq
                print("\t".join(out_fields), file=f_out)
               #if params.joint_mode:
               #   #base = fields[:5]
               #    base_out   = [fields[0], fields[1], locus, fields[2], fields[3], fields[4]]
               #    csq        = fields[5:-6]
               #    counts_out = fields[-6:]
               #    out_fields = base_out + counts_out + csq
               #else:
               #    reordered = fields
               #print("\t".join(reordered), file=f_out)

rule convert_to_excel:
    input:
        reforms = expand(
            "{bucket}/compare_pop/final_output/{ref}/{date}/sadie_negative_joint_money.{ref}.{var_type}_vars.vep_split.reform.tsv",
            bucket=config['bucket'],
            ref=config["ref"],
            var_type=["unique", "rare"],
            date=config['date']
        )
    output:
        excel_sheet = "{bucket}/compare_pop/final_output/{ref}/{date}/sadie_negative_joint_money.{date}.tables.xlsx",
    threads: 1
    resources:
         time   = 30,
         mem_mb = 60000
    run:
        import pandas as pd

        def drop_empties(df):
            missing_mask = df.isin(['.', '..', ''])
            drop_cols = df.columns[missing_mask.all()]
            return df.drop(columns=drop_cols)

        dfs = {}
        # read each input reform vep split into a df
        for i in input.reforms:
            var_type = os.path.basename(i).split(".")[-4]
            print(f"processing file: {i}, type: {var_type}")
            try:
                df = pd.read_csv(i, sep="\t", dtype=str, engine="python")
                if df.columns[0].startswith("chrom"):
                    print(f"headers are correctly aligned for {var_type}.")
                else:
                    raise ValueError(f"headers misaligned in {i}. Check the file format.")
                
                # drop variants in unassembled contigs
                df = df[~df['chrom'].str.startswith('chrUn_')]
                print("dropped variants on unassembled contigs")

                # drop empty columns
                orig_cols = df.shape[1]
                df = drop_empties(df)
                removed_cols = orig_cols - df.shape[1]
                print(f"dropped {removed_cols} columns with all missing values")
                
                if "csq_allele" in df.columns:
                    df = df.drop(columns=["csq_allele"])
                    print("dropped column: csq_allele")

                print(f"df for {var_type} head:\n{df.head()}")
            except Exception as e:
                print(f"error reading file {i}: {e}")
    
            dfs[var_type] = df

        print(f"loaded DataFrames: {list(dfs.keys())}")
        # write to excel tabs
        with pd.ExcelWriter(output.excel_sheet, engine="xlsxwriter") as writer:
            for var_type, df in dfs.items():
                print(f"writing tab: {var_type} with {len(df)} rows.")
                try:
                    df.to_excel(writer, sheet_name=var_type, index=False)
                except Exception as e:
                    print(f"Error writing tab {var_type}: {e}")
                    raise
##################################
# NOTE REDO THE NAMING OUTPUTS USING CONFIG DEFINED NAME
##################################
localrules: manifest_and_archive
rule manifest_and_archive:
    input:
        vep_vcf         = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz",
        vep_vcf_tbi     = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/joint_call.{ref}.{date}.vep.vcf.gz.tbi",
        unique_vars_gz  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.unique_vars.vcf.gz",
        unique_vars_tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.unique_vars.vcf.gz.tbi",
        rare_vars_gz    = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.rare_vars.vcf.gz",
        rare_vars_tbi   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{date}/joint_money.{ref}.rare_vars.vcf.gz.tbi",
        excel_sheet     = "{bucket}/compare_pop/final_output/{ref}/{date}/sadie_negative_joint_money.{date}.tables.xlsx",
    output:
        manifest     = "{bucket}/wgs/pipeline/{ref}/{date}/money/sadie_negative_joint_money.manifest.txt",
        money_tar_gz = "{bucket}/wgs/pipeline/{ref}/{date}/money/sadie_negative_joint_money.{ref}.{date}.K9MM.tar.gz",
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
                manifest.write(f"Manifest for joint money (tar)ball\n\n")
                manifest.write("Included files:\n")
                for file, description in file_map.items():
                    manifest.write(f"{os.path.basename(file)}: {description}\n")

            # create money (tar)ball
            shell(f"tar -czvf {output.money_tar_gz} -C {staging_dir} .")

            # copy the manifest to final destination
            shutil.copy(manifest_path, output.manifest)

