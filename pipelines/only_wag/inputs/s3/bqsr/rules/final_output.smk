
rule final_output:
    input:
        unique_vars_vep_split = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vep_split.txt",
        rare_vars_vep_split   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vep_split.txt"
    output:
        excel_sheet = "{bucket}/compare_pop/final_output/{ref}/{breed}_{sample_name}.{ref}.tables.xlsx",
    params:
        base_name = "{bucket}/compare_pop/final_output/{ref}/"
    threads: 1
    resources:
         time   = 30,
         mem_mb = 32000
    run:
        header = [
            'chrom', 'pos', 'ref', 'alt', 'ac',
            'allele', 'consequence', 'impact', 'symbol', 'gene', 
            'feature_type', 'feature', 'biotype', 'exon', 'intron', 
            'hgvsc', 'hgvsp', 'cdna_position', 'cds_position', 
            'protein_position', 'amino_acids', 'codons', 'existing_variation',
            'distance', 'strand', 'flags', 'variant_class', 'symbol_source', 
            'hgnc_id', 'canonical', 'mane_select', 'mane_plus_clinical', 
            'tsl', 'appris', 'ccds', 'ensp', 'swissprot', 'trembl', 'uniparc', 
            'uniprot_isoform', 'source', 'gene_pheno', 'sift', 'polyphen', 
            'domains', 'mirna', 'hgvs_offset', 'af', 'afr_af', 'amr_af', 
            'eas_af', 'eur_af', 'sas_af', 'aa_af', 'ea_af', 'gnomad_af', 
            'gnomad_afr_af', 'gnomad_amr_af', 'gnomad_asj_af', 
            'gnomad_eas_af', 'gnomad_fin_af', 'gnomad_nfe_af', 
            'gnomad_oth_af', 'gnomad_sas_af', 'max_af', 'max_af_pops', 
            'clin_sig', 'somatic', 'pheno', 'pubmed', 'check_ref', 
            'motif_name', 'motif_pos', 'high_inf_pos', 
            'motif_score_change', 'transcription_factors', 
            f'{wildcards.ref}_gene.gtf.gz'
        ]
     
        # for each table, add header and print the rest of table   
        dfs = {}
        for f in [input.unique_vars_vep_split, input.rare_vars_vep_split]:
            tmp = os.path.basename(f).replace(".txt",".reform.txt")
            reform = os.path.join(params.base_name,tmp)
            # get name of each table as key for dict
            name = tmp.split(".vep_split")[0]
            # open table, add header, and write to reform 
            with open(f, "r") as infile, open(reform, "w") as outfile:
                print("\t".join(header), file=outfile)
                for line in infile:
                    print(line.strip(), file=outfile)
            # read fixed table to dataframe
            dfs[name] = pd.read_csv(reform, sep="\t")

        # write each df to the same excel sheet
        with pd.ExcelWriter(
                output.excel_sheet,
                engine="xlsxwriter",
                options={"strings_to_formulas": False}
            ) as writer:
            for k,v in dfs.items():
                v.to_excel(
                    writer,
                    sheet_name=k[:30],
                    index=False
                )
            writer.save()

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
    threads: 1
    resources:
         time   = 30,
         mem_mb = 16000
    run:
        with open(output.manifest, "w") as outfile:
            s = textwrap.dedent(
                f'''
                    Included files:
                    1) results/vep_final_vcf/{wildcards.breed}_{wildcards.sample_name}.{wildcards.ref}.vep.vcf.gz - annotated VCF and index (.tbi) for all {wildcards.sample_name} variants
                    2) results/select_vars_to_table/{wildcards.sample_name}.{wildcards.ref}.unique_vars.vcf.gz - annotated VCF and index (.tbi) for all unique {wildcards.sample_name} variants
                    3) results/select_vars_to_table/{wildcards.sample_name}.{wildcards.ref}.rare_vars.vcf.gz - annotated VCF and index (.tbi) for all rare {wildcards.sample_name} variants
                    4) results/final_output/{wildcards.sample_name}_{wildcards.ref}.tables.xlsx - contains two sheets
                    \t i) {wildcards.sample_name}_{wildcards.ref}.unique_vars - table version of item 2
                    \t ii) {wildcards.sample_name}_{wildcards.ref}.rare_vars - table version of item 3
                '''
            )
            print(s, file=outfile)

        shell('''
            tar -czvf {output.money_tar_gz} \
                {input.vep_vcf} \
                {input.vep_vcf_tbi} \
                {input.unique_vars_gz} \
                {input.unique_vars_tbi} \
                {input.rare_vars_gz} \
                {input.rare_vars_tbi} \
                {input.excel_sheet} \
                {output.manifest}
        ''')
