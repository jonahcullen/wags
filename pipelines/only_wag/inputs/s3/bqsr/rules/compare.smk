# check if vcf with only common variants (AF > 0.005) has been generated for 
# cohort VCF and generate if not
if not os.path.isfile(config['common_vcf']):
    rule select_common_vars:
        output:
            common_vcf = config['common_vcf'],
        params:
            pop_vcf     = config['pop_vcf'],
            allele_freq = config['allele_freq']
        threads: 8
        resources:
             time   = 2160,
             mem_mb = 32000
        shell:
            '''
                bcftools view -Oz \
                    -i 'AF[*]>{params.allele_freq}' \
                    {params.pop_vcf} \
                    -o {output.common_vcf}

                tabix -p vcf {output.common_vcf}
            '''

rule select_variants_to_table:
    input:            
        final_vcf  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz"),
        final_tbi  = S3.remote("{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz.tbi"),
        common_vcf = config['common_vcf'],
    output:
        unique_vars           = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf",
        rare_and_common_vars  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_and_common_vars.vcf",
        rare_vars             = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf",
        unique_vars_vep_split = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vep_split.txt",
        rare_vars_vep_split   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vep_split.txt"
    params:
        tmp_dir   = config['tmp_dir']['select_variants_to_table'],
        ref_fasta = config['ref_fasta'],
        pop_vcf   = config['pop_vcf']
    threads: 4
    resources:
         time   = 2880,
         mem_mb = 16000
    shell:
        '''
            mkdir -p {params.tmp_dir} results/select_vars_to_table

            # all unique variants
            gatk SelectVariants \
                -R {params.ref_fasta} \
                -V {input.final_vcf} \
                -disc {params.pop_vcf} \
                --tmp-dir {params.tmp_dir} \
                -O {output.unique_vars}

            # all common and  rare variants
            gatk SelectVariants \
                -R {params.ref_fasta} \
                -V {input.final_vcf} \
                -disc {input.common_vcf} \
                -O {output.rare_and_common_vars}

            # rare only variants
            gatk SelectVariants \
                -R {params.ref_fasta} \
                -V {output.rare_and_common_vars} \
                -disc {output.unique_vars} \
                -O {output.rare_vars}

            # vars to table - all unique
            bcftools +split-vep \
                {output.unique_vars} \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
                -d -A tab \
                -o {output.unique_vars_vep_split}

            # vars to table - rare only
            bcftools +split-vep \
                {output.rare_vars} \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%CSQ\n' \
                -d -A tab \
                -o {output.rare_vars_vep_split}
        '''

rule bgzip_and_tabix:
    input:            
        unique_vars = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf",
        rare_vars   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf",
    output:
        unique_vars_gz  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz",
        unique_vars_tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz.tbi",
        rare_vars_gz    = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz",
        rare_vars_tbi   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz.tbi",
    threads: 8
    resources:
         time   = 60,
         mem_mb = 16000
    shell:
        '''
            # bgzip and tabix the unique vcf
            bgzip -c {input.unique_vars} > {output.unique_vars_gz}
            tabix -p vcf {output.unique_vars_gz}
        
            # and the rare vcf
            bgzip -c {input.rare_vars} > {output.rare_vars_gz}
            tabix -p vcf {output.rare_vars_gz}
        '''

