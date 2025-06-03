
rule select_variants_to_table:
    input:            
        final_vcf  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz",
        final_tbi  = "{bucket}/wgs/{breed}/{sample_name}/{ref}/money/final_gather/{breed}_{sample_name}.{ref}.vep.vcf.gz.tbi",
        interval   = "{bucket}/wgs/{breed}/{sample_name}/{ref}/gvcf/hc_intervals/scattered/{split}-scattered.interval_list",
        common_vcf = config['common_vcf'],
       #common_tbi = config['common_vcf']+'.tbi'
    output:
        unique_vars           = "{bucket}/compare_pop/select_vars_to_table/{ref}/{split}/{breed}_{sample_name}.{ref}.{split}.unique_vars.vcf",
        rare_and_common_vars  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{split}/{breed}_{sample_name}.{ref}.{split}.rare_and_common_vars.vcf",
        rare_vars             = "{bucket}/compare_pop/select_vars_to_table/{ref}/{split}/{breed}_{sample_name}.{ref}.{split}.rare_vars.vcf",
        unique_vars_vep_split = "{bucket}/compare_pop/select_vars_to_table/{ref}/{split}/{breed}_{sample_name}.{ref}.{split}.unique_vars.vep_split.txt",
        rare_vars_vep_split   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{split}/{breed}_{sample_name}.{ref}.{split}.rare_vars.vep_split.txt"
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
            mkdir -p {params.tmp_dir}

            # all unique variants
            gatk SelectVariants \
                -L {input.interval} \
                -R {params.ref_fasta} \
                -V {input.final_vcf} \
                -disc {params.pop_vcf} \
                --tmp-dir {params.tmp_dir} \
                -O {output.unique_vars}

            # all common and  rare variants
            gatk SelectVariants \
                -L {input.interval} \
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

def get_uniq_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"{split}-scattered.interval_list"))
    SPLIT = sorted(SPLIT, key=lambda x: int(x))
    # variants to vcfs to tables
    comp_dir = Path(rules.select_variants_to_table.output.unique_vars).parent
    # return list of split intervals
    return expand(
        os.path.join(
            comp_dir,
            "{breed}_{sample_name}.{ref}.{split}.unique_vars.vcf"
        ),bucket=config["bucket"],breed=breed,sample_name=sample_name,ref=config['ref'],split=SPLIT
    )

rule merge_uniq_vcfs:
    input:
        get_uniq_vcfs
    output:
        vcf = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz",
        tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz.tbi",
    params:
        vcfs     = lambda wildcards, input: " -INPUT ".join(map(str,input)),
        java_opt  = "-Xmx2000m",
    threads: 1
    resources:
         time   = 120,
         mem_mb = 4000
    shell:
        '''
            gatk --java-options "-Xmx2000m" \
                MergeVcfs \
                --INPUT {params.vcfs} \
                --OUTPUT {output.vcf} \
                --CREATE_INDEX true
        '''

def get_rare_vcfs(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"{split}-scattered.interval_list"))
    SPLIT = sorted(SPLIT, key=lambda x: int(x))
    # variants to vcfs to tables
    comp_dir = Path(rules.select_variants_to_table.output.rare_vars).parent
    # return list of split intervals
    return expand(
        os.path.join(
            comp_dir,
            "{breed}_{sample_name}.{ref}.{split}.rare_vars.vcf"
        ),bucket=config["bucket"],breed=breed,sample_name=sample_name,ref=config['ref'],split=SPLIT
    )

rule merge_rare_vcfs:
    input:
        get_rare_vcfs
    output:
        vcf = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz",
        tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz.tbi",
    params:
        vcfs     = lambda wildcards, input: " -INPUT ".join(map(str,input)),
        java_opt  = "-Xmx2000m",
    threads: 1
    resources:
         time   = 120,
         mem_mb = 4000
    shell:
        '''
            gatk --java-options "-Xmx2000m" \
                MergeVcfs \
                --INPUT {params.vcfs} \
                --OUTPUT {output.vcf} \
                --CREATE_INDEX true
        '''

def get_uniq_txts(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"{split}-scattered.interval_list"))
    SPLIT = sorted(SPLIT, key=lambda x: int(x))
    # variants to vcfs to tables
    comp_dir = Path(rules.select_variants_to_table.output.unique_vars_vep_split).parent
    # return list of split intervals
    return expand(
        os.path.join(
            comp_dir,
            "{breed}_{sample_name}.{ref}.{split}.unique_vars.vep_split.txt"
        ),
        bucket=config["bucket"],
        breed=breed,
        sample_name=sample_name,
        ref=config['ref'],
        split=SPLIT
    )

localrules: merge_uniq_txts
rule merge_uniq_txts:
    input:
        get_uniq_txts
    output:
        "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vep_split.txt"
    threads: 1
    resources:
         time   = 120,
         mem_mb = 4000
    shell:
        '''
            LC_ALL=C sort -V -m {input} > {output}
        '''

def get_rare_txts(wildcards):
    # interval dir from split intervals
    ivals_dir = checkpoints.split_intervals.get(**wildcards).output[0]
    SPLIT, = glob_wildcards(os.path.join(ivals_dir,"{split}-scattered.interval_list"))
    SPLIT = sorted(SPLIT, key=lambda x: int(x))
    # variants to vcfs to tables
    comp_dir = Path(rules.select_variants_to_table.output.rare_vars_vep_split).parent
    # return list of split intervals
    return expand(
        os.path.join(
            comp_dir,
            "{breed}_{sample_name}.{ref}.{split}.rare_vars.vep_split.txt"
        ),
        bucket=config["bucket"],
        breed=breed,
        sample_name=sample_name,
        ref=config['ref'],
        split=SPLIT
    )

localrules: merge_rare_txts
rule merge_rare_txts:
    input:
        get_rare_txts
    output:
        "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vep_split.txt"
    threads: 1
    resources:
         time   = 120,
         mem_mb = 4000
    shell:
        '''
            LC_ALL=C sort -V -m {input} > {output}
        '''

#rule bgzip_and_tabix:
#    input:            
#        unique_vars = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf",
#        rare_vars   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf",
#    output:
#        unique_vars_gz  = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz",
#        unique_vars_tbi = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.unique_vars.vcf.gz.tbi",
#        rare_vars_gz    = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz",
#        rare_vars_tbi   = "{bucket}/compare_pop/select_vars_to_table/{ref}/{breed}_{sample_name}.{ref}.rare_vars.vcf.gz.tbi",
#    threads: 8
#    resources:
#         time   = 60,
#         mem_mb = 16000
#    shell:
#        '''
#            # bgzip and tabix the unique vcf
#            bgzip -c {input.unique_vars} > {output.unique_vars_gz}
#            tabix -p vcf {output.unique_vars_gz}
#        
#            # and the rare vcf
#            bgzip -c {input.rare_vars} > {output.rare_vars_gz}
#            tabix -p vcf {output.rare_vars_gz}
#        '''
#
