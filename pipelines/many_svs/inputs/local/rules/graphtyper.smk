
# graphtype gets annoyed when the reference file used to make the cram is different from the reference file
# used by graphtyper. I can get around the problem with directory of 'cache' files made by samtools.
# The grumpiness can probably be avoided if the container reference file is used.
rule graphtyper:
  input:
    cram = f"{config['cram_dir']}/{{sample}}.{config['suffix']}.cram",
    vcf = "{bucket}/survivor/filter/survivor.filter.vcf.gz"
  output:
    vcf = "{bucket}/graphtyper/individuals/{sample}/{sample}.genotype.vcf.gz"
  params:
    out_dir = "{bucket}/graphtyper/individuals/{sample}/chrms",
    fasta = config['ref_fasta'],
    region = config['chr_list'], 
    cache = config['cache']
  threads: 10
  resources:
    time = 480,
    mem_mb = 40000
  shell:
    '''
      source activate svu
      
      export REF_PATH={params.cache}/%2s/%2s/%s
      export REF_CACHE={params.cache}/%2s/%2s/%s

      graphtyper genotype_sv \
        {params.fasta} \
        {input.vcf} \
        --force_no_filter_zero_qual \
        --threads {threads} \
        --sam {input.cram} \
        --region_file {params.region} \
        --output {params.out_dir}

      set +e

      cat {params.region} \
        |while read chrom; do \
           if [[ ! -d {params.out_dir}/${{chrom}} ]]; then continue; fi; \
           find {params.out_dir}/${{chrom}} -name "*.vcf.gz" | sort ; \
         done > {params.out_dir}/vcf.list

      bcftools concat -n -f {params.out_dir}/vcf.list \
        -Oz -o {output.vcf}

      tabix {output.vcf}
    '''

rule graph_merge:
  input:
    expand("{bucket}/graphtyper/individuals/{sample}/{sample}.genotype.vcf.gz", 
      bucket = config['bucket'],
      sample = cram_samples)
  output:
    vcf = f"{{bucket}}/graphtyper/main_raw/{config['species']}.square.vcf.gz"
  threads: 1
  resources:
    time = 240,
    mem_mb = 60000
  shell:
    '''
      source activate svu 

      graphtyper vcf_merge \
        --sv \
        {input} |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule graph_filter:
  input:
    f"{{bucket}}/graphtyper/main_raw/{config['species']}.square.vcf.gz"
  output:
    f"{{bucket}}/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vcf.gz"
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  shell:
    '''
      bcftools view \
        -f "PASS" \
        -i "SVMODEL='AGGREGATED'" \
        -Oz \
        -o {output} \
        {input}

      tabix {output}
    '''

rule graph_vep:
  input:
    vcf = f"{{bucket}}/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vcf.gz"
  output:
    out_tmp = temp(f"{{bucket}}/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vep.vcf"),
    final = f"{{bucket}}/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vep.vcf.gz"
  params:
    fasta = config['ref_fasta'],
    gtf = config['ref_gtf']
  threads: 4
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      set -e
      source activate ensembl-vep

      vep \
        -i {input.vcf} \
        -o {output.out_tmp} \
        -gtf {params.gtf} \
        --fasta {params.fasta} \
        --format vcf \
        --vcf \
        --everything \
        --dont_skip \
        --fork {threads}

      bgzip -c {output.out_tmp} > {output.final}

      tabix {output.final}
    '''

