
rule gridss_filter:
  input:
    vcf = f"{config['svar_dir']}/gridss/{{sample}}.gridss.{{ref}}.vcf.gz"
  output:
    vcf = "{bucket}/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      bcftools view -f "PASS" {input.vcf} > {output.vcf}
      sed -i '100i\##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV_length">' {output.vcf}
    '''

rule gridss_anno:
  input:
    vcf = "{bucket}/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
  output:
    vcf = "{bucket}/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz"
  params:
    tmp = "{bucket}/gridss/anno/{sample}.gridss.anno.tmp.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      ./src/simple-event-annotation.R {input.vcf} {params.tmp}

      bcftools view \
        -e 'SVTYPE="CTX"' \
        -Oz \
        -o {output.vcf} \
        {params.tmp}

      tabix {output.vcf}
    '''

rule make_gridss_list:
  input:
    expand("{bucket}/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz", 
      bucket = config['bucket'],
      sample = gridss_samples, 
      ref = config['ref'])
  output:
    "output/vcf_lists/gridss.vcf.list"
  shell:
    "ls {input} > {output}"

#needs help with svimmer
rule gridss_svimmer:
  input:
    "{bucket}/vcf_lists/gridss.vcf.list"
  output:
    "{bucket}/gridss/gridss_merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  threads: 6
  resources:
    time = 240,
    mem_mb = 40000
  shell:
    '''
      source activate svu

      ./src/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output}
    '''

rule gridss_join:
  input:
    expand("{bucket}/gridss/gridss_merged/separate/{chrm}.svimmer.vcf", 
      bucket = config['bucket'],
      chrm = chrms)
  output:
    "{bucket}/gridss/gridss_merged/gridss_merged.sites.vcf.gz"
  threads: 1
  resources:
    time = 120,
    mem_mb = 20000
  shell:
    '''
      bcftools concat \
        -Oz \
        -o {output} \
        {input}
    '''

