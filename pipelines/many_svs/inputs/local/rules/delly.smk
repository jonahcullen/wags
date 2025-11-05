
rule delly_merge:
  input:
    expand("{svar_dir}/delly/{sample}.delly.{ref}.vcf.gz",
      svar_dir = config['svar_dir'],
      sample = delly_samples,
      ref = config['ref'])
  output:
    vcf = "{bucket}/delly/delly_merged/delly_merged.sites.vcf.gz"
  params:
    tmp = "{bucket}/delly/delly_merged/delly_merged.sites.bcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
    source activate delly

    delly merge {input} \
      -o {params.tmp}

    bcftools sort -Oz -o {output.vcf} {params.tmp}

    tabix {output.vcf}
    '''
