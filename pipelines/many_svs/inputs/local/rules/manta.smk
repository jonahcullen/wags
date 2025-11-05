
rule manta_inv:
  input:
    vcf = f"{config['svar_dir']}/manta/{{sample}}.manta.diploidSV.{{ref}}.vcf.gz"
  output:
    vcf = "{bucket}/manta/inv/{sample}.manta.inv.{ref}.vcf.gz"
  params:
    fasta = config['ref_fasta'] 
  threads: 1
  resources:
    time = 60,
    mem_mb = 20000
  shell:
    ''' 
      set +eu
      source activate manta
      set -e

      convertInversion.py \
        /opt/conda/bin/samtools \
        {params.fasta} \
        {input.vcf} |sed '/^\[/d' |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule make_manta_list:
  input:
    expand("{bucket}/manta/inv/{sample}.manta.inv.{ref}.vcf.gz", 
      bucket = config['bucket'],
      sample = manta_samples, 
      ref = config['ref'])
  output:
    "{bucket}/vcf_lists/manta.vcf.list"
  shell:
    "ls {input} > {output}"

#this needs help, svimmer should be in the container
rule manta_svimmer:
  input:
    "{bucket}/vcf_lists/manta.vcf.list"
  output:
    vcf = "{bucket}/manta/manta_merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}" 
  threads: 6
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      source activate svu

      ./src/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output.vcf}
    '''

rule manta_join:
  input:
    expand("{bucket}/manta/manta_merged/separate/{chrm}.svimmer.vcf", 
      bucket = config['bucket'],
      chrm = chrms)
  output:
    vcf = "{bucket}/manta/manta_merged/manta_merged.sites.vcf.gz"
  threads: 2
  resources:
    time = 120,
    mem_mb = 20000
  shell:
    '''
      bcftools concat \
        -O z \
        -o {output.vcf} \
        {input}

      tabix {output.vcf}
    '''

