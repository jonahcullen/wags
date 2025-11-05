
 
rule smoove_merge:
  input:
    expand("{svar_dir}/smoove/{sample}.smoove.{ref}.vcf.gz",
	svar_dir = config['svar_dir'],
	sample = smoove_samples,
	ref = config['ref'])
  output:
    out = "{bucket}/smoove/smoove_merged/smoove_merged.sites.vcf.gz"
  params:
    fasta = config['ref_fasta'],
    name = "smoove_merged",
    outdir = "{bucket}/smoove/smoove_merged"
  threads: 12
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      source activate smoove

      smoove merge \
        --name {params.name} \
        -f {params.fasta} \
        --outdir {params.outdir} \
        {input}
    '''

