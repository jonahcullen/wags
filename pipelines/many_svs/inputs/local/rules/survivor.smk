
rule copy_vcfs:
  input:
    "{bucket}/manta/manta_merged/manta_merged.sites.vcf.gz",
    "{bucket}/smoove/smoove_merged/smoove_merged.sites.vcf.gz",
    "{bucket}/delly/delly_merged/delly_merged.sites.vcf.gz",
    "{bucket}/gridss/gridss_merged/gridss_merged.sites.vcf.gz"
  output:
    "{bucket}/survivor/callers/manta_merged.sites.vcf",
    "{bucket}/survivor/callers/smoove_merged.sites.vcf",
    "{bucket}/survivor/callers/delly_merged.sites.vcf",
    "{bucket}/survivor/callers/gridss_merged.sites.vcf"
  threads: 1
  resources:
    time = 20,
    mem_mb = 20000
  shell:
    '''
      mkdir -p {wildcards.bucket}/survivor/callers
      for f in {input}; do
        base=$(basename $f .gz)
        cp "$f" {wildcards.bucket}/survivor/callers/
        gunzip -f {wildcards.bucket}/survivor/callers/$(basename $f)
      done 
    '''

rule survivor:
  input:
    "{bucket}/survivor/callers/manta_merged.sites.vcf",
    "{bucket}/survivor/callers/smoove_merged.sites.vcf",
    "{bucket}/survivor/callers/delly_merged.sites.vcf",
    "{bucket}/survivor/callers/gridss_merged.sites.vcf"
  output:
    "{bucket}/survivor/merged/survivor.raw.vcf"
  params:
    caller_list = "{bucket}/survivor/survivor.list",
    distance = config['sv_overlap'],
    callers = 1,
    sv_same_type = 1,
    sv_strand = 0,
    sv_dist = 0,
    sv_len = config['sv_length']
  threads: 1 
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    ''' 
      source activate svu
      
      ls {input} > {params.caller_list}

      SURVIVOR \
        merge \
        {params.caller_list} \
        {params.distance} \
        {params.callers} \
        {params.sv_same_type} \
        {params.sv_strand} \
        {params.sv_dist} \
        {params.sv_len} \
        {output}
    '''

rule filter_survivor:
  input:
    "{bucket}/survivor/merged/survivor.raw.vcf"
  output:
    vcf = "{bucket}/survivor/filter/survivor.filter.vcf" 
  params: 
    manta = "{bucket}/survivor/callers/manta_merged.sites.vcf",
    caller_number = config['sv_callers']
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  run:
    import vcf
    
    def format_line(record):
    
      info_fields = []
      for key, value in record.INFO.items():
        if isinstance(value, list):
          value_str = ','.join(map(str,value))
        else:
          value_str = str(value)
        info_fields.append('{}={}'.format(key, value_str))
    
      vcf_line = '\t'.join([
            record.CHROM,
            str(record.POS),
            record.ID if record.ID else '.',
            record.REF,
            ','.join(str(alt) for alt in record.ALT),
            str(record.QUAL) if record.QUAL is not None else '.',
            ','.join(record.FILTER) if record.FILTER else 'PASS',
            ';'.join(info_fields)])
    
      return vcf_line
    
    vcf_r = vcf.Reader(filename = str(input[0]))
    manta = params.manta
    
    with open(output.vcf, 'w') as out_vcf:
      for line in vcf_r._header_lines:
        print(line, file = out_vcf)
      print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file = out_vcf)
      
      for record in vcf_r:
        info = record.INFO
        if info['SVTYPE'] != 'TRA':
          if info['SUPP'] == '1' and info['SVTYPE'] == 'INS':
            manta_geno = record.genotype(manta)
            if manta_geno['TY'] == 'INS':
              print(format_line(record), file = out_vcf)
          elif int(info['SUPP']) >= params.caller_number:
             print(format_line(record), file = out_vcf) 

rule bgzip_vcf:
  input:
    "{bucket}/survivor/filter/survivor.filter.vcf"
  output:
    vcf = "{bucket}/survivor/filter/survivor.filter.vcf.gz",
    tbi =  "{bucket}/survivor/filter/survivor.filter.vcf.gz.tbi"
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  shell:
    '''
      bcftools sort -Oz -o {output.vcf} {input}
      tabix {output.vcf} 
    '''
