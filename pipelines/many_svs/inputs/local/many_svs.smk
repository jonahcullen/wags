
import pandas as pd

delly_samples = pd.read_csv(config['delly_samples'], header = None)[0]
gridss_samples = pd.read_csv(config['gridss_samples'], header = None)[0]
manta_samples = pd.read_csv(config['manta_samples'], header = None)[0]
smoove_samples =  pd.read_csv(config['smoove_samples'], header = None)[0]
cram_samples = pd.read_csv(config['cram_samples'], header = None)[0]

singularity: config['sif']
chrms = pd.read_csv(config['chr_list'], header = None)[0]


rule all:
  input:
    f"output/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vep.vcf.gz" 

include: 'rules/delly.smk' 
include: 'rules/manta.smk'
include: 'rules/gridss.smk'
include: 'rules/smoove.smk'
include: 'rules/survivor.smk'
include: 'rules/graphtyper.smk'

