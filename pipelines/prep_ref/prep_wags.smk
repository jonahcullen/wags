
import os

localrules: ref_index,
            ref_dict

singularity: config['sif']

rule all:
    input:
        bwa_indices = expand(
            "{}.{{ext}}".format(config['ref_fasta']),
            ext = ["bwt","pac","ann","amb","sa"]
        )

rule ref_index:
    input:
        ref_fasta = config['ref_fasta'],
    output:
        ref_fai = "{}.fai".format(config['ref_fasta'])
    shell:
        '''
            samtools faidx {input.ref_fasta}
        '''

rule ref_dict:
    input:
        ref_fasta = config['ref_fasta'],
    output:
        ref_dict = "{}.dict".format(os.path.splitext(config['ref_fasta'])[0])
    shell:
        '''
            gatk CreateSequenceDictionary \
                -R {input.ref_fasta}
        '''
        
rule bwa_index:
    input:
        ref_dict  = "{}.dict".format(os.path.splitext(config['ref_fasta'])[0]),,
        ref_fai   = "{}.fai".format(config['ref_fasta']),
        ref_fasta = config['ref_fasta'],
    output:
        "{}.bwt".format(config['ref_fasta']),
        "{}.pac".format(config['ref_fasta']),
        "{}.ann".format(config['ref_fasta']),
        "{}.amb".format(config['ref_fasta']),
        "{}.sa".format(config['ref_fasta']),
    threads: 4
    resources:
         time   = 480,
         mem_mb = 24000,
    shell:
        '''
            bwa index {input.ref_fasta}
        '''
