
import os

localrules: ref_index,
            ref_dict

singularity: config['sif']

rule all:
    input:
        bwa_indices = expand(
            f"{config['ref_fasta']}.{{ext}}",
            ext = ["bwt","pac","ann","amb","sa"]
        )

rule ref_index:
    input:
        ref_fasta = config['ref_fasta'],
    output:
        ref_fai = f"{config['ref_fasta']}.fai"
    shell:
        '''
            samtools faidx {input.ref_fasta}
        '''

rule ref_dict:
    input:
        ref_fasta = config['ref_fasta'],
    output:
        ref_dict = f"{os.path.splitext(config['ref_fasta'])[0]}.dict"
    shell:
        '''
            gatk CreateSequenceDictionary \
                -R {input.ref_fasta}
        '''
        
rule bwa_index:
    input:
        ref_dict  = f"{os.path.splitext(config['ref_fasta'])[0]}.dict",
        ref_fai   = f"{config['ref_fasta']}.fai",
        ref_fasta = config['ref_fasta'],
    output:
        bwa_bwt = f"{config['ref_fasta']}.bwt",
        bwa_pac = f"{config['ref_fasta']}.pac",
        bwa_ann = f"{config['ref_fasta']}.ann",
        bwa_amb = f"{config['ref_fasta']}.amb",
        bwa_sa  = f"{config['ref_fasta']}.sa",
    threads: 4
    resources:
         time   = 480,
         mem_mb = 24000,
    shell:
        '''
            bwa index {input.ref_fasta}
        '''
