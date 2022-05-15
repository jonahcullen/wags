
# check if reference exists outside of the container (eg custom reference) and
# if bwa index has already been run - run if not
rule bwa_index:
    output:
        "{bucket}/wgs/{breed}/{sample_name}/{ref}/bam/index.done"
    params:
        ref_fasta = config['ref_fasta'],
        ref_index = f"{config['ref_fasta']}.amb"
    threads: 4
    resources:
         time   = 120,
         mem_mb = 24000,
    shell:
        '''
            BWA_INDEX={params.ref_index}

            if [ -f "$BWA_INDEX" ]; then
                touch {output}
            else
                bwa index {params.ref_fasta}
                touch {output}
            fi
        '''
