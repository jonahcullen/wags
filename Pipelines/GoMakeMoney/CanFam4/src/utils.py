import pandas as pd

def get_fastq(wildcards):
    '''
    Get fastq files of given sample-unit.
    '''
    fastqs = units.loc[(wildcards.readgroup_name), ["fastq_1", "fastq_2"]].dropna()
    
    return {"r1": fastqs.fastq_1, "r2": fastqs.fastq_2}

def sequence_grouping(ref_dict):
    '''
    Placeholder
    '''
    
    with open(ref_dict,"r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    # write sequence group list to interval files
    seq_group = tsv_string.split("\n")
    if not os.path.isdir("results/seq_group/no_unmap"):
        os.makedirs("results/seq_group/no_unmap",exist_ok=True)
        for i,v in enumerate(seq_group):
            with open(os.path.join("results/seq_group/no_unmap",f"group_{i}.tsv"),"w") as f:
                print(v,file=f)

    # write sequence group list with unmapped to interval file so they are recalibrated as well
    tsv_string += '\n' + "unmapped"
    seq_group_unmapped = tsv_string.split("\n")
    if not os.path.isdir("results/seq_group/with_unmap"):
        os.makedirs("results/seq_group/with_unmap",exist_ok=True)
        for i,v in enumerate(seq_group_unmapped):
            with open(os.path.join("results/seq_group/with_unmap",f"group_{i}.tsv"),"w") as f:
                print(v,file=f)
    
