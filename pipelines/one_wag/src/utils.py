import pandas as pd

def get_fastq(wildcards):
    '''
    Get fastq files of given sample-unit.
    '''
    fastqs = units.loc[(wildcards.readgroup_name), ['fastq_1', 'fastq_2']].dropna()
    
    return {'r1': fastqs.fastq_1, 'r2': fastqs.fastq_2}

# adatped from GATKs best practices, including the human-based
# adding of a sacrificial element:w
def sequence_grouping(base,ref_dict):
    '''
    Placeholder
    '''
    d = {}    
    with open(ref_dict,"r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
                # add to d
                d[line_split[1].split("SN:")[1]] = int(line_split[2].split("LN:")[1])
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

    # write sequence group list with unmapped to interval file so they are recalibrated as well
    tsv_string += '\n' + "unmapped"
    seq_group_unmapped = tsv_string.split("\n")
    if not os.path.isdir(os.path.join(base, "seq_group", "with_unmap")):
        os.makedirs(os.path.join(base, "seq_group", "with_unmap"), exist_ok=True)
        for i,v in enumerate(seq_group_unmapped):
            with open(
                os.path.join(base, "seq_group", "with_unmap", "group_{:04d}.tsv".format(i)),
                "w"
            ) as f:
                print(v,file=f)
    
    os.makedirs(os.path.join(base, "bed_group"), exist_ok=True)
    for i,v in enumerate(tsv_string.split("\n")):
        # skip unmapped
        if 'unmapped' in v:
            continue
        # interval file name
        bed_name = os.path.join(base, "bed_group", "bed_group_{:04d}.bed".format(i))
        # only generate if not already
        if not os.path.isfile(bed_name):
            # check if multiple contigs
            split_check = v.split("\t")
            if len(split_check) > 1:
                with open(bed_name, "w") as f:
                    for j in split_check:
                        j = j.split(':')[0]
                        print(j, "0", d[j], sep="\t", file=f)
            else:
                with open(bed_name, "w") as f:
                    v = v.split(':')[0]
                    print(v, "0", d[v], sep="\t", file=f)

# functions for generating intervals suitable for the t2t throughbred assembly

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def natural_sort(l):
    return sorted(l, key=natural_sort_key)

def read_genome_dict(dict_file):
    chrom_lengths = {}
    dict_lines = []
    with open(dict_file, 'r') as file:
        lines = file.readlines()
        dict_lines = lines[:]
        for line in lines:
            if line.startswith('@SQ'):
                parts = line.split()
                chrom = parts[1].split(':')[1]
                length = int(parts[2].split(':')[1])
                chrom_lengths[chrom] = length
    return dict_lines, chrom_lengths

def process_gaps(gaps_file, chrom_lengths):
    with open(gaps_file, 'r') as file:
        lines = file.readlines()
    
    contigs = []
    output = []
    for chrom in sorted(chrom_lengths.keys(), key=natural_sort_key):
        chrom_gaps = [line.strip() for line in lines if line.split(':')[0] == chrom]
        if not chrom_gaps:
            contigs.append(
                "{}\t1\t{}\t+\tACGTmer".format(chrom, chrom_lengths[chrom])
            )
            continue

        midpoints = []
        for gap in chrom_gaps:
            positions = gap.split(':')[1]
            start, end = map(int, positions.split('-'))
            midpoint = (start + end) // 2
            midpoints.append(midpoint)
        
        previous_end = 1
        for i, midpoint in enumerate(midpoints):
            if i == 0:
                output.append("{}\t{}\t{}\t+\tACGTmer".format(chrom, previous_end, midpoint))
            else:
                output.append("{}\t{}\t{}\t+\tACGTmer".format(chrom, midpoints[i - 1] + 1, midpoint))
            previous_end = midpoint
        
        output.append("{}\t{}\t{}\t+\tACGTmer".format(chrom, previous_end + 1, chrom_lengths[chrom]))
    # add contigs for inclusion
    output.extend(contigs)

    return output

