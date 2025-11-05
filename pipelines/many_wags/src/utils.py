
import re

# https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
def natural_sort(list, key=lambda s:s):
    """
    Sort the list into natural alphanumeric order.
    """
    def get_alphanum_key_func(key):
        convert = lambda text: int(text) if text.isdigit() else text 
        return lambda s: [convert(c) for c in re.split('([0-9]+)', key(s))]
    sort_key = get_alphanum_key_func(key)
    list.sort(key=sort_key)

def collapse_intervals(l, target_length=10e6):
    """
    Collapse list of intervals into sizes less than the target length.
    """
    tmp = []
    full = []
    curr_length = 0

    # chr with single interval - return interval without length
    if len(l) == 1:
        return [l[0][:-1]]

    for i in range(len(l)-1):
        # add first interval start if empty
        if not tmp:
            tmp.append(l[i][0])

        curr_length += l[i][2]
        
        if curr_length < target_length:
            tmp.append(l[i][1])
            
            if l[i+1][2] > target_length:
                # close the current interval
                tmp.append(l[i][1])
                full.append([tmp[0], tmp[-1]])
                # reset
                curr_length = 0
                tmp = []
            
            # if second to last interval, add last and end
            if i+2 == len(l):
                if tmp:
                    tmp.append(l[i+1][1])
                    full.append([tmp[0], tmp[-1]])
                # add last interval
                else:   
                    full.append([l[i+1][0], l[i+1][1]])
                
        elif curr_length > target_length:
            tmp.append(l[i][1])
            full.append([tmp[0], tmp[-1]])
            # reset
            curr_length = 0
            tmp = []

            if i+2 == len(l):
                # add last interval
                full.append([l[i+1][0], l[i+1][1]])
                
    return full

def read_genome_dict(f):
    """
    Read in genome dictionary and return chr:length dict
    """
    d = {}
   	# read in genome dict and extract chrom names and lengths 
    with open(f, 'r') as file:
        for line in file:
            if line.startswith('@SQ'):
                parts = line.strip().split('\t')
                chrom_name = None
                chrom_length = None
                for part in parts:
                    if part.startswith('SN:'):
                        chrom_name = part[3:]  # remove 'SN:'
                    elif part.startswith('LN:'):
                        chrom_length = int(part[3:])  # remove 'LN:' and convert to int
                
                # add to dictionary if both name and length were found
                if chrom_name and chrom_length:
                    d[chrom_name] = chrom_length
    
    return d

def create_bed_from_anchors(df, d, bed_out):
    """
    Create a bed file from the identified anchors
    """
    bed_records = []
    
    # process each chromosome
    for chrom, group in df.groupby("chrom"):
        # Sort anchors within chromosome
        anchors = sorted(group["anchor"].tolist())
        chrom_length = d.get(chrom, 0)
        
        # first interval - from start of chromosome to first anchor
        if anchors:
            bed_records.append((chrom, 0, anchors[0]))
            
            # middle intervals between consecutive anchors
            for i in range(len(anchors) - 1):
                bed_records.append((chrom, anchors[i], anchors[i + 1]))
            
            # last interval - from last anchor to end of chromosome
            bed_records.append((chrom, anchors[-1], chrom_length))
    
    # write to bed file
    with open(bed_out, 'w') as f_out:
        for record in bed_records:
            f_out.write(f"{record[0]}\t{record[1]}\t{record[2]}\n")

