
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
    Create a BED file from anchors. df has columns: chrom, anchor
    d maps chrom -> chromosome length.
    """
    import csv, os
    bed_records = []

    for chrom, group in df.groupby("chrom"):
        c = str(chrom).strip()
        anchors = sorted(int(a) for a in group["anchor"].tolist())
        clen = int(d.get(c, 0))
        if not anchors or clen <= 0:
            continue

        # [0, first), [a[i], a[i+1)), [last, clen)
        bed_records.append((c, 0, anchors[0]))
        for i in range(len(anchors) - 1):
            bed_records.append((c, anchors[i], anchors[i + 1]))
        bed_records.append((c, anchors[-1], clen))

    os.makedirs(os.path.dirname(bed_out), exist_ok=True)
    with open(bed_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        for c, s, e in bed_records:
            w.writerow([c, int(s), int(e)])
