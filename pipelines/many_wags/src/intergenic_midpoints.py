import sys
from collections import defaultdict

# argv[1] == ref_dict
# argv[2] == intergenic bed
# argv[3] == intergenic midpoints bed

# get ref dictionary
ref_d = {}
with open(sys.argv[1], 'r') as in_f:
    next(in_f)
    for line in in_f:
        line = line.strip().split('\t')
        ref_d[line[1].split(':')[-1]] = line[2].split(':')[-1]

# read complement.bed
d = defaultdict(list)
with open(sys.argv[2], 'r') as in_f:
    for line in in_f:
        line = line.strip().split()
        d[line[0]].append([line[1],line[2]])

# modify intervals to use the midpoints between intergenic regions
d_mod = defaultdict(list)
for k,v in d.items():
    midp = 0
    for ind,val in enumerate(v):
        if len(v) == 1:
            d_mod[k].append(val)
        elif ind == 0: # first interval
            start = int(val[0])
            end = int( (int(v[ind+1][0]) + int(v[ind+1][1]))/2 )
            d_mod[k].append([start,end])
           #midp += end + 1
            midp += end
        elif ind < len(v) - 1: # middle intervals
            start = midp
            end = int( (int(v[ind+1][0]) + int(v[ind+1][1]))/2 )
            d_mod[k].append([start, end])
           #midp = end + 1
            midp = end
        else: # last interval
            start = midp
            end = int(val[1])
            # modify end if final interval does not end at chrom length
            if end != int(ref_d[k]):
                end = int(ref_d[k])
            d_mod[k].append([start, end])

# write midpoint intergenic intervals to file
with open(sys.argv[3], 'w') as out_f:
    for k,v in d_mod.items():
        for i in v:
            print(k,i[0],i[1],sep='\t',file=out_f)

