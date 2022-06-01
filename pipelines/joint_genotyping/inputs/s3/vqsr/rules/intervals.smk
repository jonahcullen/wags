
import pandas as pd
import os
import glob

rule scatter_intervals:
    output:
        "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list",
    params:
        ref_fasta = config['ref_fasta']
    shell:
        '''
            java -jar /opt/wags/src/picard.jar \
                ScatterIntervalsByNs \
                R={params.ref_fasta} \
                OT=ACGT \
                N=50 \
                O={output}

            # dump chrun
            sed -i '/^chrUn/d' {output}
        '''

checkpoint generate_intervals:
    input:
        ival_list = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list"
    output:
        directory("{bucket}/wgs/pipeline/{ref}/{date}/intervals/import"),
    params:
        base      = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/import",
        lengths   = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/collapsed_lengths.csv",
        ivals_all = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/wags_intervals.list",
    run:
        # process interval file and collapse intervals by chromosome    
        d = {}
        header = []

        with open(input.ival_list, "r") as f:
            for line in f:
                # NOTE: canfam specific - make option?
                if "chrUn" in line:
                    continue
                
                if line.startswith("@"):
                    header.append(line)
                
                line = line.strip().split()
                if line[0] == "@SQ":
                    
                    # get chrom from header
                    chrom = line[1].split(":")[1]
                    d[chrom] = {}
                    # get chrom length from header
                    length = int(line[2].split(":")[1])
                    d[chrom]["len"] = length
                    
                
                elif line[0] in d.keys():
                    ival_start = int(line[1])
                    ival_end = int(line[2])
                    dist = ival_end - ival_start
                    
                    if not "intervals" in d[line[0]].keys():
                        d[line[0]]["intervals"] = [[ival_start, ival_end, dist]]
                    else:
                        d[line[0]]["intervals"].append([ival_start, ival_end, dist])

            # generate collapsed intervals
            for k,v in d.items():
                d[k]["collapsed"] = collapse_intervals(v["intervals"])

        # get list of keys and apply natural sort
        sorted_keys = list(d.keys())
        natural_sort(sorted_keys)

        ival = 0
        for k in sorted_keys:
            for i,j in enumerate(d[k]["collapsed"]):
                with open(
                        os.path.join(
                            params.base,
                            f"wags_{str(ival).zfill(4)}.interval_list"),
                            "w"
                        ) as out:
                    ival += 1
                    print(*header,sep="",end="",file=out)
                    print(k,
                          j[0],j[1],
                          "+",".",
                          sep="\t",file=out)
         
        # generate file of import interval paths
        with open(params.ivals_all,"w") as out:
            for ival in glob.glob(os.path.join(params.base,"*.interval_list")):
                print(ival,file=out)
        
        # generate collapsed lengths
        with open(params.lengths, "w") as out:
            print("chr","interval_number","interval_length",sep=",",file=out)
            for k,v in d.items():
                for i,j in enumerate(v["collapsed"]):
                    print(k,f"interval_{i}",j[1]-j[0],sep=",",file=out)

