import pandas as pd
import os
import glob

rule intergenic_bed:
    output:
        genome    = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/genome.txt",
        ref_gtf   = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/ref.gtf",
        ref_bed   = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/ref.bed",
        sort_bed  = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/sort.bed",
        merge_bed = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/merge.bed",
        inter_bed = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/intergenic.bed",
    params:
        ref_dict = config['ref_dict'],
        ref_gtf  = config['ref_gtf']
    shell:
        '''
            set -e

            # generate genome.txt file
            tail -n +2 {params.ref_dict} | cut -f 2,3 > {output.genome}
            sed -i 's/SN:\|LN://g' {output.genome}

            # check and uncompress if gtf gzipped
            if file {params.ref_gtf} | grep -q gzip; then
                echo HELLO
                gunzip -c {params.ref_gtf} > {output.ref_gtf}
            else
                cp {params.ref_gtf} {output.ref_gtf}
                echo NOWAY
            fi

            # bedops to convert gtf to 
            gtf2bed < {output.ref_gtf} > {output.ref_bed}

            # sort bed by chromosome and start position
            sortBed -i {output.ref_bed} -g {output.genome} > {output.sort_bed}

            # merge adjacent intervals along the same chrom and strand
            mergeBed -i {output.sort_bed} > {output.merge_bed}

            # complement the merged intervals to generate intergenic intervals
            complementBed -i {output.merge_bed} -g {output.genome} > {output.inter_bed}
        '''

rule intergenic_midpoints:
    input:
        inter_bed = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/intergenic.bed",
    output:
        midp_bed = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/intergenic_midp.bed"
    params:
        ref_dict = config['ref_dict']
    shell:
        '''
            python src/intergenic_midpoints.py \
                {params.ref_dict} \
                {input.inter_bed} \
                {output.midp_bed}
        '''

# NOTE - NEED TO MAKE THE DROPPING ON UNPLACED CONTIGS FLEXIBLE FOR OTHER
# SPECIES
rule bed_to_interval_list:
    input:
        "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/intergenic_midp.bed"
    output:
        midp_filt = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/intergenic_midp.filt.bed",
        ival_list = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/intergenic_midp.interval_list"
    params:
        ref_dict = config['ref_dict']
    shell:
        '''
            # for canine ref dump chrun
            # NOTE - needs to be flexible...
            sed '/^chrUn/d' {input} > {output.midp_filt}

            java -jar /opt/wags/src/picard.jar \
                BedToIntervalList \
                I={output.midp_filt} \
                O={output.ival_list} \
                SD={params.ref_dict}
        '''

# chromosome-length intervals - again excluding chrUn WHICH needs to be
# adjustable in the config
if "chroms" in config['anchor_type']:
    rule chrom_intervals:
        output:
            ival = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.chrom.interval_list",
        params:
            ref_dict = os.path.basename(config['ref_dict'])
        run:
            d = {}
            header = []
            with open(params.ref_dict, 'r') as f_in, open(output.ival, 'w') as f_out:
                for line in f_in:
                    if line.startswith("@"):
                        header.append(line)
                    line = line.strip().split()
                    if line[0] == "@SQ":
                        # get chrom from header
                        chrom = line[1].split(":")[1]
                        length = int(line[2].split(":")[1])
                        d[chrom] = length
                
                # write to file
                print(*header, sep="", end="", file=f_out)
                for k,v in d.items():
                    if 'chrUn' not in k:
                        print(k, '1', v, '+\tACGTmer', sep='\t', file=f_out)

# intervals based on nruns (default as does not require the reference genome
# to have a quality annotation file)
if "nruns" in config['anchor_type']:
    rule scatter_intervals:
        output:
            "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list",
        params:
            nrun_length = config['nrun_length'],
            ref_fasta   = config['ref_fasta']
        resources:
            time = 60,
            mem_mb = 20000        
        shell:
            '''
                java -jar /opt/wags/src/picard.jar \
                    ScatterIntervalsByNs \
                    R={params.ref_fasta} \
                    OT=ACGT \
                    N={params.nrun_length} \
                    O={output}

                # dump chrun
                sed -i '/^chrUn/d' {output}
            '''

def get_interval_type(wc):
    if "nruns" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list"
    elif "chroms" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.chrom.interval_list"
    elif "intergenic" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/intergenic_midp.interval_list"
    else:
        print("Check config file for correct anchor type: nruns, chroms, intergenic")

checkpoint generate_intervals:
    input:
        ival_list = get_interval_type
       #ival_list = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list"
       #    if "nruns" in config['anchor_type'] else "{bucket}/wgs/pipeline/{ref}/{date}/intervals/intergenic_midp.interval_list"
    output:
        directory("{bucket}/wgs/pipeline/{ref}/{date}/intervals/import"),
    params:
        ival_length = config['interval_length'],
        base        = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/import",
        lengths     = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/collapsed_lengths.csv",
        ivals_all   = f"{config['bucket']}/wgs/pipeline/{config['ref']}/{config['date']}/intervals/wags_intervals.list",
    run:
        os.makedirs(params.base,exist_ok=True)
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
                d[k]["collapsed"] = collapse_intervals(
                    v["intervals"], 
                    int(params.ival_length)
                )

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

