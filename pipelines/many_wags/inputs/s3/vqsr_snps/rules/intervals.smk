import pandas as pd
import os
import glob

localrules: intergenic_bed
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
                gunzip -c {params.ref_gtf} > {output.ref_gtf}
            else
                cp {params.ref_gtf} {output.ref_gtf}
            fi

            # bedops to convert gtf to 
            awk '$3 == "gene"' {output.ref_gtf} | gtf2bed > {output.ref_bed}

            # sort bed by chromosome and start position
            sortBed -i {output.ref_bed} -g {output.genome} > {output.sort_bed}

            # merge adjacent intervals along the same chrom and strand
            mergeBed -i {output.sort_bed} > {output.merge_bed}

            # complement the merged intervals to generate intergenic intervals
            complementBed -i {output.merge_bed} -g {output.genome} > {output.inter_bed}
        '''

localrules: intergenic_midpoints
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
localrules: bed_to_interval_list
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
localrules: chrom_intervals
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
# intervals based on mapping gaps (MAPQ<1), variation deserts (positions with 
# no small variants within <15kb), and patched with intergenic midpoints
localrules: desert_and_patch
rule desert_and_patch:
    input:
        intergen_mids = rules.intergenic_midpoints.output.midp_bed
    output:
        anchors = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/anchors.deserts_patched.tsv",
        bed     = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/anchors.deserts_patched.bed"
    params:
        ref_dict    = config['ref_dict'],
        map_gaps    = config['map_gaps'],
        var_deserts = config['var_deserts']
    run:
        spacing = 5000000
        priority = {"gap": 0, "desert": 1, "midpoint": 2}

        # load and label sources
        # midpoints of map gaps (regions with MAPQ<1)
        map_gaps = pd.read_csv(params.map_gaps)
        map_gaps["anchor"] = (map_gaps["start"] + map_gaps["end"]) // 2
        map_gaps["source"] = "gap"
        # midpoints of variant deserts (no small variatns <15kb)
        deserts = pd.read_csv(params.var_deserts)
        deserts["anchor"] = deserts["start"] + (deserts["length"] // 2)
        deserts["source"] = "desert"
        # intergenic midpoints
        mids = pd.read_csv(input.intergen_mids, sep="\t", header=None, names=["chrom", "start", "end"])
        mids["anchor"] = mids["end"]
        mids["source"] = "midpoint"

        all_anchors = pd.concat([
            map_gaps[["chrom", "anchor", "source"]],
            deserts[["chrom", "anchor", "source"]],
            mids[["chrom", "anchor", "source"]],
        ], ignore_index=True)

        all_anchors = all_anchors.drop_duplicates(subset=["chrom", "anchor"])
        all_anchors["priority"] = all_anchors["source"].map(priority)
        all_anchors = all_anchors.sort_values(["chrom", "anchor", "priority"])

        selected = []

        for chrom, chrom_df in all_anchors.groupby("chrom"):
            chrom_df = chrom_df.sort_values(["anchor", "priority"])
            max_anchor = chrom_df["anchor"].max()

            current_pos = 0
            while current_pos <= max_anchor:
                # get all candidates >= current_pos
                candidates = chrom_df[chrom_df["anchor"] >= current_pos]
                if candidates.empty:
                    break

                # pick the best anchor - first by position, then by priority
                best_anchor = candidates.groupby("anchor").first().reset_index().sort_values(["anchor", "priority"]).iloc[0]
                selected.append(best_anchor)

                # move to next 5mb step from current_pos, not from the anchor
                current_pos += spacing

        # read genome dict
        chrom_lengths = read_genome_dict(params.ref_dict)
        # prepare for output
        selected_df = pd.DataFrame(selected)
        # sanitize
        selected_df = selected_df.copy()
        selected_df["chrom"]  = selected_df["chrom"].astype(str).str.strip()
        selected_df["anchor"] = pd.to_numeric(selected_df["anchor"], errors="raise")
        # write to outputs
        selected_df = selected_df.drop_duplicates(subset=["chrom", "anchor"])
        selected_df = selected_df.sort_values(["chrom", "anchor"])
        selected_df.to_csv(output.anchors, sep="\t", index=False, header=False)
        create_bed_from_anchors(selected_df, chrom_lengths, output.bed)

rule desert_bed_to_interval_list:
    input:
        rules.desert_and_patch.output.bed
    output:
        desert_filt = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/beds/desert_patched.filt.bed",
        ival_list   = "{bucket}/wgs/pipeline/{ref}/{date}/intervals/desert_patched.interval_list"
    params:
        ref_dict = config['ref_dict']
    shell:
        '''
            # only keep chroms (drop haplotypes and unassigned)
            grep -P '^chr[0-9XYM]+\t' {input} > {output.desert_filt}

            java -jar /opt/wags/src/picard.jar \
                BedToIntervalList \
                I={output.desert_filt} \
                O={output.ival_list} \
                SD={params.ref_dict}
        '''

def get_interval_type(wc):
    if "nruns" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.N50.interval_list"
    elif "chroms" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/acgt.chrom.interval_list"
    elif "intergenic" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/intergenic_midp.interval_list"
    elif "patched" in config['anchor_type']:
        return "{bucket}/wgs/pipeline/{ref}/{date}/intervals/desert_patched.interval_list"
    else:
        print("Check config file for correct anchor type: nruns, chroms, intergenic, or patched")

localrules: generate_intervals
checkpoint generate_intervals:
    input:
        ival_list = get_interval_type
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
        chrom_re = re.compile(r'^chr(?:[0-9]+|X|Y|M|MT)$', re.IGNORECASE)
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
                    # keep only chrN, chrX, Y, M/Mt
                    if not chrom_re.match(chrom):
                        continue
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

