
rule input_list:
    output:
        gvcf_list = "{bucket}/wgs/pipeline/{ref}/{date}/import_gvcfs/inputs.list"
    params:
        cohort = config['joint_cohort']
    run:
        df = pd.read_csv(params.cohort,sep='\t')
        df.to_csv(
            output.gvcf_list,
            columns=['sample','gvcf'],
            sep='\t',
            header=False,index=False
        ) 
