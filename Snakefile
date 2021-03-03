import os
from datetime import datetime
from collections import defaultdict
import pandas as pd
import numpy as np

configfile: "config.yaml"

# get dog ids from current vcf (FOR NOW FROM FILE BUT SHOULD AUTOMATE BY 
# BY USING BCFTOOLS QUERY -L)
d = defaultdict(list)

for root, dirs, files in os.walk(config["gvcf_dir"]):
    for f in files:
        if f.endswith(".g.vcf.gz"):
            # get dogid, breed, and absolute path to gvcf
            dogid = f.split(".")[0]
            breed = root.split("/")[1]
            gvcf = os.path.join(root,f)
            # add to dict
            d[dogid].append(breed)

units = pd.DataFrame.from_dict(d, orient="index", columns=["breed"])
units["dogid"] = units.index

# NOTE: CHANGE THIS TO GENERATE INTERVALS AS PART OF WORKFLWO
# get intervals
intervals, = glob_wildcards(os.path.join(config["db_intervals"],"{db_interval}.interval_list"))

rule all:
    input:
        # step00 - copy db import intervals
       #expand("results/import_gvcfs/{interval}",
       #        interval=intervals
       #)
        # step01 - import gvcfs
       #"data/inputs.list",
       #expand("results/import_gvcfs/{interval}", interval=intervals),
        # step02 - genotype gvcfs
       #expand("results/genotype_gvcfs/{interval}/output.vcf.gz", interval=intervals),
        # step03 - hard filter and make sites only
       #expand("results/fltr_make_sites_only/{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz",
       #        interval=intervals
       #),
        # step04 - gather sites only vcf
       #expand("results/sites_only_gather_vcf/gather.sites_only.vcf.gz",
       #        interval=intervals),
        # step05a - recal indels
       #"results/recal/indels/indels.recal",
       #"results/recal/indels/indels.tranches"
        # step06 - apply recalibration
       #expand("results/apply_recal/recal.{interval}.vcf.{ext}",
       #        interval=intervals, ext=["gz","gz.tbi"])
        # step07 - final gather vcfs
       #f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz"
       #f"results/final_gather_vcfs/joint_genotype.{datetime.today().strftime('%Y%m%d')}.vcf.gz"
        # step08a - vep final vcf
       #vep_vcf = f"results/vep_final_vcf/joint_genotype.{config['ref']}.vcf.gz"
        # step08b - collect metrics
        f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
        f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_summary_metrics",
        # step09 - all variants table for af analysis
        all_vars_table = os.path.join(
                config["var_to_table_dir"],
                datetime.today().strftime('%Y%m%d'),
                f"all.{datetime.today().strftime('%Y%m%d')}.table"
        )
       #read_interval_file(),
       #expand("results/{u.breed}/{u.dogid}/{u.dogid}_proc.txt",
       #        u=units.itertuples()
       #),
        # step01 - download_gvcf
       #expand("results/{u.breed}/{u.dogid}/{u.dogid}.{ref}.g.vcf.{ext}",
       #        ref=config["ref"], ext=["gz","gz.tbi"],
       #        u=units.itertuples()
       #)

# THIS RULE ONLY NEEDED TO COPY OVER PREVIOUSLY GENERATED INTERVAL DBS
#rule copy_intervals:
#    output:
#        interval_db = directory("results/import_gvcfs/{interval}")
#    resources:
#        time   = 60,
#        mem_mb = 6000, 
#        cpus   = 2
#    run:
#        import os, shutil
#
#        base = "/scratch.global/friedlab_cf3_TEST/cromwell-executions/JointGenotyping/7f320d18-c357-41ca-9ec1-7a7c4f32cab3/call-ImportGVCFs"
#        for root, dirs, files in os.walk(base):
#            if wildcards.interval in dirs:
#                dirs[:] = [wildcards.interval]
#                
#                src = os.path.join(root,dirs[0])
#               #dst = f"results/import_gvcfs/{wildcards.interval}"
#               #shutil.copytree(src, dst)
#                shell(f'''
#                    cp -r {src} {{output.interval_db}}
#                ''')

rule input_list:
    input:
        gvcfs = expand("gvcf/{u.breed}/{u.dogid}/{u.dogid}.{ref}.g.vcf.gz",
                ref=config["ref"], u=units.itertuples()
        )
    output:
        gvcf_list = "data/inputs.list"
    run:
        with open(output.gvcf_list, "w") as f:
            for i in input.gvcfs:
                dog = os.path.basename(i).split(".")[0]
                f.write(dog + "\t" + i + "\n") 

rule import_gvcfs:
    input:
        gvcf_list = "data/inputs.list",
        interval  = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/Intervals/{interval}.interval_list"
       #f"{os.path.join(config['db_intervals'],)}"
    output:
        directory("results/import_gvcfs/{interval}")
    params:
        tmp_dir = "/dev/shm/{interval}",
        gatk    = config["gatk"],
        batch   = config["batch_size"],
    resources:
         time   = 2160,
         mem_mb = 60000, 
         cpus   = 5
    shell:
        ''' 
            mkdir -p {params.tmp_dir}

            {params.gatk} --java-options "-Xms50g -Xmx50g" \
            GenomicsDBImport \
                --genomicsdb-workspace-path {output} \
                --batch-size {params.batch} \
                -L {input.interval} \
                --sample-name-map {input.gvcf_list} \
                --tmp-dir {params.tmp_dir} \
                --genomicsdb-shared-posixfs-optimizations true \
                --reader-threads 5 \
                -ip 500
        '''

rule genotype_gvcfs:
    input:
        ival_db  = "results/import_gvcfs/{interval}",
        interval = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/Intervals/{interval}.interval_list"
    output:
        vcf = "results/genotype_gvcfs/{interval}/output.vcf.gz",
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"]
    resources:
         time   = 2880,
         mem_mb = 60000, 
         cpus   = 4
    shell:
        '''
            {params.gatk} --java-options "-Xmx50g -Xms50g" \
            GenotypeGVCFs \
                -R {params.ref_fasta} \
                -O {output.vcf} \
                -G StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V gendb://{input.ival_db} \
                -L {input.interval}
        '''

rule fltr_make_sites_only:
    input:
       #ival_db  = "results/import_gvcfs/{interval}",
       #interval = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/Intervals/{interval}.interval_list",
        vcf      = "results/genotype_gvcfs/{interval}/output.vcf.gz"
    output:
        var_filtrd_vcf = "results/fltr_make_sites_only/{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        sites_only_vcf = "results/fltr_make_sites_only/{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz"
    params:
        gatk       = config["gatk"],
        excess_het = config["excess_het_threshold"],
        ref_fasta = config["ref_fasta"]
    resources:
         time   = 30,
         mem_mb = 12000, 
         cpus   = 1
    shell:
        '''
            {params.gatk} --java-options "-Xmx3g -Xms3g" \
            VariantFiltration \
                --filter-expression "ExcessHet > {params.excess_het}" \
                --filter-name ExcessHet \
                -O {output.var_filtrd_vcf} \
                -V {input.vcf}

            {params.gatk} --java-options "-Xmx3g -Xms3g" \
            MakeSitesOnlyVcf \
                --INPUT {output.var_filtrd_vcf} \
                --OUTPUT {output.sites_only_vcf}
        '''

rule sites_only_gather_vcf:
    input:
        sites_only_vcf = expand("results/fltr_make_sites_only/{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz",
                interval=intervals
        )
    output:
        gather_sites_only_vcf = "results/sites_only_gather_vcf/gather.sites_only.vcf.gz",
        gather_sites_only_tbi = "results/sites_only_gather_vcf/gather.sites_only.vcf.gz.tbi"
    params:
        gatk = config["gatk"],
    resources:
         time   = 240,
         mem_mb = 24000, 
         cpus   = 4
    run:
        vcfs = " --input ".join(map(str,input.sites_only_vcf))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {vcfs} \
                --output {{output.gather_sites_only_vcf}}

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {{output.gather_sites_only_vcf}}
        ''')

rule indels_var_recal:
    input:
        gather_sites_only_vcf = "results/sites_only_gather_vcf/gather.sites_only.vcf.gz",
    output:
        indels_recal    = "results/recal/indels/indels.recal",
        indels_tranches = "results/recal/indels/indels.tranches"
    params:
        gatk                    = config["gatk"],
        recal_tranche_values    = config["indel_recalibration_tranche_values"],
        recal_annotation_values = config["indel_recalibration_annotation_values"],
        dbsnp146_indels_vcf     = config["dbsnp146_indels_vcf"]
    resources:
         time   = 240,
         mem_mb = 24000, 
         cpus   = 4
    run:
        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
        an_values = " -an ".join(map(str,params.recal_annotation_values))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx24g -Xms24g" \
                VariantRecalibrator \
                -V {{input.gather_sites_only_vcf}} \
                -O {{output.indels_recal}} \
                --tranches-file {{output.indels_tranches}} \
                --trust-all-polymorphic \
                -tranche {tranche_values} \
                -an {an_values} \
                -mode INDEL \
                --max-gaussians 4 \
                --resource:dbsnp,known=false,training=true,truth=true,prior=10 {params.dbsnp146_indels_vcf}
        ''')


rule snps_var_recal:
    input:
        gather_sites_only_vcf = "results/sites_only_gather_vcf/gather.sites_only.vcf.gz",
    output:
        snps_recal    = "results/recal/snps/snps.recal",
        snps_tranches = "results/recal/snps/snps.tranches"
    params:
        gatk                    = config["gatk"],
        recal_tranche_values    = config["snp_recalibration_tranche_values"],
        recal_annotation_values = config["snp_recalibration_annotation_values"],
        dbsnp146_snp_vcf        = config["dbsnp146_snp_vcf"],
        broad_snp_vcf           = config["broad_snp_vcf"],
        axelsson_snp_vcf        = config["axelsson_snp_vcf"],
        illumina_snp_vcf        = config["illumina_snp_vcf"]
    resources:
         time   = 240,
         mem_mb = 16000, 
         cpus   = 4
    run:
        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
        an_values = " -an ".join(map(str,params.recal_annotation_values))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx12g -Xms3g" \
                VariantRecalibrator \
                    -V {{input.gather_sites_only_vcf}} \
                    -O {{output.snps_recal}} \
                    --tranches-file {{output.snp_tranches}} \
                    --trust-all-polymorphic \
                    -tranche {tranche_values} \
                    -an {an_values} \
                    -mode SNP \
                    --max-gaussians 6 \
                    --resource:illumina,known=true,training=true,truth=true,prior=15.0 {params.illumina_snp_vcf} \
                    --resource:broad,known=true,training=true,truth=true,prior=10.0 {params.broad_snp_vcf} \
                    --resource:axelsson,known=false,training=true,truth=true,prior=8.0 {params.axelsson_snp_vcf} \
                    --resource:dbSNP146,known=true,training=true,truth=true,prior=12.0 {params.dbsnp146_snp_vcf}
        ''')

rule apply_recal:
    input:
        input_vcf       = "results/fltr_make_sites_only/{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        indels_recal    = "results/recal/indels/indels.recal",
        indels_tranches = "results/recal/indels/indels.tranches",
        snps_recal      = "results/recal/snps/snps.recal",
        snps_tranches   = "results/recal/snps/snps.tranches"
    output:
        recal_vcf       = "results/apply_recal/recal.{interval}.vcf.gz",
        recal_vcf_index = "results/apply_recal/recal.{interval}.vcf.gz.tbi"
    params:
        gatk               = config["gatk"],
        indel_filter_level = config["indel_filter_level"],
        snp_filter_level   = config["snp_filter_level"]
    resources:
         time   = 1080,
         mem_mb = 16000, 
         cpus   = 1
    shell:
        '''
            set -e

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                    -O results/apply_recal/tmp.indel.recalibrated.vcf \
                    -V {input.input_vcf} \
                    --recal-file {input.indels_recal} \
                    --tranches-file {input.indels_tranches} \
                    --truth-sensitivity-filter-level {params.indel_filter_level} \
                    --create-output-variant-index true \
                    -mode INDEL

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                    -O {output.recal_vcf} \
                    -V results/apply_recal/tmp.indel.recalibrated.vcf \
                    --recal-file {input.snps_recal} \
                    --tranches-file {input.snps_tranches} \
                    --truth-sensitivity-filter-level {params.snp_filter_level} \
                    --create-output-variant-index true \
                    -mode SNP
        '''

rule final_gather_vcfs:
    input:
        recal_vcfs = expand("results/apply_recal/recal.{interval}.vcf.gz", interval=intervals)
    output:
        final_vcf       = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
        final_vcf_index = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz.tbi"
    params:
        gatk = config["gatk"],
    resources:
         time   = 240,
         mem_mb = 22000, 
         cpus   = 4
    run:
        vcfs = " --input ".join(map(str,input.recal_vcfs))

        shell(f'''
            set -e
            set -o pipefail

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {vcfs} \
                --output {{output.final_vcf}}

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {{output.final_vcf}}
        ''')

rule collect_metrics_on_vcf:
    input:
        final_vcf       = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
        final_vcf_index = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz.tbi"
    output:
        detail_metrics  = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
        summary_metrics = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_summary_metrics"
    params:
        gatk               = config["gatk"],
        dbsnp146_snp_vcf   = config["dbsnp146_snp_vcf"],
        ref_dict           = config["ref_dict"],
        eval_interval_list = config["eval_interval_list"]
    resources:
         time   = 360,
         mem_mb = 22000, 
         cpus   = 8
    run:
        metrics_prefix = f"joint_genotype.{config['ref']}"

        f'''
            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                CollectVariantCallingMetrics \
                    --INPUT {{input.final_vcf}} \
                    --DBSNP {{params.dbsnp146_snp_vcf}} \
                    --SEQUENCE_DICTIONARY {{params.ref_dict}} \
                    --OUTPUT {metrics_prefix} \
                    --THREAD_COUNT 8 \
                    --TARGET_INTERVALS {{params.eval_interval_list}}
        '''

rule vep_final_vcf:
    input:
        final_vcf = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
       #detail_metrics  = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
    output:
        vep_vcf = f"results/vep_final_vcf/joint_genotype.{config['ref']}.vcf.gz"
    params:
        conda_vep = config["conda_vep"]
    resources:
         time   = 1440,
         mem_mb = 60000, 
         cpus   = 12
    run:
        import os
        out_name = os.path.splitext(input.vep_vcf)[0]    

        shell(f'''
            set +eu

            eval "$(conda shell.bash hook)"
            conda activate {params.conda_vep}

            vep \
                -i {{input.final_vcf}} \
                --cache \
                --everything \
                -o {out_name} \
                --vcf \
                --no_stats \
                --species=canis_familiaris \
                --offline \
                --dont_skip

            bgzip -threads 12 -c {out_name} > {{output.vep_vcf}} &&
            tabix -p vcf {{output.vep_vcf}}
        ''')

rule all_var_to_table:
    input:
        vep_vcf = f"results/vep_final_vcf/joint_genotype.{config['ref']}.vcf.gz"
    output:
        all_vars_table = os.path.join(
                config["var_to_table_dir"],
                datetime.today().strftime('%Y%m%d'),
                f"all.{datetime.today().strftime('%Y%m%d')}.table"
        )
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
        table_dir = config["var_to_table_dir"]
    resources:
         time   = 2160,
         mem_mb = 12000, 
         cpus   = 12
    shell:
        '''
            {params.gatk} \
                VariantsToTable \
                    -R {params.ref_fasta} \
                    -V {input.vep_vcf} \
                    -F CHROM -F POS -F REF -F ALT -F FILTER -F AF -F HOM-REF -F HET -F HOM-VAR -F NO-CALL -F CSQ \
                    --showFiltered \
                    -O {output.all_vars_table}
        '''

#        bslmm_seeds = expand(
#                            "results/step02/{{breed}}/{{trait}}/{{metab}}/seeds/BSLMM_10M_seed{seed}.param.txt",
#                                        seed=[*range(10,60,5)]
#                                                )
