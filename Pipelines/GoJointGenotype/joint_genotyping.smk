import os
from datetime import datetime
from collections import defaultdict
import pandas as pd
import numpy as np

configfile: "canfam4_config.yaml"

# get dog ids from current vcf (FOR NOW FROM FILE BUT SHOULD AUTOMATE BY 
# BY USING BCFTOOLS QUERY -L)
d = defaultdict(list)

for root, dirs, files in os.walk(config["gvcf_dir"]):
    for f in files:
        if f.endswith(".g.vcf.gz"):
            # get dogid, breed, and absolute path to gvcf
            dogid = f.split(".")[0]
            breed = root.split("/")[-2]
           #gvcf = os.path.join(root,f)
            # add to dict
            d[dogid].append(breed)

units = pd.DataFrame.from_dict(d, orient="index", columns=["breed"])
units["dogid"] = units.index

# NOTE: CHANGE THIS TO GENERATE INTERVALS AS PART OF WORKFLOW
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
       #f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
       #f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_summary_metrics",
        # step09 - all variants table for af analysis
       #all_vars_table = os.path.join(
       #        "results/var_to_table/",
       #        datetime.today().strftime('%Y%m%d'),
       #        f"all.{datetime.today().strftime('%Y%m%d')}.table"
       #)
        expand("results/var_to_table/all.{date}.{ref}.table", 
                date=config["date"],ref=config["ref"])
       #read_interval_file(),
       #expand("results/{u.breed}/{u.dogid}/{u.dogid}_proc.txt",
       #        u=units.itertuples()
       #),
        # step01 - download_gvcf
       #expand("results/{u.breed}/{u.dogid}/{u.dogid}.{ref}.g.vcf.{ext}",
       #        ref=config["ref"], ext=["gz","gz.tbi"],
       #        u=units.itertuples()
       #)

#include: "rules/input_and_import.smk"
#include: "rules/genotype.smk"
#include: "rules/filter.smk"
#include: "rules/gather.smk"
#include: "rules/recal.smk"
#include: "rules/metrics.smk"
#include: "rules/vep.smk"
#include: "rules/table.smk"

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
        gvcfs = expand("{gvcf_dir}/{u.breed}/{u.dogid}/{u.dogid}.{ref}.g.vcf.gz",
                ref=config["ref"], gvcf_dir=config["gvcf_dir"],
                u=units.itertuples()
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
        interval  = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/CanFam4/Intervals/{interval}.interval_list"
    output:
        directory("results/import_gvcfs/{interval}")
    params:
        tmp_dir = "/dev/shm/{interval}",
        gatk    = config["gatk"],
        batch   = config["batch_size"],
    threads: 6
    resources:
         time   = 4800,
         mem_mb = 60000
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
        interval = "/panfs/roc/groups/0/fried255/shared/gatk4_workflow/GoDawgs/CanFam4/Intervals/{interval}.interval_list"
    output:
        vcf = "results/genotype_gvcfs/{interval}/output.vcf.gz",
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"]
    threads: 6
    resources:
         time   = 4800,
         mem_mb = 60000
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
        vcf = "results/genotype_gvcfs/{interval}/output.vcf.gz"
    output:
        var_filtrd_vcf = "results/fltr_make_sites_only/{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        sites_only_vcf = "results/fltr_make_sites_only/{interval}/filtr.{interval}.sites_only.variant_filtered.vcf.gz"
    params:
        gatk       = config["gatk"],
        excess_het = config["excess_het_threshold"],
        ref_fasta = config["ref_fasta"]
    resources:
         time   = 30,
         mem_mb = 12000
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
        gatk   = config["gatk"],
        picard = config["picard"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
    run:
        vcfs = " --input ".join(map(str,input.sites_only_vcf))

        shell(f'''
            set -e


            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                GatherVcfsCloud \
                --ignore-safety-checks \
                --gather-type BLOCK \
                --input {vcfs} \
                --output results/sites_only_gather_vcf/tmp.vcf.gz

            java -jar {{params.picard}} \
                SortVcf \
                I=results/sites_only_gather_vcf/tmp.vcf.gz \
                O={{output.gather_sites_only_vcf}}

            {{params.gatk}} --java-options "-Xmx18g -Xms6g" \
                IndexFeatureFile \
                --input {{output.gather_sites_only_vcf}}

            rm -f results/sites_only_gather_vcf/tmp.vcf.gz
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
        dbsnp_indels_vcf     = config["dbsnp_indels_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 24000
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
                --resource:dbsnp,known=false,training=true,truth=true,prior=10 {params.dbsnp_indels_vcf}
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
        dbsnp_snp_vcf           = config["dbsnp_snp_vcf"],
        broad_snp_vcf           = config["broad_snp_vcf"],
        axelsson_snp_vcf        = config["axelsson_snp_vcf"],
        illumina_snp_vcf        = config["illumina_snp_vcf"]
    threads: 4
    resources:
         time   = 240,
         mem_mb = 16000
    run:
        tranche_values = " -tranche ".join(map(str,params.recal_tranche_values))
        an_values = " -an ".join(map(str,params.recal_annotation_values))

        shell(f'''
            {{params.gatk}} --java-options "-Xmx12g -Xms3g" \
                VariantRecalibrator \
                -V {{input.gather_sites_only_vcf}} \
                -O {{output.snps_recal}} \
                --tranches-file {{output.snps_tranches}} \
                --trust-all-polymorphic \
                -tranche {tranche_values} \
                -an {an_values} \
                -mode SNP \
                --max-gaussians 6 \
                --resource:illumina,known=true,training=true,truth=true,prior=15.0 {params.illumina_snp_vcf} \
                --resource:broad,known=true,training=true,truth=true,prior=10.0 {params.broad_snp_vcf} \
                --resource:axelsson,known=false,training=true,truth=true,prior=8.0 {params.axelsson_snp_vcf} \
                --resource:dbSNP146,known=true,training=true,truth=true,prior=12.0 {params.dbsnp_snp_vcf}
        ''')

rule apply_recal:
    input:
        input_vcf       = "results/fltr_make_sites_only/{interval}/filtr.{interval}.variant_filtered.vcf.gz",
        indels_recal    = "results/recal/indels/indels.recal",
        indels_tranches = "results/recal/indels/indels.tranches",
        snps_recal      = "results/recal/snps/snps.recal",
        snps_tranches   = "results/recal/snps/snps.tranches"
    output:
        recal_vcf       = "results/apply_recal/{interval}/recal.{interval}.vcf.gz",
        recal_vcf_index = "results/apply_recal/{interval}/recal.{interval}.vcf.gz.tbi"
    params:
        gatk               = config["gatk"],
        indel_filter_level = config["indel_filter_level"],
        snp_filter_level   = config["snp_filter_level"]
    resources:
         time   = 30,
         mem_mb = 16000
    shell:
        '''
            set -e

            mkdir -p results/apply_recal/{wildcards.interval}/

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -O results/apply_recal/{wildcards.interval}/tmp.indel.recalibrated.vcf \
                -V {input.input_vcf} \
                --recal-file {input.indels_recal} \
                --tranches-file {input.indels_tranches} \
                --truth-sensitivity-filter-level {params.indel_filter_level} \
                --create-output-variant-index true \
                -mode INDEL

            {params.gatk} --java-options "-Xmx15g -Xms5g" \
                ApplyVQSR \
                -O {output.recal_vcf} \
                -V results/apply_recal/{wildcards.interval}/tmp.indel.recalibrated.vcf \
                --recal-file {input.snps_recal} \
                --tranches-file {input.snps_tranches} \
                --truth-sensitivity-filter-level {params.snp_filter_level} \
                --create-output-variant-index true \
                -mode SNP

            rm -f results/apply_recal/{wildcards.interval}/tmp.indel.recalibrated.vcf
        '''

rule final_gather_vcfs:
    input:
        recal_vcfs = sorted(expand("results/apply_recal/{interval}/recal.{interval}.vcf.gz", interval=intervals))
    output:
        final_vcf       = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
        final_vcf_index = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz.tbi"
    params:
        gatk = config["gatk"],
    threads: 4
    resources:
         time   = 240,
         mem_mb = 22000
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
        dbsnp_snp_vcf      = config["dbsnp_snp_vcf"],
        ref_dict           = config["ref_dict"],
        eval_interval_list = config["eval_interval_list"],
        metrics_prefix     = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}"
    threads: 8
    resources:
         time   = 360,
         mem_mb = 22000
    shell:
        '''
            set -e

            mkdir -p results/collect_metrics_on_vcf/

            {params.gatk} --java-options "-Xmx18g -Xms6g" \
                CollectVariantCallingMetrics \
                --INPUT {input.final_vcf} \
                --DBSNP {params.dbsnp_snp_vcf} \
                --SEQUENCE_DICTIONARY {params.ref_dict} \
                --OUTPUT {params.metrics_prefix} \
                --THREAD_COUNT 8 \
                --TARGET_INTERVALS {params.eval_interval_list}
        '''

rule vep_final_vcf:
    input:
        final_vcf = f"results/final_gather_vcfs/joint_genotype.{config['ref']}.vcf.gz",
       #detail_metrics  = f"results/collect_metrics_on_vcf/joint_genotype.{config['ref']}.variant_calling_detail_metrics",
    output:
        vep_vcf = f"results/vep_final_vcf/joint_genotype.{config['ref']}.vep.vcf.gz"
    params:
        conda_vep = config["conda_vep"],
        ref_fasta = config["ref_fasta"],
        ref_gtf   = config["ref_gtf"]
    threads: 12
    resources:
         time   = 4320,
         mem_mb = 60000
    run:
        import os
        out_name = os.path.splitext(output.vep_vcf)[0] 

        shell(f'''
            set +eu

            eval "$(conda shell.bash hook)"
            conda activate {params.conda_vep}

            vep \
                -i {{input.final_vcf}} \
                -o {out_name} \
                --gtf {{params.ref_gtf}} \
                --fasta {{params.ref_fasta}} \
                --everything \
                --force_overwrite \
                --vcf \
                --dont_skip

            bgzip --threads 12 -c {out_name} > {{output.vep_vcf}} &&
            tabix -p vcf {{output.vep_vcf}}
        ''')

#            vep \
#                -i {{input.final_vcf}} \
#                --gtf {{params.ref_gtf}} \
#                --fasta {{params.ref_fasta}} \
#                --everything \
#                --force_overwrite \
#                -o {out_name} \
#                --vcf \
#                --no_stats \
#                --offline \
#                --dont_skip

rule bcftools_stats:
    input:
        vep_vcf = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
    output:
        all_stats = "results/stats/joint_genotype.{date}.{ref}.vchk",
        prefix    = "results/stats/joint_genotype.{date}.{ref}"
    params:
        ref_fasta = config["ref_fasta"],
    threads: 1
    resources:
         time   = 720,
         mem_mb = 6000
    shell:
        '''
            bcftools stats \
                -F {params.ref_fasta} \
                -s - {input.vep_vcf} \
                > {output.all_stats}

            plot-vcfstats \
                -p {output.prefix} \
                {output.all_stats}
        '''

rule all_var_to_table:
    input:
        vep_vcf        = "results/vep_final_vcf/joint_genotype.{ref}.vep.vcf.gz",
        detail_metrics = "results/collect_metrics_on_vcf/joint_genotype.{ref}.variant_calling_detail_metrics",
       #all_stats      = "results/stats/joint_genotype.{date}.{ref}.vchk"
    output:
#        all_vars_table = os.path.join(
#                "results/var_to_table/",
#                datetime.today().strftime('%Y%m%d'),
#                f"all.{datetime.today().strftime('%Y%m%d')}.table"
#        )
        all_vars_table = "results/var_to_table/all.{date}.{ref}.table"
    params:
        gatk      = config["gatk"],
        ref_fasta = config["ref_fasta"],
        table_dir = config["var_to_table_dir"]
    threads: 12
    resources:
        #time   = 2160,
         time   = 840,
         mem_mb = 12000
    shell:
        '''
            {params.gatk} \
                VariantsToTable \
                -R {params.ref_fasta} \
                -V {input.vep_vcf} \
                -F CHROM -F POS -F REF -F ALT -F FILTER -F AF -F HOM-REF -F HET -F HOM-VAR -F NO-CALL -F CSQ \
                --show-filtered \
                -O {output.all_vars_table}
        '''

