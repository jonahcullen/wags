#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 11:39:16 2022

@author: jonahcullen
"""

import os
import math
import click
import textwrap
import numpy as np
from minio import Minio
from itertools import chain
from collections import defaultdict

# setup minio client
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')
if s3_key_id is None or s3_access_key is None:
    print('No access key is available.')

# initialize minioClient with an endpoint and access/secret keys.
s3client = Minio(
    's3.msi.umn.edu',
    access_key=s3_key_id,
    secret_key=s3_access_key,
    secure=True,
)


# https://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings
def common_prefix(strings):
    """
    Find the longest string that is a prefix of all the strings.
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix


# write to slurm submission files
top = textwrap.dedent("""\
    #!/bin/bash -l
    #SBATCH -t 24:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=6gb
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=cull0084@umn.edu
    #SBATCH --job-name JOB_NAME.slurm
    #SBATCH -o %j.JOB_NAME.out
    #SBATCH -e %j.JOB_NAME.err
    #SBATCH -p msismall,msilarge

    set -e

    cd $SLURM_SUBMIT_DIR                      
""")


@click.group()
def messages():
    pass


@click.command()
@click.option('--alias', default='s3', help='minio client alias for friedlab bucket')
@click.option('--samples', default='', help='sample IDs (comma separated) or file with one ID per row')
@click.option('--outdir', default='./', help='directory to send meta file')
@click.option('--outfile', default='dog_ids.csv', help='meta file name (default dog_ids.csv')
def meta_prep(alias, samples, outdir, outfile):
    """
    generate metadata from a sample ID (separated by commas) or
    from a text file with one ID per row.
    """

    # read in file, alias name - should use click to for input file
    d = defaultdict(dict)

    if os.path.exists(samples):
        with open(samples, 'r') as infile:
            for line in infile:
                d[line.strip()]['fastqs'] = []
    else:
        fq_list = [samples]
        for i in fq_list:
            d[i]['fastqs'] = []

    out_dir = os.path.realpath(os.path.expanduser(outdir))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # list all object paths in bucket that begin with my-prefixname.
    print('Finding all S3 objects...')
    objects = list(
        s3client.list_objects('friedlab',
                              prefix='wgs/',
                              recursive=True)
    )

    # and generate dict

    # get all fastqs associated with sample ids
    for i in objects:
        if not i.object_name.endswith('fastq.gz') or i.object_name.count('/') != 5 or 'broken' in i.object_name:
            continue
        breed, dogid = i.object_name.split('/')[1:3]
        if dogid in d.keys():
            d[dogid]['fastqs'].append(os.path.join(i.bucket_name, i.object_name))

    print(f'Writing meta file to {os.path.join(out_dir, outfile)}')
    # collapse fastq list to longest common prefix and write to meta
    with open(os.path.join(out_dir, outfile), 'w') as out:
        print('sample,breed,sex,fastq', file=out)
        for k, v in d.items():
            d[k]['breed'] = v['fastqs'][0].split('/')[2]
            d[k]['fq_prefix'] = common_prefix(list(map(os.path.basename, v['fastqs'])))
            if d[k]['fq_prefix'].endswith('_R'):  # clip the R of R1/R2 from fastq names
                d[k]['fq_prefix'] = d[k]['fq_prefix'][:-2]
            print(k, d[k]['breed'], '', d[k]['fq_prefix'], sep=',', file=out)

    # write slurm submissions to download all fastqs
    slurm_dir = os.path.join(out_dir, 'jobs')
    os.makedirs(slurm_dir, exist_ok=True)
    fastq_dir = os.path.join(out_dir, 'fastqs')
    os.makedirs(fastq_dir, exist_ok=True)
    print(f'Generating slurm submissions in {slurm_dir}')
    # get all fastqs into list
    fqs = list(chain(*[v['fastqs'] for v in d.values()]))
    for ind, i in enumerate(np.array_split(fqs, math.ceil(len(d.keys()) / 2))):
        with open(os.path.join(slurm_dir, f'fastqs_{str(ind).zfill(4)}.slurm'), 'w') as out:
            print(top.replace('JOB_NAME', f'fastqs_{str(ind).zfill(4)}'), file=out)
            for j in i:
                # get breed, id, and flowcell from j to ensure fastqs with same
                # on different runs to not get overwritten
                tmp = j.split("/")
                breed, dogid = tmp[2:4]
                flow = tmp[5]
                fq_copy = (
                    f'mc cp {os.path.join(alias, j)} '
                    f'{os.path.join(fastq_dir,breed,dogid,flow)}'
                )
                print(fq_copy, file=out)
    print('Done!')


@click.command()
@click.option('--samples', default='', help='sample ID and breed ("sample,breed") or file with one ID per row')
@click.option('--ref', default='UU_Cfam_GSD_1.0_ROSY', help='reference to check against (default: UU_Cfam_GSD_1.0_ROSY)')
def all_done(samples, ref):
    """
    check if all done dog
    """
    # get dogs and breeds into d
    d = {}

    if os.path.exists(samples):
        with open(samples, 'r') as infile:
            next(infile)
            for line in infile:
                sample, breed = line.strip().split(',')[:2]
                d[sample] = breed
    else:
        sample, breed = samples.strip().split(',')
        d[sample] = breed

    l_all_done = []
    d_not_done = {}

    for k, v in d.items():
        # expected files for an all done dog
        done_outs = [
            f'{k}.{ref}.cram',
            f'{k}.{ref}.cram.crai',
            f'{k}.{ref}.cram.md5',
            f'{k}.{ref}.g.vcf.gz',
            f'{k}.{ref}.g.vcf.gz.tbi',
            f'{k}.{ref}.analyze_cov.csv',
            f'{k}.{ref}.analyze_cov.pdf',
            f'{k}.{ref}.flagstat.txt',
            f'multiqc.log',
            f'multiqc_data.json',
           #f'multiqc_fastqc.txt',
            f'multiqc_general_stats.txt',
            f'multiqc_picard_dups.txt',
            f'multiqc_qualimap_bamqc_genome_results.txt',
            f'multiqc_samtools_flagstat.txt',
            f'multiqc_sources.txt',
            f'multiqc_report.html',
            f'genome_results.txt',
            f'{k}.delly.{ref}.vcf.gz',
            f'{k}.delly.{ref}.vcf.gz.tbi',
            f'{k}.gridss.{ref}.bam',
            f'{k}.gridss.{ref}.vcf.gz',
            f'{k}.gridss.{ref}.vcf.gz.tbi',
            f'{k}.lumpy.{ref}.vcf.gz',
            f'{k}.lumpy.{ref}.vcf.gz.tbi',
            f'{k}.manta.diploidSV.{ref}.vcf.gz',
            f'{k}.manta.diploidSV.{ref}.vcf.gz.tbi',
            'svCandidateGenerationStats.tsv',
            'svLocusGraphStats.tsv',
            'alignmentStatsSummary.txt'
        ]

        # list all object paths in bucket that begin with my-prefixname.
        objects = list(
            s3client.list_objects('friedlab',
                                  prefix=f'wgs/{v}/{k}/{ref}',
                                  recursive=True)
        )

        # get all (29) currently completed outputs
        # NOTE - 29 instead of 30 to ignore multiqc_fastqc.txt
        is_done = [
            os.path.basename(i.object_name)
            for i in objects
            if i.object_name.split('/')[4] in ['cram', 'gvcf', 'svar', 'qc']
            and 'multiqc_fastqc.txt' not in i.object_name
        ]

        # check dog is done
        if set(done_outs) == set(is_done):
            l_all_done.append(k)
        else:
            d_not_done[k] = list(set(done_outs) - set(is_done))

    # print to stdout
    print('----all done dogs----')
    print('\n'.join(l_all_done))
    if d_not_done.keys():
        print('--NOT all done dogs--')
        for k, v in d_not_done.items():
            print(k, f'({len(done_outs) - len(v)}/{len(done_outs)})')


messages.add_command(meta_prep)
messages.add_command(all_done)

if __name__ == '__main__':
    messages()

# either write something to download them in parallel or write a function to
# generate a suitable number of slurm submissions...
# s3client.fget_object("friedlab", "wgs/gldr/D06264/fastq/HLCNLDSXX/S_2007_S94_R1_001.fastq.gz",
#                      "./S_2007_S94_R1_001.fastq.gz")

# and i.obje¥ct_name.count("/") == 6

# this is not entirely safe...requires a bunch of hand annotation to fix fastq names...

# f = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/PerformEval/random_100.fastqs.list"
#
# d = defaultdict(list)
#
# with open(f, "r") as infile:
#     for line in infile:
#         line = line.strip().split("/")
#         breed, sample = line[3:5]
#         fq = re.search(r"(.*)_(R1|1|R2|2)_{0,1}(.*)\.(f.*q)(\..*){0,1}", line[-1]).group(1)
#         print(fq)
#         if "L00" in fq:
#             fq = fq[:-5]
#         key = f"{sample}_{breed}"
#         d[key].append(fq)
#
# meta = "/Users/jonahcullen/projects/friedenberg/gatk_pipeline/PerformEval/dog_ids.csv"
#
# with open(meta, "w") as out:
#     for k, v in d.items():
#         sample, breed = k.split("_")
#         print(sample, breed, "", v[0], sep=",", file=out)
