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
import pandas as pd
from minio import Minio
from itertools import chain
from datetime import datetime
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
@click.option('--alias', default='s3',
              help='minio client alias for friedlab bucket')
@click.option('--samples', default='',
              help='sample IDs (comma separated) or file with one ID per row')
@click.option('--outdir', default='./',
              help='directory to send meta file (default: ./)')
@click.option('--outfile', default='dog_ids.csv',
              help='meta file name (default: dog_ids.csv)')
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
                    f'{os.path.join(fastq_dir, breed, dogid, flow) + "/"}'
                )
                print(fq_copy, file=out)
    print('Done!')


@click.command()
@click.option('--samples', default='',
              help='sample ID and breed ("sample,breed") or file with one ID per row')
@click.option('--ref', default='UU_Cfam_GSD_1.0_ROSY',
              help='reference to check against (default: UU_Cfam_GSD_1.0_ROSY)')
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

    done_outs = ''  # added to stop referenced before assignment warning
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
            # f'multiqc_fastqc.txt',
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


@click.command()
@click.option('--samples', default='',
              help='sample ID and breed ("sample,breed") or file with one ID per row')
@click.option('--ref', default='UU_Cfam_GSD_1.0_ROSY',
              help='reference to check against (default: UU_Cfam_GSD_1.0_ROSY)')
@click.option('--outdir', default='./fetched_logs',
              help='directory to send logs and multiqc report (default: fetched_logs)')
def fetch_logs(samples, ref, outdir):
    """
    fetch slurm run logs and mean depth from multiqc text file
    """
    # get dogs and breeds into d
    d = defaultdict(dict)
    if os.path.exists(samples):
        with open(samples, 'r') as infile:
            for line in infile:
                d[line.strip()]['logs'] = []
                d[line.strip()]['depth'] = []
    else:
        fq_list = [samples]
        for i in fq_list:
            d[i]['logs'] = []
            d[i]['depth'] = []

    # list all object paths in bucket that begin with my-prefixname
    objects = list(
        s3client.list_objects('friedlab',
                              prefix='wgs/',
                              recursive=True)
    )
    # filter to include only those with ref
    ref_objects = list(filter(lambda x: ref in x.object_name, objects))

    # get all logs and depth files into d
    for i in ref_objects:
        breed, dogid = i.object_name.split('/')[1:3]
        if dogid in d:
            d[dogid]['breed'] = breed
            if i.object_name.endswith('one_wags.err'):
                d[dogid]['logs'].append(i.object_name)
            if 'multiqc_general_stats.txt' in os.path.basename(i.object_name):
                d[dogid]['depth'].append(i.object_name)

    # filter d exclude dogs with more than one run file (ie took multiple restarts)
    filt_d = {k: v for k, v in d.items() if len(v['logs']) == 1}

    # loop through filt_d and download/parse the slurm logs and qc file for depth
    depth_col = 'QualiMap_mqc-generalstats-qualimap-mean_coverage'
    for k, v in filt_d.items():
        print(k, v['breed'])
        fetched_dir = os.path.join(outdir, v['breed'], k)
        os.makedirs(fetched_dir, exist_ok=True)
        # slurm log
        slurm_err = os.path.join(fetched_dir, os.path.basename(v['logs'][0]))
        if not os.path.isfile(slurm_err):
            s3client.fget_object('friedlab', v['logs'][0], slurm_err)
        # parse for processing time
        log_times = []
        with open(slurm_err, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('['):
                    log_times.append(line)
            # time stamp format from err files
            fmt = '[%a %b %d %H:%M:%S %Y]'
            runtime = datetime.strptime(log_times[-1], fmt) - datetime.strptime(log_times[0], fmt)
            filt_d[k]['runtime'] = (runtime.days * 24 * 60) + (runtime.seconds / 60)
        # qualimap mean coverage log
        depth_txt = os.path.join(fetched_dir, 'multiqc_general_stats.txt')
        if not os.path.isfile(depth_txt):
            s3client.fget_object('friedlab', v['depth'][0], depth_txt)
        # parse for mean depth
        df = pd.read_csv(depth_txt, sep='\t')
        filt_d[k]['mean_depth'] = df[depth_col].dropna().to_list()[0]

    # write to file
    with open(os.path.join(outdir, 'fetched.tsv'), 'w') as out:
        print('breed\tsample\tmean_depth\trun_time', file=out)
        for k, v in filt_d.items():
            print(v['breed'], k, v['mean_depth'], v['runtime'], sep='\t', file=out)


@click.command()
@click.option('--alias', default='s3',
              help='minio client alias for friedlab bucket')
@click.option('--bucket', default='',
              help='bucket name')
@click.option('--samples', default='',
              help='sample ID and breed ("sample,breed") or file with one ID per row')
@click.option('--ref', default='UU_Cfam_GSD_1.0_ROSY',
              help='reference to check against (default: UU_Cfam_GSD_1.0_ROSY)')
@click.option('--outdir', default='./fetched_gvcfs',
              help='directory to send logs and multiqc report (default: fetched_gvcfs)')
@click.option('--outfile', default='gvcfs.list',
              help='list of gvcfs locations (default: gvcfs.list)')
def fetch_gvcfs(alias, bucket, samples, ref, outdir, outfile):
    """
    fetch gvcfs from samples ids, prepare slurm submissions to download, and
    output file containing path to downloaded gvcfs (following submission of
    the prepared slurm jobs)
    """
    # get dogs and breeds into d
    d = defaultdict(dict)
    if os.path.exists(samples):
        with open(samples, 'r') as infile:
            for line in infile:
                d[line.strip()]['gvcf'] = ''
    else:
        fq_list = [samples]
        for i in fq_list:
            d[i]['gvcf'] = ''

    # list all object paths in bucket that begin with my-prefixname
    objects = list(
        s3client.list_objects(bucket,
                              prefix='wgs/',
                              recursive=True)
    )
    # filter to include only those with ref
    ref_objects = list(filter(lambda x: ref in x.object_name, objects))
    # get absolute path of outdir
    out_dir = os.path.abspath(os.path.expanduser(outdir))

    # write slurm submissions to download all gvcfs
    slurm_dir = os.path.join(out_dir, 'jobs')
    os.makedirs(slurm_dir, exist_ok=True)
    gvcfs_dir = os.path.join(out_dir, 'gvcfs')
    os.makedirs(gvcfs_dir, exist_ok=True)
    print(f'Generating slurm submissions in {slurm_dir}')
    # get gvcfs and indices
    for i in ref_objects:
        if i.object_name.endswith('.g.vcf.gz'):
            breed, dogid = i.object_name.split('/')[1:3]
            if dogid in d:
                d[dogid]['breed'] = breed
                d[dogid]['gvcf'] = [i.object_name]
                d[dogid]['gvcf'].append(i.object_name.replace('g.vcf.gz', 'g.vcf.gz.tbi'))
                # keep gvcf downloads organized
                fetched_dir = os.path.join(gvcfs_dir, breed, dogid)
                os.makedirs(fetched_dir, exist_ok=True)

    # get all fastqs into list
    gvcfs = list(chain(*[v['gvcf'] for v in d.values()]))
    # prepare slurm submissions to download all gvcfs
    for ind, i in enumerate(np.array_split(gvcfs, math.ceil(len(d.keys()) / 2))):
        with open(os.path.join(slurm_dir, f'gvcfs_{str(ind).zfill(4)}.slurm'), 'w') as out:
            print(top.replace('JOB_NAME', f'gvcfs_{str(ind).zfill(4)}'), file=out)
            for j in i:
                # get breed and sample
                tmp = j.split('/')
                breed, dogid = tmp[-5:-3]
                gvcf_copy = (
                    f'mc cp {os.path.join(alias, bucket, j)} '
                    f'{os.path.join(gvcfs_dir, breed, dogid) + "/"}'
                )
                print(gvcf_copy, file=out)

    # write gvcfs.list for input to joint genotyping
    with open(os.path.join(outdir, outfile), 'w') as out:
        for k, v in d.items():
            tmp = v['gvcf'][0].split('/')
            breed, dogid = tmp[-5:-3]
            print(dogid, os.path.join(gvcfs_dir, breed, dogid, os.path.basename(tmp[-1])), sep='\t', file=out)
    print('Done!')


messages.add_command(meta_prep)
messages.add_command(all_done)
messages.add_command(fetch_logs)
messages.add_command(fetch_gvcfs)

if __name__ == '__main__':
    messages()
