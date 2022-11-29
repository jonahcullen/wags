# rescuer

Placeholder (description of the pipeline)

## Dependencies

- [Python](https://www.python.org/)
- [Mamba](https://github.com/mamba-org/mamba) or [Conda](https://conda.io/)
- [Snakemake](https://snakemake.readthedocs.io/)
- Miscellaneous python modules [pyaml](https://pyyaml.org/), [wget](https://bitbucket.org/techtonik/python-wget/), and [xlsxwriter](https://xlsxwriter.readthedocs.io/)
- [Apptainer/Singularity](https://apptainer.org/)
- [MinIO Client](https://min.io/docs/minio/linux/reference/minio-mc.html)

The `wags` Singularity container includes the following tools:
- [Placeholder]

## Initial setup

**1. Install dependencies**

```
# download Mambaforge installer (assuming Unix-like platform)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# updata Mamba
mamba update mamba -c conda-forge

# create a Snakemake environment which includes all Snakemake dependencies in
# addition to the miscellaneous modules above
mamba create \
    -c conda-forge -c bioconda \
    -n snakemake \ # name of the environment
    snakemake pyaml wget xlsxwriter
```

**2. Download the container**

Due to the size of the included reference genomes and index files (and depending on your internet speed) this should check ~5 minutes.

```
wget https://s3.msi.umn.edu/wags/wags.sif
```

**3. Clone this repo**

```
git clone git@github.com:jonahcullen/rescuer.git
```

**4. Install and setup MinIO client (optional)**

If desired, `wags` can work with MinIO Client (`mc`) to save pipeline outputs (e.g. FASTQ, BAM, GVCF), job-runner logs, and individual processing logs.

```
placeholder
```

## Reference genome

placeholder

## Prepare pipeline submissions

The executable `prep_subs.py` is responsible for generating all required pipeline inputs to be submitted to the cluster scheduler. The required input is a CSV file containing sample name, breed, sex (if known), and the associated FASTQ file prefix. The prefix should be the portion of the paired FASTQ name that distinguishes a given sample across multiple paired FASTQs (e.g. a sample split across multiple lanes). For example, if the input CSV was

```
dogid,breed,gender,fastq_id
SampleA,poodle,M,SampleA_2022
```

with FASTQs located in `/path/to/raw/data/`

```
SampleA_2022_L001_R1.fastq.gz
SampleA_2022_L001_R2.fastq.gz
SampleA_2022_L002_R1.fastq.gz
SampleA_2022_L002_R2.fastq.gz
SampleA_2000_L001_R1.fastq.gz
SampleA_2000_L001_R2.fastq.gz
```

`SampleA` would be matched with the two sets of paired FASTQs, one from lane L001 `SampleA_2022_L001_R1.fastq.gz/SampleA_2022_L001_R2.fastq.gz` and the other from L002 `SampleA_2022_L002_R1.fastq.gz/SampleA_2022_L002_R2.fastq.gz`. If however, `SampleA` should also have included the pair `SampleA_2000_L001_R1.fastq.gz/SampleA_2000_L001_R2.fastq.gz`, the input CSV should be modified as

```
dogid,breed,gender,fastq_id
SampleA,poodle,M,SampleA_20
```

The input CSV may contain as many samples as desired where `prep_subs.py` will generate pipeline inputs for each included sample. A minimal example of running `prep_subs.py` for processing using S3 storage is below (for all required and optional arguments `prep_subs.py --help`)

```
python ./prep_subs.py \
    --meta input.csv \
    --fastqs PATH/TO/FASTQ_DIR \ 
    --ref REF_GENOME \          
    --out PATH/TO/OUTDIR/ \                    
    --bucket RESULTS \                   
    --snake-env SNAKEMAKE_ENV \
    --partition PAR1,PAR2 \
    --email USER@email.com \
    --account SLURM_ACCT \
    --remote S3 \
    --alias MINIO_ALIAS
```

This will generate the following directory structure

```
PATH/TO/FASTQ_DIR/
    poodle/
        SampleA/
            REF_GENOME/
                input.tsv
                REF_GENOME_config.yaml
                poodle_SampleA.one_wags.slurm
                one_wags.smk
                rules/
                slurm_logs/
                slurm.go_wags/
                src/
```
