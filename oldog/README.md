Dogpile
=======

Usage
-----
```
$ ./Dogpile --help

usage: Dogpile -m -f -o [-i] [-h]

Dogpile generates all required input to process FASTQs to gVCF following GATK
best practices. For each sample (and associated FASTQ pair), Dogpile outputs a
directory structure organized by breed, wherein GATK pipeline input are
contained by sample ID.

required arguments:
  -m, --meta    csv of meta data and user ID to UMN ID conversions
  -f, --fastqs  list of full path of all fastqs
  -o, --out     path to out dir

optional arguments:
  -i, --ids         file containing list of UMN dog IDs to process
  -h, --help        show this help message and exit
```
