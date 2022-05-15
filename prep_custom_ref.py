
import os
import yaml
import gzip
import shutil
import textwrap
import argparse

def main():
    # prepare outdir
    ref_home = os.path.join(outdir,species,ref)
    ref_res  = os.path.join(ref_home,"resources")
    ref_fa   = os.path.join(ref_home,os.path.basename(fasta))

    # generate custom reference and resource directories
    os.makedirs(ref_res,exist_ok=True)
    
    # custom known sites
    d_sites = {}
    if sites:
        with open(os.path.abspath(sites),'r') as f:
            for line in f:
                name,vcf = line.strip().split(',')
                assert os.path.isfile(vcf), f"{name} was not found, check path is correct"
                # copy known site to outdir (default ~/.wags/SPECIES/REF/resources)
                vcf_out = os.path.join(ref_res,os.path.basename(vcf))
                shutil.copy(vcf,vcf_out)
                # add updated locatio to dictionary
                d_sites[name] = vcf_out
   
    config = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "Pipelines/GoProcessWGS/configs/custom_config.yaml"
    )

    # copy reference fasta to outdir
    shutil.copy(fasta,ref_home)

    # modify config file
    with open(config) as f:
        doc = yaml.safe_load(f)
    # update sif and other cli args
    doc['ref']       = ref
    doc['species']   = species
    doc['ref_dict']  = f"{os.path.splitext(ref_fa)[0]}.dict"
    # adjust ref_dict if fasta gzipped
   #if ref_fa.endswith(".gz"):
        
    doc['ref_fasta'] = ref_fa

    # add known site resources if provided
    doc['known_sites'] = d_sites

    # dump
    with open(os.path.join(ref_home,f"{ref}_config.yaml"),'w') as out:
        yaml.dump(doc,out,sort_keys=False)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog="wags (custom ref)",
        add_help=False,
        description=(
            "      ___           ___           ___           ___                    \n"
            "     /__/\         /  /\         /  /\         /  /\ *                 \n"
            "    _\_ \:\       /  /::\       /  /:/_       /  /:/_                  \n"
            "   /__/\ \:\     /  /:/\:\     /  /:/ /\     /  /:/ /\                 \n"
            "  _\_ \:\ \:\   /  /:/~/::\   /  /:/_/::\   /  /:/ /::\                \n"
            " /__/\ \:\ \:\ /__/:/ /:/\:\ /__/:/__\/\:\ /__/:/ /:/\:\               \n"
            " \  \:\ \:\/:/ \  \:\/:/__\/ \  \:\ /~~/:/ \  \:\/:/~/:/               \n"
            "  \  \:\ \::/   \  \::/       \  \:\  /:/   \  \::/ /:/                \n"
            "   \  \:\/:/     \  \:\        \  \:\/:/     \__\/ /:/                 \n"
            "    \  \::/       \  \:\        \  \::/        /__/:/   *prepare custom\n"
            "     \__\/         \__\/         \__\/         \__\/     reference     \n\n"
            "Prepare custom reference input for generating GVCFs from FASTQs by \n"
            "providing the reference FASTA.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument(
        "-r", "--ref",
        nargs="?",
        metavar="",
        help="reference name",
    )
    required.add_argument(
        "-s", "--species",
        nargs="?",
        metavar="",
        help="species name",
    )
    required.add_argument(
        "-f", "--fasta",
        nargs="?",
        metavar="",
        help="path to fasta to be used with --ref custom"
    )
    optional.add_argument(
        "-o", "--out",
        default=os.path.join(os.path.expanduser("~"),".wags/"),
        help="path to custom reference out dir [default: ~/.wags/SPECIES/REF]"
    )
    optional.add_argument(
        "--sites",
        help=textwrap.dedent('''\
            comma-separated file containing names (col 1) and
            paths to resource VCFs (and indices) (col 2) to be 
            used with --ref custom and --bqsr
        ''')
    )
    optional.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit"
    )
    
    args    = parser.parse_args()
    ref     = args.ref.lower()
    species = args.species
    fasta   = os.path.expanduser(args.fasta) \
        if "~" in args.fasta else os.path.abspath(args.fasta)
    outdir  = args.out
    sites   = args.sites

    # assert fasta exists
    assert os.path.isfile(fasta), "fasta not found, check path is correct"

    # gunzip if necessary
    if fasta.endswith(".gz"):
        uncomp_fasta = os.path.join(outdir,os.path.splitext(fasta)[0])
        if not os.path.isfile(uncomp_fasta):
            with gzip.open(fasta,'r') as f_in, open(uncomp_fasta,'wb') as f_out:
                shutil.copyfileobj(f_in,f_out)
        fasta = uncomp_fasta 
    # if non default outdir
    if outdir != os.path.join(os.path.expanduser("~"),".wags/"):
        if "~" in outdir:
            outdir = os.path.expanduser(outdir)
        else:
            outdir = os.path.abspath(outdir)
 
    main()
