If you are using an original wags container (as opposed to cat_wags.sif) with NCBI scaffold IDs, you will need to adjust the config names in this directory

mv Fca126_mat1.0_config.yaml Fca126_mat1.0_config.NEW.yaml
mv Fca126_mat1.0_config.OLD.yaml Fca126_mat1.0_config.yaml

to ensure the original Fca126_mat1.0 reference is used

ref_fasta : /home/refgen/cat/Fca126_mat1.0/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna

However, if you would you prefer to use updated cat_wags with "nice" scaffold IDs (ie chrA1, chrA2, chrA3)

wget https://s3.msi.umn.edu/wags/cat_wags.sif
