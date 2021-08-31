# my_blobtoolkit
install and usage of btk on ccv

## Additional install instructions
Be sure to add the following requirements via pip in the main blobtools directory, with btk_env activated
```
python -m pip install -r requirements.txt
```

## create dataset

first, activate btk_env
```
conda activate /gpfs/home/ibishop/data/ibishop/condas/btk_env
```

#### prepare files
```
mkdir -p ~/scratch/btk_tutorial/files
cd ~/scratch/btk_tutorial/files
cp ~/scratch/ref_euk.fa ~/scratch/btk_tutorial/files/ref_euk.fa
```

#### split fasta into many pieces
```
#!/bin/bash

#SBATCH --time=0:01:00
#SBATCH --mem=4G #default is 1 core with 2.8GB of memory
#SBATCH -n 1
#SBATCH --account=epscor-condo

# Specify a job name:
#SBATCH -J splitting_fasta_contigs

#----- End of slurm commands ----

#load modules
module load seqkit

REF="$1"
SPLIT_NUM="$2"

#split ref contigs into 500 pieces
cd ~/scratch/btk_tutorial/files
seqkit split2 -p $SPLIT_NUM $REF
```


#### blobtools create
```
/users/ibishop/data/ibishop/blobtoolkit/blobtools2/blobtools create \
    --fasta ~/scratch/btk_tutorial/files/ref_euk.fa \
    --meta ~/scratch/btk_tutorial/files/ref_euk.yaml \
    --taxid 49265 \
    --taxdump /users/ibishop/data/ibishop/blobtoolkit/taxdump \
    ~/scratch/btk_tutorial/datasets/ref_euk
```


## add files

first, blast hits:
```
#you may have to remove leading zeroes in sequence to make batch work
blastn -db nt \
       -query ref_euk.part_${ID}.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out ref_euk.ncbi.blastn.part_${ID}.out"
       
       blastn -db $NT_DIR/nt \
       -query $CONTIG_DIR/contig_split_${ID}.fa \
       -out $OUTPUT_DIR/fin${ID}.out \
       -evalue 1e-6 \
       -max_hsps 1\
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10
       
cat ref_euk.ncbi.blastn.part_*.out > ref_euk.blastn.out
```
