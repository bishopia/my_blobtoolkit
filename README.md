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


## Adding files

#### First, blast hits:
```
#!/bin/bash

#SBATCH --time=0:30:00
#SBATCH --mem=4G #default is 1 core with 2.8GB of memory
#SBATCH -n 1
#SBATCH --account=epscor-condo
#SBATCH --array=1-500
##SBATCH --array=1-2

ID=$SLURM_ARRAY_TASK_ID
SEQUENCE=$(printf "%03d" ${ID})

# Specify a job name:
#SBATCH -J blastn_ref_euk

#----- End of slurm commands ----

#load modules
module load blast/2.9.0+

#run blastn
blastn -db ~/data/ibishop/blobtoolkit/nt/nt \
	-query ~/scratch/btk_tutorial/files/ref_euk.fa.split/ref_euk.part_${SEQUENCE}.fa \
	-outfmt "6 qseqid staxids bitscore std" \
	-max_target_seqs 10 \
	-max_hsps 1 \
	-evalue 1e-25 \
	-num_threads 1 \
	-out ~/scratch/btk_tutorial/files/out_split/ref_euk.ncbi.blastn.part_${SEQUENCE}.out
```

#### Then combine split results
```
cd $WORKDIR/files
cat out_split/*.out > blastn_nt_resuls.out
```

#### Now do diamond blast
```
#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem=8G #default is 1 core with 2.8GB of memory
#SBATCH -n 16
#SBATCH --account=epscor-condo
##SBATCH --array=1-500
#SBATCH --array=1-2

ID=$SLURM_ARRAY_TASK_ID
SEQUENCE=$(printf "%03d" ${ID})

# Specify a job name:
#SBATCH -J diamond_ref_euk

#----- End of slurm commands ----

#load modules
module load diamond

#run blastn
#blastn -db ~/data/ibishop/blobtoolkit/nt/nt \
#	-query ~/scratch/btk_tutorial/files/ref_euk.fa.split/ref_euk.part_${SEQUENCE}.fa \
#	-outfmt "6 qseqid staxids bitscore std" \
#	-max_target_seqs 10 \
#	-max_hsps 1 \
#	-evalue 1e-25 \
#	-num_threads 1 \
#	-out ~/scratch/btk_tutorial/files/out_split/ref_euk.ncbi.blastn.part_${SEQUENCE}.out
    
diamond blastx \
        --query ~/scratch/btk_tutorial/files/ref_euk.fa.split/ref_euk.part_${SEQUENCE}.fa \
        #--db /path/to/uniprot.db.with.taxids \
        --db ~/data/ibishop/blobtoolkit/uniprot \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 16 \
        > ~/scratch/btk_tutorial/ref_euk.diamond.blastx.part_${SEQUENCE}.out
```











