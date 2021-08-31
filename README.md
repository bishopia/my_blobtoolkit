# my_blobtoolkit
install and usage of btk on ccv


### create dataset

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

#### blobtools create
```
/users/ibishop/data/ibishop/blobtoolkit/blobtools2 create \
    --fasta ~/scratch/btk_tutorial/files/ref_euk.fa \
    --meta ~/scratch/btk_tutorial/files/ref_euk.yaml \
    --taxid 49265 \
    --taxdump /users/ibishop/data/ibishop/blobtoolkit/taxdump \
    ~/scratch/btk_tutorial/datasets/ref_euk
```


