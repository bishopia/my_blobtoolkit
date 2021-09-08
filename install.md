# my_blobtoolkit
install blobtoolkit, including conda environment, github repos and all required dependencies for implementation.

## Install instructions
#### Set up conda environment
Here I am using a local installation of miniconda rather than the anaconda module available via HPC. make sure that mamba is installed in your base environment
```
#create btk_env conda environment
mamba create --prefix /gpfs/home/ibishop/data/ibishop/condas/btk_env
```

#### Activate btk_env and install software
```
conda activate /gpfs/home/ibishop/data/ibishop/condas/btk_env
mamba install -c conda-forge -y python=3.6 docopt psutil pyyaml ujson tqdm nodejs=10 yq
mamba install -c bioconda -y pysam seqtk;
mamba install -c conda-forge -y geckodriver selenium pyvirtualdisplay;
mamba install -c conda-forge minimap2 samtools
```

#### Go find directory to place git repo in, clone repos for toolkit
```
cd ~/data/ibishop
mkdir -p ~/blobtoolkit
cd blobtoolkit
git clone https://github.com/blobtoolkit/blobtools2;
git clone https://github.com/blobtoolkit/viewer;
git clone https://github.com/blobtoolkit/specification;
git clone https://github.com/blobtoolkit/insdc-pipeline;
git clone https://github.com/blobtoolkit/blobtools2.git

cd blobtools2
python -m pip install -r requirements.txt #see github issue #31
python -m pip install fastjsonschema
cd ..

cd viewer;
npm install;
cd ..;
```

#### Fetch relevant databases

##### get NCBI taxonomy
```
cd ~/data/ibishop/blobtoolkit
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
```

##### get NCBI nt database
For this step, be sure to be logged into a VNC node. If you're connected via ssh, you will be timed out and the download will fail.
```
cd ~/data/ibishop/blobtoolkit

mkdir -p nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P nt/ && \
        for file in nt/*.tar.gz; \
            do tar xf $file -C nt && rm $file; \
        done
```


##### get Uniprot database
```
cd ~/data/ibishop/blobtoolkit

mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e accession'\t'accession.version'\t'taxid'\t'gi > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
```

##### get BUSCO and lineages of interest
```
cd ~/data/ibishop/blobtoolkit
mkdir -p busco
wget -q -O eukaryota_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz" && tar xf eukaryota_odb10.gz -C busco
wget -q -O stramenopiles_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/stramenopiles_odb10.2020-08-05.tar.gz" && tar xf stramenopiles_odb10.gz -C busco
```

##### Deactivate environment and log out
```
conda activate
exit
```
