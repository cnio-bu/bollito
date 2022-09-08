if [ -z "$1" ]
then
    PYTHON_VERSION=3.6
else
    PYTHON_VERSION=$1
fi

git clone http://gitlab.com/bu_cnio/sc-test-data data
conda create -y -q -n bollito_snakemake snakemake python=${PYTHON_VERSION}
source activate bollito_snakemake
cd data
snakemake --use-conda
cd -
conda create -y -q -n bollito_prep star=2.7.1a
conda activate bollito_prep
mkdir data/ref/genome.chr19.fa_idx
STAR --runMode genomeGenerate --genomeDir data/ref/genome.chr19.fa_idx --genomeFastaFiles data/ref/genome.chr19.fa --genomeSAindexNbases 11

mkdir -p res/whitelists
curl -L "https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz?raw=true" | gunzip -c | cut -f1 > res/whitelists/3M-february-2018.txt
