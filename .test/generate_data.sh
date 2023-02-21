git clone http://github.com/cnio-bu/sc-test-data data

cd data
snakemake --use-conda
cd ..
mkdir data/ref/genome.chr19.fa_idx
STAR --runMode genomeGenerate --genomeDir data/ref/genome.chr19.fa_idx --genomeFastaFiles data/ref/genome.chr19.fa --genomeSAindexNbases 11

mkdir -p res/whitelists
curl -L "https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz?raw=true" | gunzip -c | cut -f1 > res/whitelists/3M-february-2018.txt
