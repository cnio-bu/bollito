curl --parallel --parallel-immediate $(cat files_reads) -o data/reads
curl --parallel --parallel-immediate $(cat files_ref) -o data/ref

for f in data/ref/*
do
    bunzip2 $f
done

mkdir data/ref/genome.chr19.fa_idx
STAR --runMode genomeGenerate --genomeDir data/ref/genome.chr19.fa_idx --genomeFastaFiles data/ref/genome.chr19.fa --genomeSAindexNbases 11

mkdir -p res/whitelists
curl -L "https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz?raw=true" | gunzip -c | cut -f1 > res/whitelists/3M-february-2018.txt
