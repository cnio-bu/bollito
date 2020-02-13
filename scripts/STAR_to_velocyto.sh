#!/usr/bin/python
head -n 2 ${1}/matrix.mtx > ${1}/mtx_header.txt
head -n 3 ${1}/matrix.mtx | tail -n 1 > ${1}/mtx_summary.txt

tail -n +4 ${1}/matrix.mtx | cut -d " " -f 1-3 > ${1}/ms.txt
tail -n +4 ${1}/matrix.mtx | cut -d " " -f 1-2,4 > ${1}/mu.txt
tail -n +4 ${1}/matrix.mtx | cut -d " " -f 1-2,5 > ${1}/ma.txt

#rm -r ${1}/spliced ${1}/unspliced ${1}/ambiguous
mkdir -p ${1}/spliced ${1}/unspliced ${1}/ambiguous

cat ${1}/mtx_header.txt ${1}/mtx_summary.txt ${1}/ms.txt > ${1}/spliced/matrix.mtx
cat ${1}/mtx_header.txt ${1}/mtx_summary.txt ${1}/mu.txt > ${1}/unspliced/matrix.mtx
cat ${1}/mtx_header.txt ${1}/mtx_summary.txt ${1}/ma.txt > ${1}/ambiguous/matrix.mtx

cp ${1}/features.tsv ${1}/genes.tsv
cp ${1}/genes.tsv ${1}/barcodes.tsv ${1}/spliced/
cp ${1}/genes.tsv ${1}/barcodes.tsv ${1}/unspliced/
cp ${1}/genes.tsv ${1}/barcodes.tsv ${1}/ambiguous/

rm ${1}/*.txt

