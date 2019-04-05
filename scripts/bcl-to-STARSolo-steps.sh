#! /bin/env bash

# los genomas de referencia están en: /opt/bioinfo/cellranger/references/3.0.0/
# los datos están en /home/shared/single_cell
# Cell ranger: /opt/bioinfo/cellranger/software/3.0.2/
# Muestras JL: /home/cfustero/Projects/20190320_single-cell-test-run/HNWNLBGX9/jl_High

# Step 1: If using bcl files. Cell ranger mkfastq
/opt/bioinfo/cellranger/software/3.0.2/cellranger mkfastq --run /home/shared/single_cell/data/cnio_test_run/run_folder/ --csv /home/shared/single_cell/data/cnio_test_run/samplesheet.csv --output-dir /home/cfustero/Projects/20190320_single-cell-test-run;

# Step 2: Combine fastqs from different lanes into one
cat jl_High_S7_L00?_R1_001.fastq.gz > jl_High_S7_R1.fastq.gz;
cat jl_High_S7_L00?_R2_001.fastq.gz > jl_High_S7_R2.fastq.gz;

# Step 3: Create mm10 genome Reference for STAR 2.7.0e (run just once.)
sbatch -o scJL.log -e scJL.err -J Index-mm10 -c 6 --mem=32G --wrap "STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /home/cfustero/NGS/STARSolo-mm10 --genomeFastaFiles /opt/bioinfo/cellranger/references/3.0.0/refdata-cellranger-mm10-3.0.0/fasta/genome.fa --sjdbGTFfile /opt/bioinfo/cellranger/references/3.0.0/refdata-cellranger-mm10-3.0.0/genes/genes.gtf";

# Step 4: Run alignment: Gene. V3 chemistry: UMI length = 12nt, CB = 16nt
sbatch -o scJL.log -e scJL.err -J scJL -c 6 --mem=50G --wrap "STAR --soloType Droplet
--readFilesCommand gunzip -c --readFilesIn /home/cfustero/Projects/20190320_single-cell-test-run/HNWNLBGX9/jl_High/jl_High_S7_R2.fastq.gz /home/cfustero/Projects/20190320_single-cell-test-run/HNWNLBGX9/jl_High/jl_High_S7_R1.fastq.gz 
--genomeDir /home/cfustero/NGS/STARSolo-mm10 --soloCBwhitelist /home/cfustero/Projects/20190320_single-cell-test-run/3M-february-2018.txt 
--soloFeatures Gene --outFilterMultimapNmax 50 --winAnchorMultimapNmax 50 --alignEndsType EndToEnd 
--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --soloUMIlen 12" 