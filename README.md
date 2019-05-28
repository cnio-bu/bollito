# Bollito SC pipeline

## Test data

The system is pre-configured to run an example based on sample data available from 10X Genomics. The required datasets can be found at these URLS. Please update the "units.tsv" file to point at the data as needed.

* https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/neuron_10k_v3
* https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/heart_10k_v3

## FASTQ\_SCREEN

The pipeline can optionally run FASTQ\_SCREEN to check the samples for contamination.

To disable it use the config file option ```config["rules"]["fastq_screen"]["disabled"]```.

Config file pointing to indexes should be placed in a directory named ```config['rules']['fastq_screen_indexes']['outdir']/FastQ_Screen_Genomes```.

If the rule is enabled and no config file provided, default indexes will be downloaded using the command ```fastq_screen --get_genomes```.

