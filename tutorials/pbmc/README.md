# PBMC tutorial
## Tutorial overview

In this tutorial we will analyze the *PBMC* dataset. This dataset is composed of various types of *human* blood cells including T cells, B cells, NK cells, and monocytes.

The sample is part of [10x Genomics example datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets). In this case the dataset has been
downsampled to 1 million reads in order to reduce the execution time required to complete the tutorial.

The tutorial will analyse the dataset using the following steps:
- sequencing QC
- single-cell QC
- single-cell normalization
- clustering
- differential expression analysis
- functional analysis.

## Input file configuration

The first step is to download the data we will be analysing. In this case it will be composed of two pairs of paired-end, FASTQ-formatted files.

You should download these files and put then in a directory of your choice. We will use the example directory `/my/data/files` throughout the tutorial. You should instead use the real path to your files.

First pair of files, corresponding to lane 1 of the flowcell:
- [3M_pbmc_1k_v3_S1_L001_R1_001.fastq.gz](https://bioinformatics.cnio.es/data/pipelines/single_cell/10x_sample_data/pbmc_1M/3M_pbmc_1k_v3_S1_L001_R1_001.fastq.gz)
- [3M_pbmc_1k_v3_S1_L001_R2_001.fastq.gz](https://bioinformatics.cnio.es/data/pipelines/single_cell/10x_sample_data/pbmc_1M/3M_pbmc_1k_v3_S1_L001_R2_001.fastq.gz)

Second pair of files, corresponding to lane 2 of the flowcell:
- [3M_pbmc_1k_v3_S1_L002_R1_001.fastq.gz](https://bioinformatics.cnio.es/data/pipelines/single_cell/10x_sample_data/pbmc_1M/3M_pbmc_1k_v3_S1_L002_R1_001.fastq.gz)
- [3M_pbmc_1k_v3_S1_L002_R2_001.fastq.gz](https://bioinformatics.cnio.es/data/pipelines/single_cell/10x_sample_data/pbmc_1M/3M_pbmc_1k_v3_S1_L002_R2_001.fastq.gz)

> The files with *R1* (for *read 1*) in their names store the unique mollecular identifiers (UMI) and the cell barcodes.
> The files with *R2* (for *read 2*) in their names store the cDNA data.

### Samples file (samples.tsv)
We will now configure the samples file, *samples.tsv*, to adjust it to our data. Since we are dealing with a single sample, so we will need to enter a single line:

- In the *sample* column we will enter the sample name. We will use the name of the dataset, *PBMC*.

Your sample.tsv should look like this:

| sample |
| ------ |
| PBMC   |

> TIP: you can add arbitrary columns to the sample.tsv file. These will then be added
> to the output Seurat objects as metadata associated to the sample.

### Units file (units.tsv)
We will now configure the units file, *units.tsv*, to indicate which data files belong to each sample.

- In the *sample* column we will indicate the sample name (NOTE: this should match the name in *samples.tsv*)
- In the *unit* column we will indicate a unique name for each set of files that belong to the same sample.
Let's use *l1* and *l2* to indicate thay they come from different lanes.
- In the *fq1* column we will indicate the path to our file containing the UMIs and barcodes.
- In the *fq2* column we will indicate the path to our file containing the cDNA data.

Your units.tsv file should look like this:

| sample | unit | fq1 | fq2 |
| ------ | ---- | ---------- | ----- |
| PBMC   | l1   | /my/data/files/3M_pbmc_1k_v3_S1_L001_R1_001.fastq.gz      | /my/data/files/3M_pbmc_1k_v3_S1_L001_R2_001.fastq.gz |
| PBMC   | l2   | /my/data/files/3M_pbmc_1k_v3_S1_L002_R1_001.fastq.gz      | /my/data/files/3M_pbmc_1k_v3_S1_L002_R2_001.fastq.gz |

> NOTE: remember to change the file paths to wherever you stored the downloaded fastq files.

### Main configuration file (config.yaml)

We have now successfully configured our input data. Let's now configure resources and parameters
that the pipeline needs in order to correctly process the data.

The configuration file features a large amount of parameters. Out of these we will only be adjusting the ones
needed to correctly process our current dataset, and to obtain the results we are after.

#### Required files

Since the sample we are analysing is from Homo sapiens, you will need to download a copy of the human genome sequence, and
the corresponding annotation file:

- [Human genome sequence](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz)
- [Annotation file](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz)
- [Barcode whitelist file](https://bioinformatics.cnio.es/data/pipelines/single_cell/10x_sample_data/pbmc_1M/3M-february-2018.txt.gz)

> You should descompress these files and place them in a directory of your choice. Again, we will use */my/data/files* here.

#### General parameters

The following parameters affect the pipeline as a whole. Use the following table to set the correct values.

| Name               | Value | Notes                                             |
|--------------------|-------|--------------------------------------------------------|
| input_type         | fastq | Input files are FASTQ in fastq format.                                 |
| technology         | 10x   | Input files are generated with 10x Genomics technology. |
| case               | uppercase | Human genes are generally written using upper case letters.              |
| random_seed        | 4848  | A random seed. Using the same seed allows us to have replicable results between pipeline runs.                                  |
| annotation         | /my/data/files/gencode.v34.annotation.gtf | The GENCODE annotation file you just downloaded.  |
| fasta              | /my/data/files/GRCh38.primary_assembly.genome.fa | The GENCODE sequence file you just downloaded.         |
| idx                | /my/data/files/GRCh38.primary_assembly.genome_index  | The path where to store the genome index.                  |
| whitelist           | /my/data/files/3M-February-2018.txt | Barcode whitelist for 10x Chromium v3 chemistry. |

### Step-specific parameters

#### Alignment step (star)
STAR extra parameter is configured depending on the sample. In this case, due to the Chromium V3 chemistry, the string is:
"--soloType Droplet --soloFeatures Gene Velocyto --outFilterMultimapNmax 50 --winAnchorMultimapNmax 50 --alignEndsType EndToEnd --outReadsUnmapped Fastx --soloUMIlen 12 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM"

#### Single-cell pre-QC (seurat_qc)
For the purpose of this tutorial we want to maximise the number of genes and cell we use. We will therefore set the required number of
cells in which a gene must be detected to just 1.

| Parameter        | Value |  Commentary                                             |
|------------------|-------|---------------------------------------------------------|
| min_cells_per_gene | 1 | Require a gene to be detected in at least 1 cell for it to be included in the analysis |

#### Single-cell post-QC (seurat_postqc)
In this section we will establish some thresholds for the inclusion and exclusion of cells, based on results from the previous rules.

| Name             | Value |  Commentary                                             |
|------------------|-------|---------------------------------------------------------|
| min_feat         | 1000  | Bottom limit for the number of expressed genes.         |
| max_feat         | 4000  | Upper limit for the number of expressed genes.          |
| mit_pct          | 15    | Upper limit for the mitochondrial percentage of counts. |

#### Single-cell normalization (seurat_normalization)

The normalization method chosen was SCT. The QC parameters (features, counts, mitochondrial percent and ribosomal percent)
does not affect heavily the sample, so the regression was not needed. Similar conclusions were obtained by the cell cycle study.

| Name             | Value |  Commentary                                             |
|------------------|-------|---------------------------------------------------------|
| normalization         | SCT   | SCT method           |
| regress_out           | False | Upper limit for the number of expressed genes.          |
| vars_to_regress       | Empty | Unused filtering.                                       |
| regress_cell_cycle    | Empty | Unused filtering.                                       |
| regress_merge_effect  | Empty | Upper limit for the mitochondrial percentage of counts. |

#### Single-cell clustering - “seurat_clustering“
A low set of resolution were selected, since few clusters were wanted. k parameter was left by default.
As the result of the previous rule indicated, the significant PCA components were the first 7.

| Name                  | Value                   | Commentary                          |
|-----------------------|-------------------------|-------------------------------------|
| principal_components  | 7                       | Significant principal components.   |
| resolutions           | [0.2, 0.4, 0.6, 0.8, 1] | Set of tested resolutions.          |
| k_neighbors           | 20                      | Number of *k* neighbors.            |

#### Single-cell functional analysis Seurat-based - “seurat_gs“
In this step, it is required a GMT file storing the molecular signatures and a threshold value.
This value reflects the minimum ratio (expressed genes / total genes) for a geneset to be tested. In this case, some PBMC related signature were tested.

| Name             | Value |  Commentary                                             |
|----------------|-------|---------------------------------------------------------|
| geneset_collection  | PBMC_related_msigs.gmt   | Molecular signatures GMT file.    |
| geneset_percentage  | 0.2 |Minimum ratio (expressed genes / total genes) for a geneset to be tested |

#### Single-cell functional analysis Vision-based - “vision“
To obtain a wider interpretation of the sample the molecular signatures tested were the hallmarks from MSigDB.
Also, metadata variables from the Seurat object were added and the clustering resolution chosen is 0.4.

| Name             | Value |  Commentary                                             |
|------------------|-------|---------------------------------------------------------|
| perform         | True   | Set to "True" to execute the analysis.          |
| mol_signatures  | Hallmarks.gmt | Upper limit for the number of expressed genes.          |
| meta_colums     | ["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"] | Metadata variables from Seurat. |
| n_cores         | 8 | Threads provided to Vision.                                      |
| selected_res    | 0.4 | Clustering resolution chosen. |

#### Single-cell differential expression analysis - “seurat_degs“
The differential expression analysis only needs a statistical test, which, in this case, is the Wilcoxon test, and cluster distribution related to resolution.

| Name             | Value |  Commentary                                             |
|------------------|-------|---------------------------------------------------------|
| selected_res  | 0.4   | Clustering resolution chosen.  |
| test   | "wilcox" | Statistical test to use for the DE analysis   |

#### Resources per rule specification.

| Rule                        | RAM memory | Threads    |
|-----------------------------|------------|------------|
| star                        | 64000      | 8 |
| fastqc                      |            |  | 
| rseqc_junction_saturation   |            |           | 
| rseqc_readdup               |                        |           | 
| multiqc                     |               |   | 
| seurat_qc                   |                   |   |
| seurat_postqc               |                        |   |
| seurat_normalization        | 64000         |  |
| seurat_find_clusters        |               |  | 
| seurat_degs                 |              |           | 
| seurat_gs                   |                        |           | 
| vision                      |             |  | 

## Execute the pipeline
Once all the input files are ready, the pipeline is executed. The number of available cores is The command line used is:

snakemake --use-conda -j N

At this point, Snakemake decides the execution order based on the rules resources specification and triggers the first rules.

## Review the results

In this tutorial, the pipeline will be reviewed step by step, focusing on the main results and the parameter values that need to be specified in the configuration file.
