**Post-filtering violin plot**

Violin plot representing the number of reads, number of genes, mitochondrial percentage and ribosomal percentage in each sample after the filtering.

Filter paremeters used were:

- **Min genes expressed**: {{ snakemake.config["rules"]["seurat_postqc"]["min_feat"] }}
- **Max genes expressed**: {{ snakemake.config["rules"]["seurat_postqc"]["max_feat"] }}
- **Min counts**: {{ snakemake.config["rules"]["seurat_postqc"]["min_count"] }}
- **Max counts**: {{ snakemake.config["rules"]["seurat_postqc"]["max_count"] }}
- **Max mitochondrial percentage**: {{ snakemake.config["rules"]["seurat_postqc"]["mit_pct"] }} %
- **Max ribosomal percentage**: {{ snakemake.config["rules"]["seurat_postqc"]["ribo_pct"] }} %
