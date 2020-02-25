**Post-filtering violin plot**

Violin plot representing the number of reads, number of genes and mitochondrial percentage in each sample after the filtering.

Filter paremeters used were:

- **Min genes expressed**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["min_feat"] }}
- **Max genes expressed**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["max_feat"] }}
- **Min counts**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["min_count"] }}
- **Max counts**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["max_count"] }}
- **Max mitochondrial percentage**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["mit_pct"] }} %
- **Max ribosomal percentage**: {{ snakemake.config["rules"]["seurat_postqc"]["params"]["ribo_pct"] }} %
