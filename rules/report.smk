rule scQC_report:
    input:
        data=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds"
    output:
        report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/1_vlnplot_ngene_numi_pctmit_beforefilt.pdf", caption="../report/conf/beforefilt.rst", category="Single-cell Quality Control"),
        report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/4_vlnplot_ngene_numi_pctmit_afterfilt.pdf", caption="../report/conf/afterfilt.rst",category="Single-cell Quality Control"),
        report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/3_geneplot_numi_vs_pctmit_ngene.pdf", caption="../report/conf/geneplot.rst", category="Quality control")

