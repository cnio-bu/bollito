rule seurat:
    input:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/seurat_final.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}"
    threads: get_resource("seurat","threads")
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat","mem")
    script: 
        "../scripts/seurat.R"
