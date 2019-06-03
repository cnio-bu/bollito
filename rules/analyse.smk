rule seurat:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    output:
        "{}/seurat/{{sample}}.done".format(OUTDIR)
    log:
        "{}/seurat/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/seurat/{{sample}}.bmk".format(LOGDIR)
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample)
    threads: get_resource("seurat","threads")
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat","mem")
    script: 
        "../scripts/seurat.R"
