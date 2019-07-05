rule qc:
    input:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-QC.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}.preQC.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}.preQC.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}"
    threads: get_resource("seurat","threads")
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat","mem")
    script: 
        "../scripts/step1_qc.R"

