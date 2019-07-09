rule qc:
    input:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-QC.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preQC.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}"
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_qc","mem")
    script: 
        "../scripts/step1_qc.R"

rule post_qc:
    input:
        f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-QC.rds"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR, wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        min = config["rules"]["seurat_postqc"]["params"]["min_genes"],
        max = config["rules"]["seurat_postqc"]["params"]["max_genes"],
        mit = config["rules"]["seurat_postqc"]["params"]["mit_pct"],
        ribo = config["rules"]["seurat_postqc"]["params"]["ribo_pct"],
        filter_out = config["rules"]["seurat_postqc"]["params"]["filter_out"],
        filter_threshold = config["rules"]["seurat_postqc"]["params"]["filter_threshold"],
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_postqc","mem")
    script:
        "../scripts/step2_postqc.R"
