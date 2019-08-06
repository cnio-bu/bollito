rule seurat_qc:
    input:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preQC.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        project_name = config["rules"]["seurat_qc"]["params"]["project_name"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_qc","mem")
    script: 
        "../scripts/step1_qc.R"

rule seurat_post_qc:
    input:
        f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds"
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
rule seurat_normalization:
    input:
        f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR, wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_normalization","mem")
    script:
        "../scripts/step3_normalization.R"
rule seurat_find_clusters:
    input:
        f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR, wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        seed =  config["rules"]["seurat_find_clusters"]["params"]["random_seed"],
        pc = config["rules"]["seurat_find_clusters"]["params"]["principal_components"],
        res = config["rules"]["seurat_find_clusters"]["params"]["resolutions"],
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_find_clusters","mem")
    script:
        "../scripts/step4_find-clusters.R"
rule seurat_degs:
    input:
        f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/4_degs/seurat_degs.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/4_degs/{{sample}}.seurat_degs.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/4_degs/{{sample}}.seurat_degs.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR, wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        selected_res = config["rules"]["seurat_degs"]["params"]["selected_res"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_degs","mem")
    script:
        "../scripts/step5_degs.R"
rule seurat_gs:
    input:
	    f"{OUTDIR}/seurat/{{sample}}/4_degs/seurat_degs.rds"
    output:
        data=f"{OUTDIR}/seurat/{{sample}}/5_gs/seurat_complete.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/5_gs/{{sample}}.seurat_complete.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/5_gs/{{sample}}.seurat_complete.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR, wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        gs_collection = config["rules"]["seurat_gs"]["params"]["geneset_collection"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_gs","mem")
    script:
        "../scripts/step6_gs-scoring.R"

