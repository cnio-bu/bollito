rule seurat_qc:
    input: 
        seurat_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds",
        before_filt_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/1_vlnplot_ngene_numi_pctmit_beforefilt.pdf", caption="../report/conf/before_filt_plot.rst", category="2_Single-cell QC"),
        gene_umi_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/3_geneplot_numi_vs_pctmit_ngene.pdf", caption="../report/conf/gene_umi_plot.rst", category="2_Single-cell QC")
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preQC.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        samples_path = config["samples"],
        units_path = config["units"],
        input_type = config["input_type"],
        random_seed = config["random_seed"],
        sample = f"{{sample}}",
        min_cells_per_gene = config["rules"]["seurat_qc"]["params"]["min_cells_per_gene"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_qc","mem"),
        walltime=get_resource("seurat_qc","walltime")
    script: 
        "../scripts/step1_qc.R"

rule seurat_post_qc:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds",
        after_filt_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/4_vlnplot_ngene_numi_pctmit_afterfilt.pdf", caption="../report/conf/after_filt_plot.rst", category="2_Single-cell QC"),
        stats_table=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/6_pre_vs_post_stats.tsv", caption="../report/conf/stats_table.rst", category="2_Single-cell QC")
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        min_feat = config["rules"]["seurat_postqc"]["params"]["min_feat"],
        max_feat = config["rules"]["seurat_postqc"]["params"]["max_feat"],
        min_count = config["rules"]["seurat_postqc"]["params"]["min_count"],
        max_count = config["rules"]["seurat_postqc"]["params"]["max_count"],
        mit = config["rules"]["seurat_postqc"]["params"]["mit_pct"],
        ribo = config["rules"]["seurat_postqc"]["params"]["ribo_pct"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_postqc","mem"),
        walltime=get_resource("seurat_postqc","walltime")
    script:
        "../scripts/step2_postqc.R"

rule seurat_filter:
    input:
        lambda wc: f"{OUTDIR}/seurat/{wc.sample}/1_preprocessing/seurat_post-qc.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc-filtered.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.filter.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.filter.bmk"
    params: 
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        gene = config["rules"]["seurat_filter"]["params"]["gene"],
        filter_out = config["rules"]["seurat_filter"]["params"]["filter_out"],
        threshold = config["rules"]["seurat_filter"]["params"]["threshold"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_filter", "mem"),
        walltime=get_resource("seurat_filter", "walltime")
    script:
        "../scripts/step2.1_filter.R"

def norm_input(wc):
    if wc.sample == "integrated":
        return ""
    else:
        if config["rules"]["seurat_filter"]["params"]["gene"]:
            return f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc-filtered.rds"
        else:
            return f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds"

rule seurat_normalization:
    input:
        seurat_obj=norm_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds",
        elbow_plot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/3_elbowplot.pdf", caption="../report/conf/elbow_plot.rst", category="3_Normalization"),
        dimplot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/2_dimplot.pdf", caption="../report/conf/dimplot.rst", category="3_Normalization"),
        cellcycle_plot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/6_cell_cycle_dimplot.pdf", caption="../report/conf/cellcycle_plot.rst", category="3_Normalization")
    log:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        normalization = config["rules"]["seurat_normalization"]["params"]["normalization"],
        regress_out = config["rules"]["seurat_normalization"]["params"]["regress_out"],
        vars_to_regress = config["rules"]["seurat_normalization"]["params"]["vars_to_regress"],
        regress_cell_cycle = config["rules"]["seurat_normalization"]["params"]["regress_cell_cycle"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_normalization","mem"),
        walltime=get_resource("seurat_normalization","walltime")
    script:
        "../scripts/step3_normalization.R"

rule seurat_integration:
    input:
        data=get_integration_input_sm
    output:
        seurat_obj=f"{OUTDIR}/seurat/integrated/2_normalization/seurat_normalized-pcs.rds"
    log:
        f"{LOGDIR}/seurat/integrated/2_normalization/integrated.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/integrated/2_normalization/integrated.normalization.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/integrated",
        random_seed = config["random_seed"],
        norm_type = config["rules"]["seurat_integration"]["params"]["norm_type"],
        vars_to_regress = config["rules"]["seurat_integration"]["params"]["vars_to_regress"]    
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_integration","mem"),
        walltime=get_resource("seurat_integration","walltime")
    script:
        "../scripts/step3.2_integration.R"

rule seurat_find_clusters:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        pc = config["rules"]["seurat_find_clusters"]["params"]["principal_components"],
        res = config["rules"]["seurat_find_clusters"]["params"]["resolutions"],
        k_neighbors = config["rules"]["seurat_find_clusters"]["params"]["k_neighbors"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_find_clusters","mem"),
        walltime=get_resource("seurat_find_clusters","walltime")
    script:
        "../scripts/step4_find-clusters.R"

rule seurat_degs:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/4_degs/seurat_degs.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/4_degs/{{sample}}.seurat_degs.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/4_degs/{{sample}}.seurat_degs.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["rules"]["seurat_degs"]["params"]["selected_res"],
        test = config["rules"]["seurat_degs"]["params"]["test"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_degs","mem"),
        walltime=get_resource("seurat_degs","walltime")
    script:
        "../scripts/step5_degs.R"

rule seurat_gs:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/5_gs/seurat_complete.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/5_gs/{{sample}}.seurat_complete.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/5_gs/{{sample}}.seurat_complete.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        gs_collection = config["rules"]["seurat_gs"]["params"]["geneset_collection"]
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("seurat_gs","mem"),
        walltime=get_resource("seurat_gs","walltime")
    script:
        "../scripts/step6_gs-scoring.R"

rule slingshot:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        seurat_obj=f"{OUTDIR}/slingshot/{{sample}}/6_traj_in/slingshot_sce_objects.RData"
    log:
        f"{LOGDIR}/slingshot/{{sample}}/6_traj_in/{{sample}}.slingshot.log"
    benchmark:
        f"{LOGDIR}/slingshot/{{sample}}/6_traj_in/{{sample}}.slingshot.bmk"
    params:
        output_dir = f"{OUTDIR}/slingshot/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["rules"]["slingshot"]["params"]["selected_res"],
        start_clus = config["rules"]["slingshot"]["params"]["start_clus"],
        end_clus = config["rules"]["slingshot"]["params"]["end_clus"],
        n_var_genes = config["rules"]["slingshot"]["params"]["n_var_genes"],
        n_plotted_genes = config["rules"]["slingshot"]["params"]["n_plotted_genes"],
        pc = config["rules"]["seurat_find_clusters"]["params"]["principal_components"]
    conda: "../envs/slingshot.yaml"
    resources:
        mem=get_resource("slingshot","mem"),
        walltime=get_resource("slingshot","walltime")
    script:
        "../scripts/step7_traj_in.R"

rule vision:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        seurat_obj=f"{OUTDIR}/vision/{{sample}}/7_func_analysis/vision_object.rds"
    log:
        f"{LOGDIR}/vision/{{sample}}/7_func_analysis/{{sample}}.vision.log"
    benchmark:
        f"{LOGDIR}/vision/{{sample}}/7_func_analysis/{{sample}}.vision.bmk"
    params:
        output_dir = f"{OUTDIR}/vision/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["rules"]["vision"]["params"]["selected_res"],
        mol_signatures = config["rules"]["vision"]["params"]["mol_signatures"],
        meta_columns = config["rules"]["vision"]["params"]["meta_columns"],
        n_cores = config["rules"]["vision"]["params"]["n_cores"]
    conda: "../envs/vision.yaml"
    resources:
        mem=get_resource("vision","mem"),
        walltime=get_resource("vision","walltime")
    script:
        "../scripts/step8_func_analysis.R"
