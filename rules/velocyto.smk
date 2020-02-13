rule STAR_to_velocyto:
    input:
        f"{OUTDIR}/star/{{sample}}/Solo.out/Velocyto/Summary.csv"
    output: 
        f"{OUTDIR}/star/{{sample}}/Solo.out/Velocyto/raw/spliced/matrix.mtx"
    log:
        f"{LOGDIR}/star/{{sample}}/Solo.out/Velocyto/raw/{{sample}}.STAR_to_velocyto.log"
    benchmark:
        f"{LOGDIR}/star/{{sample}}/Solo.out/Velocyto/raw/spliced/{{sample}}.STAR_to_velocyto.bmk"
    params:
        input_dir= lambda wc: "{}/star/{}/Solo.out/Velocyto/raw".format(OUTDIR,wc.sample)
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("STAR_to_velocyto","mem"),
        walltime=get_resource("STAR_to_velocyto","walltime")
    shell:"""
        bash scripts/STAR_to_velocyto.sh {params.input_dir} 2> {log}
    """

def get_velocyto_dirs(wc):
     paths = f"{OUTDIR}/star/{{sample}}/Solo.out/Velocyto/raw/"
     paths = [i for i in paths if not ('/integrated/' in i)]
     return paths


rule velocyto:
    input:
        data = f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output: 
        f"{OUTDIR}/velocyto/{{sample}}/8_RNA_velocity/seurat_velocity.rds"
    log:
        f"{LOGDIR}/velocyto/{{sample}}/8_RNA_velocity/{{sample}}.velocyto.log"
    benchmark:
        f"{LOGDIR}/velocyto/{{sample}}/8_RNA_velocity/{{sample}}.velocyto.bmk"
    params:
        velocyto_dir = f"{OUTDIR}/star/{{sample}}/Solo.out/Velocyto/raw/",
        output_dir = f"{OUTDIR}/velocyto/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["rules"]["velocyto"]["params"]["selected_res"] 
    conda: "../envs/velocyto.yaml"
    resources:
        mem=get_resource("velocyto","mem"),
        walltime=get_resource("velocyto","walltime")
    script:
        "../scripts/step9_RNA_velocity.R"
