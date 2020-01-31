rule STAR_to_velocyto:
    input:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    output: 
        f"{OUTDIR}/star/{{sample}}/Solo.out/Velocyto/raw/spliced/matrix.mtx"
    log:
        f"{LOGDIR}/star/{{sample}}/Solo.out/Velocyto/raw/spliced/{{sample}}.STAR_to_velocyto.log"
    benchmark:
        f"{LOGDIR}/star/{{sample}}/Solo.out/Velocyto/raw/spliced/{{sample}}.STAR_to_velocyto.bmk"
    params:
        input_dir= lambda wc: "{}/star/{}/Solo.out/Velocyto/raw".format(OUTDIR,wc.sample)
    conda: "../envs/seurat.yaml"
    resources:
        mem=get_resource("STAR_to_velocyto","mem"),
        walltime=get_resource("STAR_to_velocyto","walltime")
    shell:"""
        head -n 2 {params.input_dir}/matrix.mtx > {params.input_dir}/mtx_header.txt
        head -n 3 {params.input_dir}/matrix.mtx | tail -n 1 > {params.input_dir}/mtx_summary.txt

        tail -n +4 {params.input_dir}/matrix.mtx | cut -d " " -f 1-3 > {params.input_dir}/ms.txt
        tail -n +4 {params.input_dir}/matrix.mtx | cut -d " " -f 1-2,4 > {params.input_dir}/mu.txt
        tail -n +4 {params.input_dir}/matrix.mtx | cut -d " " -f 1-2,5 > {params.input_dir}/ma.txt

        #rm -r {params.input_dir}/spliced {params.input_dir}/unspliced {params.input_dir}/ambiguous
        mkdir -p {params.input_dir}/spliced {params.input_dir}/unspliced {params.input_dir}/ambiguous

        cat {params.input_dir}/mtx_header.txt {params.input_dir}/mtx_summary.txt {params.input_dir}/ms.txt > {params.input_dir}/spliced/matrix.mtx
        cat {params.input_dir}/mtx_header.txt {params.input_dir}/mtx_summary.txt {params.input_dir}/mu.txt > {params.input_dir}/unspliced/matrix.mtx
        cat {params.input_dir}/mtx_header.txt {params.input_dir}/mtx_summary.txt {params.input_dir}/ma.txt > {params.input_dir}/ambiguous/matrix.mtx

        cp {params.input_dir}/features.tsv {params.input_dir}/genes.tsv
        cp {params.input_dir}/genes.tsv {params.input_dir}/barcodes.tsv {params.input_dir}/spliced/
        cp {params.input_dir}/genes.tsv {params.input_dir}/barcodes.tsv {params.input_dir}/unspliced/
        cp {params.input_dir}/genes.tsv {params.input_dir}/barcodes.tsv {params.input_dir}/ambiguous/

        rm {params.input_dir}/*.txt
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
