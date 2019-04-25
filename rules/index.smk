rule star_index:
    shadow:"shallow"
    input:
        config["ref"]["sequence"]
    output:
        f="{}/star_index/genomeParameters.txt".format(OUTDIR)
    params:
        d="{}/star_index".format(OUTDIR)
    conda: "../envs/star.yaml".format(OUTDIR)
    resources:
        mem=get_resource("star_index","mem")
    threads: get_resource("star_index","threads")
    benchmark:
        "{}/star_index/star_index.bmk".format(LOGDIR)
    log:
        stdout="{}/star_index/log.out".format(LOGDIR),
        stderr="{}/star_index/log.err".format(LOGDIR)
    shell:"""
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.d} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}
    """
