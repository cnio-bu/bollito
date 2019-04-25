rule index:
    shadow:"shallow"
    input:
        config["ref"]["sequence"]
    output:
        f="{}/index/genomeParameters.txt".format(OUTDIR)
    params:
        d="{}/index".format(OUTDIR)
    threads: 6
    conda: "../envs/star.yaml".format(OUTDIR)
    resources:
        mem = 64000
    benchmark:
        "{}/index/index.bmk".format(LOGDIR)
    log:
        stdout="{}/index/log.out".format(LOGDIR),
        stderr="{}/index/log.err".format(LOGDIR)
    shell:"""
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.d} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}
    """
