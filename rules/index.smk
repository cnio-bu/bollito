rule index:
    shadow:"shallow"
    input:
        config["ref"]["sequence"]
    output:
        f="out/index/genomeParameters.txt"
    params:
        d="out/index"
    threads: 6
    conda: "../envs/star.yaml"
    resources:
        mem = 64000
    log:
        stdout="log/index/log.out",
        stderr="log/index/log.err"
    shell:"""
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.d} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}
    """
