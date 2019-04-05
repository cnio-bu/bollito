rule star_idx:
    shadow:"shallow"
    input:
        config["ref"]["sequence"]
    output:
        f="out/star_idx/genomeParameters.txt",
        d=directory("out/star_idx")
    threads: 1
    conda: "envs/star.yaml"
    resources:
        mem = 64000
    log:
        stdout="log/star_idx/{genome}.out",
        stderr="log/star_idx/{genome}.err"
    shell:"""
        STAR --runMode genomeGenerate --genomeDir {out.d} --genomeFastaFiles {input} > {log.stdout} 2> {log.stderr}
    """
