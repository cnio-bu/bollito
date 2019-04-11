def get_fastq(wildcards):
    return list(units.loc[(wildcards.sample), "fq" + wildcards.group])

rule merge:
    input:
        get_fastq
    output:
        "out/merged/{sample}.r{group}.fastq.gz"
    log:
        "log/merge/{sample}.r{group}.log"
    shell:"""
        cat {input} > {output} 2> {log}
    """

rule cutadapt:
    input:
        fastq1="out/merged/{sample}.r1.fastq.gz",
        fastq2="out/merged/{sample}.r2.fastq.gz"
    output:
        fastq1="out/trimmed/{sample}.r1.fastq.gz",
        fastq2="out/trimmed/{sample}.r2.fastq.gz",
        qc="out/trimmed/{sample}.qc.txt"
    params:
        config["params"]["cutadapt"]
    log:
        "log/cutadapt/{sample}.log"
    wrapper:
        "0.31.1/bio/cutadapt/pe"
