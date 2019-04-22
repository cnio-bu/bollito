def get_fastq(wildcards):
    return list(units.loc[(wildcards.sample), "fq" + wildcards.group])

rule merge:
    input:
        get_fastq
    output:
        "out/merged/{sample}.r{group}.fastq.gz"
    log:
        "log/merge/{sample}.r{group}.log"
    benchmark:
        "log/merge/{sample}.r{group}.bmk"
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
        qc="log/cutadapt/{sample}.out"
    params:
        config["params"]["cutadapt"]
    log:
        "log/cutadapt/{sample}.err"
    benchmark:
        "log/cutadapt/{sample}.bmk"
    wrapper:
        "0.32.0/bio/cutadapt/pe"
