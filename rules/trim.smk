def get_fastq(wildcards):  
    return list(units.loc[(wildcards.sample), "fq" + wildcards.group])

rule merge:
    input:
        get_fastq
    output:
        f"{OUTDIR}/merged/{{sample}}.r{{group}}.fastq.gz"
    log:
        f"{OUTDIR}/merge/{{sample}}.r{{group}}.log"
    benchmark:
        f"{OUTDIR}/merge/{{sample}}.r{{group}}.bmk"
    threads: get_resource("merge","threads")
    resources:
        mem=get_resource("merge","mem"),
        walltime=get_resource("merge","walltime")
    shell:"""
        cat {input} > {output} 2> {log}
    """

rule cutadapt:
    input:
        fastq1=f"{OUTDIR}/merged/{{sample}}.r1.fastq.gz",
        fastq2=f"{OUTDIR}/merged/{{sample}}.r2.fastq.gz"
    output:
        fastq1=f"{OUTDIR}/trimmed/{{sample}}.r1.fastq.gz",
        fastq2=f"{OUTDIR}/trimmed/{{sample}}.r2.fastq.gz",
        qc=f"{OUTDIR}/cutadapt/{{sample}}.out"
    params:
        config["rules"]["cutadapt"]["params"]
    log:
        f"{LOGDIR}/cutadapt/{{sample}}.err"
    benchmark:
        f"{LOGDIR}/cutadapt/{{sample}}.bmk"
    threads: get_resource("cutadapt","threads")
    resources:
        mem=get_resource("cutadapt","mem"),
        walltime=get_resource("cutadapt","walltime")
    wrapper:
        "0.35.0/bio/cutadapt/pe"
