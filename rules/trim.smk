def get_fastq(wildcards):
    return list(units.loc[(wildcards.sample), "fq" + wildcards.group])

rule merge:
    input:
        get_fastq
    output:
        "{}/merged/{{sample}}.r{{group}}.fastq.gz".format(OUTDIR)
    log:
        "{}/merge/{{sample}}.r{{group}}.log".format(LOGDIR)
    benchmark:
        "{}/merge/{{sample}}.r{{group}}.bmk".format(LOGDIR)
    threads: get_resource("merge","threads")
    resources:
        mem=get_resource("merge","mem")
    shell:"""
        cat {input} > {output} 2> {log}
    """

rule cutadapt:
    input:
        fastq1="{}/merged/{{sample}}.r1.fastq.gz".format(OUTDIR),
        fastq2="{}/merged/{{sample}}.r2.fastq.gz".format(OUTDIR)
    output:
        fastq1="{}/trimmed/{{sample}}.r1.fastq.gz".format(OUTDIR),
        fastq2="{}/trimmed/{{sample}}.r2.fastq.gz".format(OUTDIR),
        qc="{}/cutadapt/{{sample}}.out".format(LOGDIR)
    params:
        config["params"]["cutadapt"]
    log:
        "{}/cutadapt/{{sample}}.err".format(LOGDIR)
    benchmark:
        "{}/cutadapt/{{sample}}.bmk".format(LOGDIR)
    threads: get_resource("cutadapt","threads")
    resources:
        mem=get_resource("cutadapt","mem")
    wrapper:
        "0.32.0/bio/cutadapt/pe"
