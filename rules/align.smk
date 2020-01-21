rule star_index:
    input:
        fasta = config["ref"]["fasta"]
    output:
        directory(config["ref"]["idx"])
    threads: get_resource("star_index","threads")
    params:
        extra = ""
    log:
        f"{LOGDIR}/star_index/index.log"
    wrapper:
        "0.49.0/bio/star/index"

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

rule star:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1=f"{OUTDIR}/merged/{{sample}}.r2.fastq.gz",
        fq2=f"{OUTDIR}/merged/{{sample}}.r1.fastq.gz",
        whitelist=config["whitelist"],
        index=config["ref"]["idx"]
    output:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    log:
        f"{LOGDIR}/star/{{sample}}.log"
    benchmark:
        f"{LOGDIR}/star/{{sample}}.bmk"
    params:
        index=config["ref"]["idx"],
        extra="--sjdbGTFfile {} --soloCBwhitelist {} {}".format(config["ref"]["annotation"], config["whitelist"], config["rules"]["star"]["params"])
    threads: get_resource("star","threads")
    resources:
        mem=get_resource("star","mem"),
        walltime=get_resource("star","walltime")
    wrapper: 
        "0.49.0/bio/star/align"
