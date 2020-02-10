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

def input_merge(wildcards):  
    return list(units.loc[(wildcards.sample), "fq" + wildcards.group])

rule merge:
    input:
        input_merge
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

def input_star(sample,group):
    '''call the merge rule if there's more than one fastq, otherwise return the fastq directly'''
    fastqs = list(units.loc[(sample), f"fq{group}"])
    if len(fastqs) == 1:
        return fastqs[0]
    else:
        return f"{OUTDIR}/merged/{sample}.r{group}.fastq.gz"


def extra_params(wc):
    if config["technology"] == "10X":
        extra="--sjdbGTFfile {} --soloCBwhitelist {} {}".format(config["ref"]["annotation"], config["whitelist"], config["rules"]["star"]["params"])
        return extra
    elif config["technology"] == "Drop-seq":
        extra = "--sjdbGTFfile {} --soloCBwhitelist None {}".format(config["ref"]["annotation"], config["rules"]["star"]["params"])
        return extra
    else:
        raise ValueError('Specified technology is not valid.')

rule star:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1=lambda wc: input_star(wc.sample,2),
        fq2=lambda wc: input_star(wc.sample,1),
        index=config["ref"]["idx"]
    output:
        f"{OUTDIR}/star/{{sample}}/Aligned.sortedByCoord.out.bam"
    log:
        f"{LOGDIR}/star/{{sample}}.log"
    benchmark:
        f"{LOGDIR}/star/{{sample}}.bmk"
    params:
        index=config["ref"]["idx"],
        extra=extra_params
    threads: get_resource("star","threads")
    resources:
        mem=get_resource("star","mem"),
        walltime=get_resource("star","walltime")
    wrapper: 
        "0.49.0/bio/star/align"
