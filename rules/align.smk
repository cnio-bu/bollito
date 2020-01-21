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

rule star:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1="{}/{}/{{sample}}.r2.fastq.gz".format(OUTDIR,"merged" if config["rules"]["cutadapt"]["disabled"] else "trimmed"),
        fq2="{}/{}/{{sample}}.r1.fastq.gz".format(OUTDIR,"merged" if config["rules"]["cutadapt"]["disabled"] else "trimmed"),
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
