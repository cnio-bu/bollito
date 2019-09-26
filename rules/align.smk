rule star:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1="{}/{}/{{sample}}.r2.fastq.gz".format(OUTDIR,"merged" if config["rules"]["cutadapt"]["disabled"] else "trimmed"),
        fq2="{}/{}/{{sample}}.r1.fastq.gz".format(OUTDIR,"merged" if config["rules"]["cutadapt"]["disabled"] else "trimmed"),
        whitelist=config["whitelist"]
    output:
        # see STAR manual for additional output files
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    log:
        "{}/star/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/star/{{sample}}.bmk".format(LOGDIR)
    params:
        # path to STAR reference genome index
        index=config["ref"]["idx"],
        # optional parameters
        extra="--sjdbGTFfile {} --soloCBwhitelist {} {}".format(config["ref"]["annotation"], config["whitelist"], config["rules"]["star"]["params"])
    threads: get_resource("star","threads")
    resources:
        mem=get_resource("star","mem"),
        walltime=get_resource("star","walltime")
    wrapper: 
        "68ce581/bio/star/align"
