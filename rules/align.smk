rule align:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1="{}/trimmed/{{sample}}.r2.fastq.gz".format(OUTDIR),
        fq2="{}/trimmed/{{sample}}.r1.fastq.gz".format(OUTDIR),
        idx="{}/index/genomeParameters.txt".format(OUTDIR)
    output:
        # see STAR manual for additional output files
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    log:
        "{}/star/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/star/{{sample}}.bmk".format(LOGDIR)
    params:
        # path to STAR reference genome index
        index="out/index",
        # optional parameters
        extra="--sjdbGTFfile {} --soloCBwhitelist {} {}".format(config["ref"]["annotation"], config["whitelist"], config["params"]["star"])
    threads: 24
    resources:
        mem=64000
    conda: "../envs/star.yaml"
    wrapper: 
        "file:wrappers/star/align"
