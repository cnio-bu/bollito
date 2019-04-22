rule align:
    input:
        #star needs the barcoding read (read 1) to be in the second position
        fq1="out/trimmed/{sample}.r2.fastq.gz",
        fq2="out/trimmed/{sample}.r1.fastq.gz",
        idx="out/index/genomeParameters.txt"
    output:
        # see STAR manual for additional output files
        "out/star/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        "log/star/{sample}.log"
    benchmark:
        "log/star/{sample}.bmk"
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
