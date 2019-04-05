def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule align:
    input:
        sample=get_trimmed,
        idx="out/index/genomeParameters.txt"
    output:
        # see STAR manual for additional output files
        "out/star/{sample}-{unit}/Aligned.out.bam"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index="out/index",
        # optional parameters
        extra="--sjdbGTFfile {} --soloCBwhitelist {} {}".format(config["ref"]["annotation"], config["whitelist"], config["params"]["star"])
    threads: 24
    wrapper:
        "0.31.1/bio/star/align"
