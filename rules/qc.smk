# RSEQC

rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed="out/qc/rseqc/annotation.bed",
        db=temp("out/qc/rseqc/annotation.db")
    log:
        "log/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"
       
rule rseqc_junction_annotation:
    input:
        bam="out/star/{sample}/Aligned.sortedByCoord.out.bam",
        bed="out/qc/rseqc/annotation.bed"
    output:
        "out/qc/rseqc/{sample}.junctionanno.junction.bed"
    priority: 1
    log:
        "log/rseqc/rseqc_junction_annotation/{sample}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="out/qc/rseqc/{sample}.junctionanno"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"
        
rule rseqc_junction_saturation:
    input:
        bam="out/star/{sample}/Aligned.sortedByCoord.out.bam",
        bed="out/qc/rseqc/annotation.bed"
    output:
        "out/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "log/rseqc/rseqc_junction_saturation/{sample}.log"
    params:
        extra=r"-q 255", 
        prefix="out/qc/rseqc/{sample}.junctionsat"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"

rule rseqc_stat:
    input:
        "out/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "out/qc/rseqc/{sample}.stats.txt"
    priority: 1
    log:
        "log/rseqc/rseqc_stat/{sample}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="out/star/{sample}/Aligned.sortedByCoord.out.bam",
        bed="out/qc/rseqc/annotation.bed"
    output:
        "out/qc/rseqc/{sample}.infer_experiment.txt"
    priority: 1
    log:
        "log/rseqc/rseqc_infer/{sample}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
        
rule rseqc_innerdis:
    input:
        bam="out/star/{sample}/Aligned.sortedByCoord.out.bam",
        bed="out/qc/rseqc/annotation.bed"
    output:
        "out/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "log/rseqc/rseqc_innerdis/{sample}.log"
    params:
        prefix="out/qc/rseqc/{sample}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"

rule rseqc_readdis:
    input:
        bam="out/star/{sample}/Aligned.sortedByCoord.out.bam",
        bed="out/qc/rseqc/annotation.bed"
    output:
        "out/qc/rseqc/{sample}.readdistribution.txt"
    priority: 1
    log:
        "log/rseqc/rseqc_readdis/{sample}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_readdup:
    input:
        "out/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "out/qc/rseqc/{sample}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "log/rseqc/rseqc_readdup/{sample}.log"
    params:
        prefix="out/qc/rseqc/{sample}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"
        
rule rseqc_readgc:
    input:
        "out/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "out/qc/rseqc/{sample}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "log/rseqc/rseqc_readgc/{sample}.log"
    params:
        prefix="out/qc/rseqc/{sample}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"
        
rule multiqc:
    input:
        expand("out/star/{unit.sample}/Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.junctionanno.junction.bed", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.infer_experiment.txt", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.stats.txt", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.readdistribution.txt", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("out/qc/rseqc/{unit.sample}.readgc.GC_plot.pdf", unit=units.itertuples()),
        expand("log/rseqc/rseqc_junction_annotation/{unit.sample}.log", unit=units.itertuples())
    output:
        "out/qc/multiqc_report.html"
    log:
        "log/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"
