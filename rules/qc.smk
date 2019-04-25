# RSEQC
#
rule fastqc:
    input:
        lambda wc: units.loc[(wc.sample,wc.unit)]['fq' + wc.read]
    output:
        html="{}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.html".format(OUTDIR),
        zip="{}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.zip".format(OUTDIR)
    threads: 4
    params: "-t 4"
    log:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.bmk".format(LOGDIR)
    wrapper:
        "0.32.0/bio/fastqc"

rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR),
        db=temp("{}/qc/rseqc/annotation.db".format(OUTDIR))
    log:
        "{}/rseqc/rseqc_gtf2bed.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_gtf2bed.bmk".format(LOGDIR)
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"
       
rule rseqc_junction_annotation:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.junctionanno.junction.bed".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_junction_annotation/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_junction_annotation/{{sample}}.bmk".format(LOGDIR)
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="{}/qc/rseqc/{{sample}}.junctionanno".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"
        
rule rseqc_junction_saturation:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.junctionsat.junctionSaturation_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_junction_saturation/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_junction_saturation/{{sample}}.bmk".format(LOGDIR)
    params:
        extra=r"-q 255", 
        prefix="{}/qc/rseqc/{{sample}}.junctionsat".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem=8000
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"

rule rseqc_stat:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
    output:
        "{}/qc/rseqc/{{sample}}.stats.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_stat/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_stat/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.infer_experiment.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_infer/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_infer/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
        
rule rseqc_innerdis:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.inner_distance_freq.inner_distance.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_innerdis/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_innerdis/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.inner_distance_freq".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"

rule rseqc_readdis:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readdistribution.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readdis/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readdis/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_readdup:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readdup.DupRate_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readdup/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readdup/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.readdup".format(OUTDIR)
    resources:
        mem=32000
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"
        
rule rseqc_readgc:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readgc.GC_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readgc/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readgc/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.readgc".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"
        
rule multiqc:
    input:
        expand("{OUTDIR}/qc/fastqc/{unit.sample}.{unit.unit}.r{read}_fastqc.zip", unit=units.itertuples(), read=('1','2'), OUTDIR=OUTDIR),
        expand("{LOGDIR}/cutadapt/{unit.sample}.out", unit=units.itertuples(), LOGDIR=LOGDIR),
        expand("{OUTDIR}/star/{unit.sample}/Aligned.sortedByCoord.out.bam", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.junctionanno.junction.bed", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.infer_experiment.txt", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.stats.txt", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.inner_distance_freq.inner_distance.txt", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.readdistribution.txt", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.readdup.DupRate_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{OUTDIR}/qc/rseqc/{unit.sample}.readgc.GC_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR),
        expand("{LOGDIR}/rseqc/rseqc_junction_annotation/{unit.sample}.log", unit=units.itertuples(), LOGDIR=LOGDIR)
    output:
        "{}/qc/multiqc_report.html".format(OUTDIR)
    params:
        "--config res/config/multiqc_config.yaml"
    benchmark:
        "{}/multiqc.bmk".format(LOGDIR)
    log:
        "{}/multiqc.log".format(LOGDIR)
    wrapper:
        "0.32.0/bio/multiqc"
