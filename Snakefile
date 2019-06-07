import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)

units = pd.read_csv(config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

##### target rules #####

rule all:
    input:
        f"{OUTDIR}/qc/multiqc_report.html",
        expand("{OUTDIR}/seurat/{unit.sample}/seurat_final.rds", unit=units.itertuples(),OUTDIR=OUTDIR)


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/analyse.smk"
