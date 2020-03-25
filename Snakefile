import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.10.0")

##### load config and sample sheets #####

try:
    configfile: "config.yaml"
except WorkflowError:
    quit(f"ERROR: the config file config.yaml does not exist. Please see the README file for details.")

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

try:
    samples = pd.read_csv(config["samples"], sep="\t", comment="#").set_index("sample", drop=False)
except FileNotFoundError:
    quit(f"ERROR: the samples file ({config['samples']}) does not exist. Please see the README file for details.")

try:
    units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(["sample", "unit"], drop=False)
except FileNotFoundError:
    quit(f"ERROR: the units file ({config['units']}) does not exist. Please see the README file for details.")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index



def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

## get rule all inputs depending on the enabled fields ## 
def get_input_degs(wc):
    if config["rules"]["seurat_degs"]["params"]["selected_res"]:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        file = expand("{OUTDIR}/seurat/{sample}/4_degs/seurat_degs.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_gs(wc):
    if config["rules"]["seurat_gs"]["params"]["geneset_collection"]:
        samples = [u.sample for u in units.itertuples()]
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        file = expand("{OUTDIR}/seurat/{sample}/5_gs/seurat_complete.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_ti(wc):
    if config["rules"]["slingshot"]["params"]["perform"] == True:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        file = expand("{OUTDIR}/slingshot/{sample}/6_traj_in/slingshot_sce_objects.RData", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_fa(wc):
    if config["rules"]["vision"]["params"]["perform"] == True:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        file = expand("{OUTDIR}/vision/{sample}/7_func_analysis/vision_object.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_integration(wc):
    samples =list(set([u.sample for u in units.itertuples()]))
    if len(samples) == 1 and config["rules"]["seurat_integration"]["params"]["perform"] == True:
        raise ValueError("\nATTENTION: Seurat integration must be deactivated when there is only one sample.\nChange integration perform value to 'False' in config.yaml file.")
    if config["rules"]["seurat_integration"]["params"]["perform"] == True:
        file = expand("{OUTDIR}/seurat/integrated/2_normalization/seurat_normalized-pcs.rds",OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_integration_input_sm(wc):
    file = expand("{OUTDIR}/seurat/{unit.sample}/1_preprocessing/seurat_post-qc.rds", unit=units.itertuples(),OUTDIR=OUTDIR)
    file = list(set(file))
    return file

def get_velocity_matrices(wc):
    if config["input_type"] == "matrix":
        file = []
    elif config["input_type"] == "fastq":
        if config["rules"]["velocyto"]["params"]["perform"] == True:
            file = expand("{OUTDIR}/star/{unit.sample}/Solo.out/Velocyto/raw/spliced/matrix.mtx", unit=units.itertuples(),OUTDIR=OUTDIR)
        else:
            file = []
    else: 
        file = []
    return file 

def do_velocity(wc):
    if config["input_type"] == "matrix":
        file = []
    elif config["input_type"] == "fastq":
        if config["rules"]["velocyto"]["params"]["perform"] == True:
            samples = [u.sample for u in units.itertuples()]
            file = expand("{OUTDIR}/velocyto/{sample}/8_RNA_velocity/seurat_velocity.rds", sample=samples,OUTDIR=OUTDIR)
        else:
            file = []
    else:
        file = []
    return file 

def get_multiqc(wc):
    if config["input_type"] == "fastq": 
        file = f"{OUTDIR}/qc/multiqc_report.html"
    else: 
        file = []
    return file
    
def seurat_input(wc): 
    if config["input_type"] == "matrix": 
        file = []
    elif config["input_type"] == "fastq":
        file = f"{OUTDIR}/star/{{sample}}/Solo.out/Gene/Summary.csv"
    return file

def get_input_find_clus(wc):
    samples = [u.sample for u in units.itertuples()] 
    if config["rules"]["seurat_integration"]["params"]["perform"] == True:
        samples = samples + ['integrated']
    file = expand("{OUTDIR}/seurat/{sample}/2_normalization/seurat_normalized-pcs.rds", sample=samples,OUTDIR=OUTDIR)
    return file

rule all:
    input:
        get_multiqc,
        get_input_degs, 
        get_input_gs,
        get_input_ti,
        get_input_fa,
        get_integration,
        do_velocity,
        get_velocity_matrices

rule expression_matrix:
    input:
        get_velocity_matrices,
        get_multiqc

rule normalized_expression_matrix:
    input:
        get_input_find_clus,
        get_integration,
        get_multiqc

rule differential_expression:
    input:
        get_input_degs,
        get_multiqc

rule functional_analysis:
    input:
        get_input_fa,
        get_multiqc

rule trajectory_inference:
    input:
        get_input_ti,
        get_multiqc

rule geneset_analysis:
    input:
        get_input_gs,
        get_multiqc

rule RNA_velocity:
    input:
        do_velocity,
        get_velocity_matrices,
        get_multiqc


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/conf/report.rst"

##### load rules #####

include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/analyse.smk"
include: "rules/velocyto.smk"
