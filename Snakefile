import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.10.0")

##### load config and sample sheets #####
#

def warning(msg):
    FAIL = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    print(f"\n{BOLD}{FAIL}{msg}{ENDC}\n",file=sys.stderr)

try:
    configfile: "config.yaml"
except WorkflowError:
    warning("ERROR: config.yaml does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

try:
    samples = pd.read_csv(config["samples"], sep="\t", comment="#").set_index("sample", drop=False)
except FileNotFoundError:
    warning(f"ERROR: the samples file ({config['samples']}) does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

try:
    units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(["sample", "unit"], drop=False)
except FileNotFoundError:
    warning(f"ERROR: the units file ({config['units']}) does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

def confirm():
    prompt = "Continue anyway? [y/N] "
    while 1:
        sys.stdout.write(prompt)
        choice = input().lower()
        if choice in ['yes','y']:
            return True
        elif choice in ['','no','n']:
            return False
        else:
            sys.stdout.write("Please respond with 'yes' or 'no': ") 

def get_resource(rule,resource):
    try:
        return config["rules"][rule]["res"][resource]
    except KeyError:
        return config["rules"]["default"]["res"][resource]

def get_input_degs(wc):
    if config["rules"]["seurat_degs"]["params"]["selected_res"]:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        if config["rules"]["seurat_merge"]["params"]["perform"] == True:
            samples = samples + ['merged']
        file = expand("{OUTDIR}/seurat/{sample}/4_degs/seurat_degs.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_gs(wc):
    if config["rules"]["seurat_gs"]["params"]["geneset_collection"]:
        samples = [u.sample for u in units.itertuples()]
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        if config["rules"]["seurat_merge"]["params"]["perform"] == True:
            samples = samples + ['merged']
        file = expand("{OUTDIR}/seurat/{sample}/5_gs/seurat_complete.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_ti(wc):
    if config["rules"]["slingshot"]["params"]["perform"] == True:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        if config["rules"]["seurat_merge"]["params"]["perform"] == True:
            samples = samples + ['merged']
        file = expand("{OUTDIR}/slingshot/{sample}/6_traj_in/slingshot_sce_objects.RData", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_input_fa(wc):
    if config["rules"]["vision"]["params"]["perform"] == True:
        samples = [u.sample for u in units.itertuples()] 
        if config["rules"]["seurat_integration"]["params"]["perform"] == True:
            samples = samples + ['integrated']
        if config["rules"]["seurat_merge"]["params"]["perform"] == True:
            samples = samples + ['merged']
        file = expand("{OUTDIR}/vision/{sample}/7_func_analysis/vision_object.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_integration(wc):
    file = []

    samples =list(set([u.sample for u in units.itertuples()]))

    if config["rules"]["seurat_integration"]["params"]["perform"] == True:
        if len(samples) == 1:
            print("WARNING: found only one sample. Deactivating Seurat integration.")
            print("You can remove this warning by disabling Seurat integration in config.yaml.")
        else:
            file = expand("{OUTDIR}/seurat/integrated/2_normalization/seurat_normalized-pcs.rds",OUTDIR=OUTDIR)

    return file

def get_merge(wc):
    file = []

    samples =list(set([u.sample for u in units.itertuples()]))

    if config["rules"]["seurat_merge"]["params"]["perform"] == True:
        if len(samples) == 1:
            print("\nWARNING: found only one sample. Deactivating Seurat merge.")
            print("You can remove this warning by disabling Seurat merge in config.yaml.")
            if not confirm():
                quit()
        else:
            file = expand("{OUTDIR}/seurat/merged/1_preprocessing/seurat_post-qc.rds",OUTDIR=OUTDIR)

    return file

def get_integration_input_sm(wc):
    file = expand("{OUTDIR}/seurat/{unit.sample}/1_preprocessing/seurat_post-qc.rds", unit=units.itertuples(),OUTDIR=OUTDIR)
    file = list(set(file))
    return file

def get_merge_input(wc):
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
            if config["rules"]["seurat_integration"]["params"]["perform"] == True:
                samples = samples + ['integrated']
            if config["rules"]["seurat_merge"]["params"]["perform"] == True:
                samples = samples + ['merged']
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
    if config["rules"]["seurat_merge"]["params"]["perform"] == True:
        samples = samples + ['merged']
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
        get_merge,
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
