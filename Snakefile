import pandas as pd
from snakemake.utils import min_version, validate
##### set minimum snakemake version #####
min_version("5.10.0")

##### load config and sample sheets #####
#
class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)

try:
    configfile: "config.yaml"
except WorkflowError:
    warning("ERROR: config.yaml does not exist or is incorrectly formatted. Please see the README file for details. Quitting now.")
    sys.exit(1)

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

try:
    samples = pd.read_csv(config["samples"], sep="\t", comment="#").set_index("sample", drop=False)
    validate(samples, schema="schemas/samples.schema.yaml")
except FileNotFoundError:
    warning(f"ERROR: the samples file ({config['samples']}) does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

try:
    if config["input_type"] == "fastq":
    	units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index(["sample", "unit"], drop=False)
    	validate(units, schema="schemas/units_fastq.schema.yaml")
        units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
    if config["input_type"] == "matrix":
        units = pd.read_csv(config["units"], dtype=str, sep="\t", comment="#").set_index("sample", drop=False)
        if config["technology"] == "10x": 
            validate(units, schema="schemas/units_matrix_10x.schema.yaml")
        elif config["technology"] == "standard":   
            validate(units, schema="schemas/units_matrix_standard.schema.yaml")
except FileNotFoundError:
    warning(f"ERROR: the units file ({config['units']}) does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

if config["single_samples"] == True and any([config["parameters"]["seurat_merge"]["enabled"], config["parameters"]["seurat_integration"]["enabled"]]) == False:
    warning("WARNING: 'single_samples' option is set as True while not asking for merge or integration.\nSetting 'single_samples' option to False.")


def confirm():
    prompt = f"{ansitxt.BOLD}Continue anyway? [y/N]{ansitxt.ENDC} "
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
        return config["resources"][rule][resource]
    except KeyError:
        return config["resources"]["default"][resource]

def process_samples(samples):
    integrated = []
    merged = []

    if config["single_samples"] == True:
        if any([config["parameters"]["seurat_merge"]["enabled"], config["parameters"]["seurat_integration"]["enabled"]]) == False:
            return samples
        else:
            if config["parameters"]["seurat_merge"]["enabled"] == True:
                merged = ["merged"]
            if config["parameters"]["seurat_integration"]["enabled"] == True:
                integrated = ["integrated"]
            return merged + integrated
    else:
        return samples


def get_output_degs(wc):
    if config["parameters"]["seurat_degs"]["enabled"] == True:
        samples = process_samples([u.sample for u in units.itertuples()])
        file = expand("{OUTDIR}/seurat/{sample}/4_degs/seurat_degs.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file


def get_output_gs(wc):
    if config["parameters"]["seurat_gs"]["enabled"] == True:
        samples = process_samples([u.sample for u in units.itertuples()])
        file = expand("{OUTDIR}/seurat/{sample}/5_gs/seurat_complete.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_output_ti(wc):
    if config["parameters"]["slingshot"]["enabled"] == True:
        samples = process_samples([u.sample for u in units.itertuples()])
        file = expand("{OUTDIR}/slingshot/{sample}/6_traj_in/slingshot_sce_objects.RData", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_output_fa(wc):
    if config["parameters"]["vision"]["enabled"] == True:
        samples = process_samples([u.sample for u in units.itertuples()])
        file = expand("{OUTDIR}/vision/{sample}/7_func_analysis/vision_object.rds", sample=samples,OUTDIR=OUTDIR)
    else:
        file = []
    return file

def get_integration(wc):
    file = []

    samples =list(set([u.sample for u in units.itertuples()]))

    if config["parameters"]["seurat_integration"]["enabled"] == True:
        if len(samples) == 1:
            print(f"\nWARNING: found only one sample. {ansitxt.RED}{ansitxt.BOLD}Deactivating Seurat integration.{ansitxt.ENDC}")
            print("You can remove this warning by disabling Seurat integration in config.yaml.")
            if not confirm():
                sys.exit(1)
            else:
                config["parameters"]["seurat_integration"]["enabled"] = False
        else:
            file = expand("{OUTDIR}/seurat/integrated/2_normalization/seurat_normalized-pcs.rds",OUTDIR=OUTDIR)

    return file

def get_merge(wc):
    file = []

    samples =list(set([u.sample for u in units.itertuples()]))

    if config["parameters"]["seurat_merge"]["enabled"] == True:
        if len(samples) == 1:
            print(f"\nWARNING: found only one sample. {ansitxt.RED}{ansitxt.BOLD}Deactivating Seurat merge.{ansitxt.ENDC}")
            print("You can remove this warning by disabling Seurat merge in config.yaml.")
            if not confirm():
                sys.exit(1)
            else:
                config["parameters"]["seurat_merge"]["enabled"] = False
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
        if config["parameters"]["velocyto"]["enabled"] == True:
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
        if config["parameters"]["velocyto"]["enabled"] == True:
            samples = process_samples([u.sample for u in units.itertuples()])
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

def get_output_normalization(wc):
    samples = process_samples([u.sample for u in units.itertuples()])
    file = expand("{OUTDIR}/seurat/{sample}/2_normalization/seurat_normalized-pcs.rds", sample=samples,OUTDIR=OUTDIR)
    return file

def get_output_find_clus(wc):
    samples = process_samples([u.sample for u in units.itertuples()])
    file = expand("{OUTDIR}/seurat/{sample}/3_clustering/seurat_find-clusters.rds", sample=samples,OUTDIR=OUTDIR)
    return file

def get_anndata_check(wc):
    samples = process_samples([u.sample for u in units.itertuples()])
    file = expand("{OUTDIR}/scanpy/{sample}/check.txt", sample=samples,OUTDIR=OUTDIR)
    return file

def get_output_qc(wc):
    samples = [u.sample for u in units.itertuples()] 
    if config["parameters"]["seurat_merge"]["enabled"] == True:
        samples = samples + ['merged']
    file = expand("{OUTDIR}/seurat/{sample}/1_preprocessing/seurat_post-qc.rds", sample=samples,OUTDIR=OUTDIR)
    return file



rule all:
    input:
        get_integration,
        get_merge,
        get_multiqc,
        get_output_find_clus,
	    get_output_degs, 
        get_output_gs,
        get_output_ti,
        get_output_fa,
        do_velocity,
        get_velocity_matrices,
        get_anndata_check


rule expression_matrix:
    input:
        get_velocity_matrices,
        get_multiqc

rule qc_expression_matrix:
    input:
        get_output_qc,
        get_merge,
        get_multiqc

rule normalized_expression_matrix:
    input:
        get_output_normalization,
        get_integration,
        get_merge,
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
