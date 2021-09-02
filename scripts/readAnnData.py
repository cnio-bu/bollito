import sys
sys.stderr=open(snakemake.log[0], "a+")
sys.stdout=open(snakemake.log[0], "a+")

import scanpy
import warnings

anndata_path = snakemake.input["anndata_obj"]
outdir = snakemake.params["output_dir"]

a = scanpy.read_h5ad(filename = anndata_path)
check_path = outdir + "/check.txt"
print(check_path)
check = open(check_path, 'w')

'''
try: 

except:
    warnings.warn("AnnData can not be read.")
    sys.exit(1)

'''
sys.stderr.close()
sys.stdout.close()

