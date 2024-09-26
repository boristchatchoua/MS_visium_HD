import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt

import scanpy as sc 


sample = "SC54_probe2_long"
data_path = f"/home/t/tandrew6/btchatch/.$SCRATCH/spatial/Kerfoot_HD_Visium/{sample}/outs/binned_outputs"

out_dir = f"/home/t/tandrew6/btchatch/.$SCRATCH/spatial/Kerfoot_HD_Visium/pp/{sample}"


if no os.path.exist(out_dir):
    os.makedirs(out_dir)

# Create 'raw_data' subdirectory inside the output directory
raw_data_dir = os.path.join(out_dir, "raw_data")
if not os.path.exists(raw_data_dir):
    os.makedirs(raw_data_dir)    
    
    
    
# convert parquete files to csv 
def convert_parquet_to_csv(input_dir):
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".parquet"):
                df = pd.read_parquet(os.path.join(root, file))
                df.to_csv(os.path.join(root, file.replace(".parquet", "_list.csv")), index=False)
                print(f"Converted {os.path.join(root, file)} to csv")
                
convert_parquet_to_csv(data_path)

    
# list of different resolutions 
sdata = {}

# Load different resolutions and read data
for dir_name in os.listdir(data_path):
    # Check if the directory name has the expected format
    if "_" in dir_name:
        # Extract the resolution part
        res = dir_name.split("_")[1][1:3]
        # Construct the full path to the directory
        dir_path = os.path.join(data_path, dir_name)
        # Read the data from the directory
        adata = sc.read_visium(dir_path)
        # Store the data in the dictionary with the resolution as the key
        sdata[res] = adata



adata_02 = sdata["02"]
adata_08 = sdata["08"]
adata_16 = sdata["16"]

# save raw data
if adata_02 is not None:
    adata_02.write(os.path.join(raw_data_dir, f"raw_{sample}_02.h5ad"))
if adata_08 is not None:
    adata_08.write(os.path.join(raw_data_dir, f"raw_{sample}_08.h5ad"))
if adata_16 is not None:
    adata_16.write(os.path.join(raw_data_dir, f"raw_{sample}_16.h5ad"))

# preprocessing pipeline 
# pp pipline 

def pp(adata):
    """
    Preprocess adata for clustering.
    Perform normalization, log1p, highly variable genes selection, PCA, UMAP, and Leiden clustering.

    Parameters
    ----------
    adata : AnnData object 
        Annotated data matrix.
    id : str
    """
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.tl.pca(adata, use_highly_variable=True)
    print("PCA done!")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    res = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
#     res = [ 0.6, 0.7]
    for r in res:
        sc.tl.leiden(adata, resolution=r, key_added= f"GE_leiden_{r}")
    print("Leiden done!")
    return adata

print("starting pp for 02")
adata_02 = pp(adata_02)

# save pp data 
adata_02.write(os.path.join(out_dir, f"pp_{sample}_02.h5ad"))

print("starting pp for 08")
adata_08 = pp(adata_08)
adata_08.write(os.path.join(out_dir, f"pp_{sample}_08.h5ad"))


print("starting pp for 16")
adata_16 = pp(adata_16)
adata_16.write(os.path.join(out_dir, f"pp_{sample}_16.h5ad"))

print("done!")