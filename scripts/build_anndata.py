"""Combine the 22 GSE221553 count matrices into one AnnData with sample/patient/response annotations."""
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "GSE221553"
OUT = ROOT / "data" / "processed" / "cd45_raw.h5ad"

patients = pd.read_csv(ROOT / "data" / "metadata" / "patients.csv").set_index("patient")


def read_counts(path):
    """Dense genes x cells TSV -> sparse cells x genes, read in row chunks to bound memory."""
    blocks, genes = [], []
    for chunk in pd.read_csv(path, sep="\t", index_col=0, chunksize=2000):
        genes.extend(chunk.index)
        blocks.append(sp.csr_matrix(chunk.values.astype(np.int32)))
        cells = chunk.columns
    return sp.vstack(blocks).T.tocsr(), pd.Index(genes), cells


adatas = []
for f in sorted(RAW.glob("*.tsv.gz")):
    gsm, patient, tp = re.match(r"(GSM\d+)_patient(\d+)_(T\d+)_CD45", f.name).groups()
    X, genes, cells = read_counts(f)
    a = ad.AnnData(X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    a.obs["gsm"] = gsm
    a.obs["patient"] = int(patient)
    a.obs["sample"] = f"patient{patient}_{tp}"
    a.obs["timepoint"] = "post" if tp == "T30" else "pre"
    a.obs["biopsy_site"] = 2 if tp == "T02" else 1
    adatas.append(a)
    print(f"{a.obs['sample'].iloc[0]:<18} {a.n_obs:>6} cells x {a.n_vars} genes")

same_genes = all(a.var_names.equals(adatas[0].var_names) for a in adatas)
print("identical gene sets across samples:", same_genes)
adata = ad.concat(adatas, join="inner" if not same_genes else "outer", merge="same")
adata.obs = adata.obs.join(patients[["best_overall_response", "response", "tls"]], on="patient")
adata.obs["patient"] = adata.obs["patient"].astype(str).astype("category")
for c in ("sample", "timepoint", "response", "best_overall_response", "gsm", "tls"):
    adata.obs[c] = adata.obs[c].astype("category")
adata.obs_names_make_unique()
adata.write_h5ad(OUT, compression="gzip")
print(adata)
print(adata.obs.groupby(["timepoint", "response"], observed=True).size())
