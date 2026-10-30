# Tumor immune microenvironment and response to TIL therapy in melanoma

Single-cell RNA-seq re-analysis of CD45+ immune cells from metastatic melanoma tumors sampled before and after tumor-infiltrating lymphocyte (TIL) adoptive cell therapy, asking which immune cell states at baseline differ between responders and non-responders.

- **Data:** GEO [GSE221553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221553) (22 libraries, 13 patients; part of SuperSeries GSE222448)
- **Source study:** Barras et al., *Science Immunology* 2024, [doi:10.1126/sciimmunol.adg7995](https://doi.org/10.1126/sciimmunol.adg7995). Patient response calls come from that paper's Data file S1 (best overall response per patient).
- **Stack:** Python, Scanpy, scVI (scvi-tools), PyDESeq2

## Pipeline

| Notebook | Step |
|---|---|
| [01_qc](notebooks/01_qc.ipynb) | Per-sample QC, MAD-based outlier flagging, Scrublet doublets. 93,522 cells in, 87,486 kept |
| [02_integration](notebooks/02_integration.ipynb) | scVI batch integration (batch = library), UMAP, Leiden clustering (22 clusters) |
| [03_annotation](notebooks/03_annotation.ipynb) | Marker-based cluster annotation with stated confidence per cluster |
| [04_responder_comparison](notebooks/04_responder_comparison.ipynb) | Responder vs non-responder composition, signatures and pseudobulk DE |

Key figures are in [results/figures](results/figures) and tables in [results/tables](results/tables).

## Design decisions worth knowing

- **Patient is the unit of analysis.** All tests compare 6 responders against 7 non-responders at patient level; cells are never treated as independent replicates. Differential expression is pseudobulk (PyDESeq2).
- **Responder definition (this project's grouping):** best overall response of CR or PR (patients 2, 3, 7, 8, 9, 13) is called responder; SD or PD is called non-responder. The paper's own responder definition was not checked against this grouping. Patients 7, 8 and 9 had a partial response but progressed by month 6, so the grouping is not clear-cut.
- **Composition is computed on immune cells only.** Several post-therapy libraries contain many melanoma, fibroblast and endothelial cells (e.g. patient 5 post-ACT is 96% non-immune), which would otherwise distort proportions. Libraries with fewer than 200 immune cells are excluded (patient 5 post-ACT).
- **Pre-specified signatures.** The gene sets tested in section 3 of notebook 04 were written from the paper's stated findings and standard marker genes before any responder vs non-responder comparison was run, so those four tests are a replication check, not a search.

## Results

With 13 patients the study is underpowered, and **no comparison reaches significance after multiple-testing correction.**

- **Baseline composition:** no cell type differs (lowest FDR 0.73). The largest raw differences are cytotoxic/exhausted-like CD8 T cells (median 26% of immune cells in responders vs 6% in non-responders) and C1Q+ macrophages (2% vs 29%), both highly variable between patients.
- **Signatures reported in the paper (CD8 T cells; myeloid cells):** exhaustion, costimulation and type I interferon scores are higher in responders, matching the paper's direction (unadjusted p = 0.07, 0.10, 0.05; FDR 0.135 each). CD8 cytotoxicity does not differ (p = 0.53).
- **Pre vs post therapy:** only six patients have usable paired libraries and changes are heterogeneous, so the comparison is descriptive.
- **CD8 pseudobulk DE:** 17 genes at FDR < 0.1; treated as exploratory.
- **Myeloid pseudobulk DE is confounded and should not be interpreted as myeloid biology.** Genes called higher in responders include T/NK genes (CCL5, GZMH, PRF1, TRBV6-6). A T-cell gene score computed inside myeloid cells is about 3.4x higher in responders (median 0.21 vs 0.06, p = 0.004), consistent with T-cell RNA contamination of myeloid cells in T-cell-rich tumors.

## Limitations and next steps

- Ambient-RNA removal (SoupX, DecontX or CellBender) is needed before the myeloid analysis can be trusted; it may also affect the myeloid interferon signature.
- Cluster annotation is manual at cluster level. Three clusters are low-confidence (low-quality lymphocytes, n=1,112; a 51-cell CD8 cluster; an 11-cell myeloid cluster) and are excluded from the statistics.
- Post-therapy samples are few, and several contain many non-immune cells, suggesting they were not purely CD45-sorted (the paper's supplement notes that some 5' libraries used viable-cell rather than CD45+ sorting). This limits conclusions about treatment effects.
- The paper's TCR-seq data (GSE222750) could add clonotype-level analysis.

## Reproduce

Requires [uv](https://docs.astral.sh/uv/) and about 4 GB of free disk space (environment ~1.3 GB, data ~2 GB).

```bash
uv sync
uv run python scripts/download_data.py      # ~165 MB from NCBI GEO
uv run python scripts/build_anndata.py      # combines the 22 libraries
cd notebooks
uv run jupyter lab                          # run 01 -> 04 in order
```

Notebook 02 trains scVI on CPU and takes roughly 30-45 minutes on a laptop. Intermediate data objects are written to `data/processed/` and are not tracked by git.

Notebook 02 uses a fixed random seed, but scVI training on CPU is not guaranteed to be identical between runs or machines, so a rerun may renumber the Leiden clusters. Notebook 03 maps cluster numbers to cell-type labels, so after rerunning 02, re-check the marker evidence in `results/tables/cluster_markers.csv` and update the label table before trusting notebooks 03-04. The results shown in this repository come from a single complete run.
