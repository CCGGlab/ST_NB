Computational analyses underlying the results reported in *Siaw et al., Spatial transcriptomics exploration of the primary neuroblastoma microenvironment unveils novel paracrine interactions, 2024.* (available in [preprint]())

# Environment

Analysis was performed in a Conda environment. See *ST_NB.yml* for details.

# Experimental setup

Visium spatial transcriptomics (ST) was performed on six 5 µM tumor sections obtained from 3 FFPE embedded tumors from 2 high-risk NB patients:

-   **NB1**: 2-year-old boy, no MYCN-amplification, no ALK mutations
    -   Treatment-naïve tumors (*NB1Pre*): 2 sections
    -   Chemotherapy-treated tumors (*NB1Post*): 2 sections
    -   Analysed using Visium v1 workflow - NB2:
-   **NB2**: 2-year-old girl, MYCN amplfied, ALK-R1278Q mutation that disappeared after chemotherapy
    -   Chemotherapy-treated tumors (*NB2Post*): 2 sections
    -   Analysed using Visium v2-CytAssist workflow

# Data

## Data availability

Space Ranger output files (raw expression matrices and spatial metadata) and Seurat objects are available at [Zotero]()

## ST data processing

Quantification and processing of Visium fastqc files using Space Ranger

``` shell
scripts/other/process_spaceranger_count.sh
```

Adjusting for spot swapping in ST data using spotclean:

``` r
scripts/process_spotclean.R
```

Create seurat object

``` r
scripts/process_seurat.R
```

## Gene sets

[Reactome gene sets](http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/msigdb_v2023.1.Hs_files_to_download_locally.zip) were downloaded from MSigDB (v2023.1)

[PanglaoDB gene sets](https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz) were downloaded from PanglaoDB

Data available in **data/ST_NB.RData**:

-   *PDB_Hs_ls*: Panglao DB gene sets
-   *Rea_ls*: Reactome gene sets

## Experimental proliferation and migration data

Different NB cell lines were treated with CCL18 and NRTN and migration and proliferation was measured as described in the manuscript.

Data available in object *exp_val_ls* from **data/ST_NB.RData**

# Manuscript

To run the code below, first download the seurat R object *ST_NB_seurat.rds* from Zenodo in the data/ folder

## Sample overview and clustering

Seurat clusters (Fig. S2A-B)

``` r
source("scripts/manuscript_basic_clustering.R")
```

Annotated clusters

-   Plot H&E histology (Fig. 1B)
-   Spatial plot of low resolution clusters (Fig. 1C)
-   Plot UMAP: low resolution clusters, copy number score & adrenergic/mesenchymal cell state (Fig. 1E)
-   Spatial plot & UMAP plots of high resolution clusters (Fig. S4A-B)
-   Adrenergic/mesenchymal cell states: spatial plots & correlation plots (Fig. S5)

``` r
source("scripts/manuscript_clustering.R")
```

## DGE analysis

Dotplots (Fig. 1D, Fig. S2C, Fig. S4C) & DGE table (table S1)

``` r
# Perform DGE
source("scripts/manuscript_do_DGE.R")
# Create dotplots
source("scripts/manuscript_DGE_dotplots.R")
```

fGSEA analysis - Running score GSEA plots of Reactome steroid pathways (Fig. 2A) - fGSEA plots (Fig. S2D) & table (table S2)

``` r
# fGSEA analysis
source("scripts/manuscript_DGE_fGSEA.R")
# Plot
source("scripts/manuscript_DGE_fGSEA_visualize.R")
```

## NB2Post cluster analysis

-   Spatial UCell Reactome steroid pathway scores (Fig. 2B)
-   Spatial plots of steroid gene expression (Fig. 2C)

``` r
source("scripts/manuscript_AC.R")
```

## NB2Post ALKAL2-ALK & cell-cell communication analysis

-   Spatial plot of AC-like, NE1 and NE2-CA clusters in NB2Post (Fig. 3A)
-   Spatial plots of ALKAL2/ALK (Fig. 3B), NRTN/RET/GFRA2 (Fig. 4A) & PNMT/DLK1 expression (Fig. S4D)
-   Circus plot of predicted ligand-receptor interactions (Fig. 4A)
-   CellChat results (Table S3)

``` r
# CellChat analysis
source("scripts/manuscript_run_cellchat.R")
# Create plots
source("scripts/manuscript_NB2_interaction.R")
```

## Adrenal development analysis

-   Analysis of ALK & ALKAL2 expression in the single cell human adrenal development data of [Jansky et al](https://www.nature.com/articles/s41588-021-00806-1)

``` r
# Process Jansky data
source("scripts/manuscript_process_jansky.R")
# Analyse Jansky data
source("scripts/manuscript_jansky_ALK_ALKAL.R")
```

-   Analysis of ALK & ALKAL2 expression evolution during human adrenal development in in the adrenal atlas of [Del Valle et al](https://doi.org/10.12688/wellcomeopenres.11253.2)

``` r
# Process del Valle data
source("scripts/manuscript_process_DelValle.R")
# Analysis
source("scripts/manuscript_adrenal_development.R")
```

## NB1Post analysis of cluster NE4 & surrounding macrophages

-   Spatial plot of macrophages & NE4 cells in NB1Post section 2 (Fig. 5A)
-   Top 10 DE genes in NE4 (Fig. 5B) & top 10 DE ligand encoding genes in macrophage cluster (Fig. 5E)
-   Spatial plot of 11p gains
-   Spatial plots of WEE1/NKX6-1/NTRK1/CCL18/PITPNM/FGF1/FGFR1 (Fig. 5C,G),
-   Circus plot of predicted ligand-receptor interactions (Fig. 5F)

``` r
source("scripts/manuscript_NB1_interaction.R")
```

-   ST-derived copy number profiles for chr. 11 in NB1Post section 2 (Fig. 5D) & the complete genome in all samples (Fig. S6)

``` r
source("scripts/manuscript_CNV.R")
```

## Clinical validation of AC-like associations & survival analysis

The following analysis were performed on the RNA-Seq & associated clinical data of [Cangelosi et al., 2020](https://www.mdpi.com/2072-6694/12/9/2343):

-   Association of AC-like signature with
    -   ALKAL2 expression (Fig. 3F)
    -   NRTN expression (Fig. 4C)
-   Survival analysis of CCL18-PITPNM3 (Fig. 6D-F)

``` r
# Calculate AC-like score on cangelosi data
source("scripts/manuscript_process_cangelosi.R")
# Analysis
source("scripts/manuscript_clin_val.R")
```

## Experimental validation of NRTN & CCL18 effects on NB cells

-   NRTN effects on NB cell migration & proliferation (Fig. 4E-F)
-   CCL18 effects on NB cell migration & proliferation (Fig. 6B-C)

``` r
source("scripts/manuscript_exp_val.R")
```
