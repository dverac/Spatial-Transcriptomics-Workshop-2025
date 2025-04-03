# GeoMx
### Author: Diana Vera Cruz

## Dataset

The dataset is the same used for the GeomxTools tutorial and is used 
throughout all the scripts for this workshop.
FFPE and FF tissue sections from 4 diabetic kidney disease (DKD) and 3
healthy kidney sections [Merritt et al.,
2020](https://pubmed.ncbi.nlm.nih.gov/32393914/), processed using the
GeoMx Digital Spatial Profiler (DSP) platform. The ROIs were selected to
focus on tubules or glomeruli regions.

- Glomeruli: Each ROI defined as glomeruli contains a single
  glomerulus.

- Tubular: Each ROI contains multiple tubules and are further
  classified into distal (PanCK+) or proximal (PanCK-) AOIs.

**For this workshop, we simulated 3 columns: The number of
nuclei and X and Y coordinates per ROI, since the initial data did not
have these fields.**

Box-folder with the raw dataset (DCCs, PKC and annotation): [Kidney dataset](https://uchicago.box.com/s/4lpumf7ekhw8711192frmzt9ndf0cofb)

Counts, features, and metadata tables post-QC are in the **results** folder.

NanoStringGeoMxSet and SpatialExperiment  objects are in the **env** folder.

## Steps

For each script, there is an md and Rmd file.

### 1. QC

[1_geomx_setup_qc](codes/1_geomx_setup_qc.md): Segment and Probe QC.

### 2. Normalization and DE analysis. 

[2_geomx_norm_DE](codes/2_geomx_norm_DE.md) DE analysis using GeomxTools, and LMM to adjust for batch. 

[2B_geomx_standR_norm_DE](codes/2B_geomx_standR_norm_DE.md): Extended version, various normalization methods, batch removal methods and DE analysis using Limma, adjusting for batch in model or using batch-corrected data.

[2C_geomx_limma_norm_DE](codes/2C_geomx_limma_norm_DE.md): DE analysis using Limma, adjusting for batch in model. 

### Spatial Deconvolution
[geomx_SpatialDecon](codes/geomx_SpatialDecon.md): Cell type abundances per segment inferred from a reference set. 

### Conversion to Seurat
[geomx_Seurat](codes/geomx_Seurat.md): Convert your data to Seurat.


## R packages

### GeoMx key packages
  * `GeomxTools`
  * `SpatialExperiment`
  * `limma`
  * `standR`
  * `SpatialDecon`
  * `SingleCellExperiment`
  * `Seurat`

### Other packages used

  * `tidyverse`
  * `readxl`
  * `ggalluvial`
  * `ComplexHeatmap`
  * `networkD3`
  * `PCAtools`
  * `umap`
  * `limma`
  * `edgeR`
  * `ggrepel`
  * `harmony`


## Relevant sites

  * [Nanostring Spatial Organ Atlas](https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/)
  Contains a list of the organs and tissues that have been studied. The datasets can be downloaded from this site.
  
  * [Nanostring Spatial tissue book](https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/)
    Contains the description of various studies, a study summary, and the purpose, including an overview of the slides, the number of genes, and the separation between segments. 


