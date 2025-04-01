# GeoMx
### Author: Diana Vera Cruz

## Dataset

The dataset is the same used for the GeomxTools tutorial, and is used 
trhoughout all the scripts for this workshop.
FFPE and FF tissue sections from 4 diabetic kidney disease (DKD) and 3
healthy kidney sections [Merritt et al.,
2020](https://pubmed.ncbi.nlm.nih.gov/32393914/), processed using the
GeoMx Digital Spatial Profiler (DSP) platform. The ROI were profiled to
focus on tubules or glomeruli regions.

- Glomeruli: Each ROIs defined as glomeruli contains a single
  glomerulus.

- Tubular: Each ROI contains multiple tubules, and are further
  classified into distal (PanCK+) or proximal (PanCK-) AOIs.

**For the purposes of the tutorial, we simulated 3 columns: Number of
nuclei, and X and Y coordinate per ROI, since the initial data did not
have this fields.**

Box folder with dataset: [Kidney dataset](https://uchicago.box.com/s/4lpumf7ekhw8711192frmzt9ndf0cofb)

## Steps

### 1. QC

### 2. Normalization and DE analysis

### Spatial Deconvolution

### Coersion to Seurat


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
  Contains the description of various studies, study summary and purpose, including an overview of the slides, the number of genes and overview of separation between segments. 


