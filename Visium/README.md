# Visium
### Author: Jason Shapiro

## Datasets

### Example 1: Mouse Brain
This dataset uses the anterior brain sample also used in the [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette) for spatial data.
The data are available with SeuratData, from [10X](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0), or
from a [Box folder](https://uchicago.box.com/s/zlslf6av9flfnx8b9gjw976qw9cpttu1) provided for this workshop. We also use the same scRNA reference data
from [Tasic et al (2016)](https://www.nature.com/articles/nn.4216).

### Topics covered in Example 1
#### 1. Loading Data
#### 2. Visualizing Features
#### 3. Normalization
#### 4. Dimensionality Reduction
#### 5. Integration with scRNA reference data

### R Packages required for Example 1
  - `Seurat`
  - `dplyr`
  - `RColorBrewer`
  - `ggplot2`
  - `patchwork`


### Example 2: Tumor Microenvironment
This dataset includes two samples taken from [Batf3+ DCs and the 4-1BB/4-1BBL axis are required at the effector phase in the tumor microenvironment for PD-1/PD-L1 blockade efficacy](https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00469-8) by Ziblat et al (2024).
The original data are available from GEO ([GSE238145](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238145)), and the two samples used here are provided through an accompanying [Box folder](https://uchicago.box.com/s/pbsklrnue31wvq10qohlh7x9hev1nfaq).
As part of this analysis, we also created a scRNA reference dataset using data from [Zheng et al (2024)](https://pubmed.ncbi.nlm.nih.gov/38428409/). Our reference is generated from a downsample of a new integration and
annotation of the original data, and is also available in the same Box folder.

### Topics covered in Example 2
#### 1. Merging Data
#### 2. Distance Analysis with Giotto
#### 3. Spatial correlations between features
#### 4. Deconvolution
#### 5. Spatial correlations between cell types

### R Packages required for Example 1
  - `Seurat`
  - `Giotto`
  - `dplyr`
  - `RColorBrewer`
  - `ggplot2`
  - `patchwork`
  - `corrplot`	 
