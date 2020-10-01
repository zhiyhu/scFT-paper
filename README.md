# scFT-paper

## Citation

This repo is associated with the publication:

Zhiyuan Hu, Mara Artibani, Abdulkhaliq Alsaadi, ..., Tatjana Sauka-Spengler, Christopher Yau \*, Ahmed A. Ahmed \*. (2020). **The repertoire of serous ovarian cancer non-genetic heterogeneity revealed by single-cell sequencing of normal fallopian tube epithelial cells**. *Cancer Cell*. Volume 37, Issue 2: p226-242. doi: https://doi.org/10.1016/j.ccell.2020.01.003

## Interactive visualisation and cell annotations

Please visit our [**Cell Browser**](https://ovariancancercell.github.io)

**Notes:** The **cell annotation files** can also be downloaded from our Cell Browser, e.g. [here](https://ovariancancercell.github.io/sampleSecretoryCells/meta.tsv) for secretory cell annotations. 

## File description

* Rmd is the code file
* html is the Rmarkdown report

There are four file folders corresponding to different parts in the manuscript.

### 1. Culture effects

We used differential expression (DE) analysis and pseudotime analysis to compare the freshly dissociated cells and cultured cells. 

* In `culture_effect/culturing_180918.Rmd` and `culture_effect/M_S4_0913_culturing.html` we used **DE analysis** and **pathway analysis** to investigate the difference between cells from various sources.

* In `culture_effect/PhenoPath_181010.Rmd` and `culture_effect/M_S4_1010_PhenoPath.html` we used **psuedotime ananlysis** ([PhenoPath](https://github.com/kieranrcampbell/phenopath), Campbell and Yau) to dig deeper.


### 2. QC by CNVs

Before we entering the "true" analysis, we must do some QC steps to avoid the inclusion of tumour cells into our analysis. A key characteristics of HGSOC cells is the frequent copy number variants (CNVs), which is similar to the glioblastoma cells.

* In `cnvQC/HoneyBadger_fresh_secretory_exprs20180706.Rmd` and its report `cnvQC/HoneyBadger_fresh_secretory_exprs20180706.html` you can see how we use [HoneyBadger](https://jef.works/HoneyBADGER/) (Fan et al., 2018) to infer the CNV from cells dissociated from FT of cancer patients. 

* In `cnvQC/P11528_tumour_FTE_SNPsCNVs20180711.Rmd` and its report `cnvQC/P11528_tumour_FTE_SNPsCNVs20180711.html`, we revealed some results that were not included in the manuscript. By comparing the SNVs called from the scRNA-seq data and the ones called from WES data, we found that the cells from pt11528 carrying CNVs also harhoured the pathological p53 mutation, indicating that they are either early lesion or metastasis. 


### 3. Clustering

The part contains the some key coding for the manuscript. 

* In `clustering/Github_clusteing_all_data.Rmd` and its report `clustering/Github_clusteing_all_data.html`, we first clustered all the FT cells from cancer patients, identifying major FTE cell types.

* In `clustering/Github_manuscript_clustering.Rmd` and its report `clustering/Github_manuscript_clustering.html`, we further clustered the secretory cells into fine-grained subtypes.

* In `clustering/visualisation_secretory.Rmd` and its report `clustering/visualisation_secretory.html`, you will see how the plots were produced for the manuscript.

* In `clustering/data_integration.Rmd` and its report `clustering/data_integration.html`, we used [Seurat v3](https://satijalab.org/seurat/v3.1/integration.html) to integrate the secretory cells from cancer patients and from benign donors, in which the existence of the secretory subtypes was valdiated.


### 4. Deconvolution

In the last part, we used the information obtained by scRNA-seq to deconvolute TCGA, AOCS and other datasets from [CuratedOvarianData](http://bioconductor.org/packages/release/data/experiment/html/curatedOvarianData.html).

* In `deconvolution/deconvolution_analysis.Rmd` and its report `deconvolution/deconvolution_analysis.html`, we performed deconvolution and survival analysis. The deconvolution was conducted by using [Cibersort](https://cibersort.stanford.edu/) (Newman et al.). 

* In `deconvolution/DEanalysis_EMThigh_TCGA.Rmd` and its report `deconvolution/d DEanalysis_EMThigh_TCGA.html`, we studied the molecular characteristics of those EMT-high tumours that had worse prognosis.


## Video

A video explaining the bio finding of our work at https://youtu.be/AwKZVEtzjhs
