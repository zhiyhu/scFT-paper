# scFT-paper

## Citation

This repo is associated with the preprint:


Zhiyuan Hu, Mara Artibani, Abdulkhaliq Alsaadi, Nina Wietek, Matteo Morotti, Laura Santana Gonzalez, Salma El-Sahhar, Mohammad KaramiNejadRanjbar, Garry Mallett, Tingyan Shi, Kenta Masuda, Yiyan Zheng, Kay Chong, Stephen Damato, Sunanda Dhar, Riccardo Garruto Campanile, Hooman Soleymani majd, Vincenzo Cerundolo, Tatjana Sauka-Spengler, Christopher Yau \*, Ahmed A. Ahmed \*. **The repertoire of serous ovarian cancer non-genetic heterogeneity revealed by single-cell sequencing of normal fallopian tube epithelial cells**. *bioRxiv* 672626; doi: https://doi.org/10.1101/672626


## File structure


* Rmd is the code file
* html is the Rmarkdown report

```
clustering
├── Github_clusteing_all_data.Rmd # Clustering all FT cells
├── Github_clusteing_all_data.html
├── Github_manuscript_clustering.Rmd # Clustering secretory cells
├── Github_manuscript_clustering.html
├── clincluster
│   ├── clincluster_functions.R # Main functions
│   └── function_getMarkers.R
├── visualisation_secretory.Rmd # Visualisation for the manuscript
└── visualisation_secretory.html



cnvQC # filtering cells by CNV
├── HoneyBadger_fresh_secretory_exprs20180706.Rmd
├── HoneyBadger_fresh_secretory_exprs20180706.html
├── P11528_tumour_FTE_SNPsCNVs20180711.Rmd
└── P11528_tumour_FTE_SNPsCNVs20180711.html



culture_effect # studying the effect of culture
├── M_S4_0913_culturing.html
├── M_S4_1010_PhenoPath.html
├── PhenoPath_181010.Rmd
└── culturing_180918.Rmd



deconvolution # deconvolution
├── DEanalysis_EMThigh_TCGA.Rmd
├── DEanalysis_EMThigh_TCGA.html
├── deconvolution_analysis.Rmd
└── deconvolution_analysis.html


```

More details are on the way.


