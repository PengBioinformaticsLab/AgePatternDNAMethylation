# **Age-Related Patterns of DNA Methylation Changes**

*Click [here](https://www.biorxiv.org/content/10.1101/2024.12.10.627727v2) to view the latest version of the manuscript!*

## Abstract:

Epigenetic clocks have achieved significant success in aging research, but they often assume linear methylation changes with age and lack biological interpretability. Using data from 4,641 samples across 23 GEO datasets, we analyzed 1,557 CpGs from nine widely used clocks with minimal overlap, and identified consistent age-associated methylation patterns. We then identified 19,432 age-associated CpGs (aaCpGs) that were strongly correlated with age and showed high consistency between sexes, with faster methylation changes observed in males. Most aaCpGs were identified during early and late life stages, indicating accelerated epigenetic changes during development and aging. No specific genomic enrichment was observed. Clustering analysis revealed four distinct, non-linear age-related methylation trajectories. These findings underscore the complexity of epigenetic aging and suggest that current clocks may overlook important dynamic patterns, particularly after age 65. Incorporating these insights could improve the accuracy and biological relevance of future epigenetic clocks, especially for use across diverse age ranges and populations.

### Brief Overview of Pipeline:

1.  PreprocHelpers.Rmd -\> Preproc.Rmd *(Supp. Fig 1, Supp. Fig 2, Supp. Fig 5)*

2.  FilterHelpers.Rmd -\> Filtering.Rmd *(Supp. Fig 3, Supp. Fig 4, Table 1, Fig. 3, Fig. 4)*

3.  ClusteringHelpers.Rmd -\> Clustering.Rmd *(Fig. 2, Fig. 5)*

Each "Helper" file should be fun before running the main chunk to ensure dependencies are all loaded

Code is written in chunks such that intermediate files are saved to disk to allow for resuming from any point

## Contact:

[kchen24\@stanford.edu](mailto:kchen24@stanford.edu){.email}

Last updated: 6/16/2025
