# GlycoNorm

This repository contains the R code necessary to replicate the findings reported in the paper [Benedetti et al., "_Systematic Evaluation of Normalization Methods for Glycomics Data Based on Performance of Network Inference_", Metabolites, 10:7 (2020)](https://www.mdpi.com/2218-1989/10/7/271).

This code was created with R version 4.0.1 and Rstudio version 1.3.959.

Upon sourcing, the script Main.R will automatically download the [unprocessed glycomics data](http://dx.doi.org/10.6084/m9.figshare.12581735) used in the paper from the figshare repository and generate the following files into the working directory:

- **Barplot_Korcula2.pdf** -> corresponding to **Figure 3** in the paper 
- **Barplot_Korcula.pdf** -> corresponding to **Figure S5** in the paper
- **Barplot_Split.pdf** -> corresponding to **Figure S6** in the paper
- **Barplot_Vis.pdf** -> corresponding to **Figure S7** in the paper
- **Barplot_CRC.pdf** -> corresponding to **Figure 4** in the paper

- **AgeAssociation_Korcula2.pdf** / **AgeAssociation_Korcula.pdf** / **AgeAssociation_Split.pdf** / **AgeAssociation_Vis.pdf** / **AgeAssociation_CRC.pdf** -> corresponding to **Figure S8** in the paper

- **AgeAssociation_Results.xlsx** -> corresponding to **Table S1** in the paper

- **Validation_Summary.pdf** / **Validation_Summary.xlsx** -> corresponding to **Table 3** in the paper


Notes:

1. The current version of the code performes _nboot=100_ bootstrapping to compute the confidence intervals of the barplots, while in the paper we used _nboot=1000_. Increasing the number of bootstrapping will substantially increase the runtime, which for _nboot=100_ is roughly 9 minutes on a MacBook Pro (macOS version 10.15.1) with a 2.3 GHz Quad-Core Intel Core i5 processor and 16GB of RAM. Results obtained with the default _nboot_ value might slightly differ from the ones reported in the paper due to the difference in bootstrapping samples, but will be qualitatively equivalent.

2. One of the glycomics datasets used in the paper is not included in the figshare file (LLS cohort). In order to get access to it, users need to sign a Data Transfer Agreement (see paper for details).
