---
title: 'Recan: R-based tool for detection of recombination in viral genomes'
tags:
  - R
  - virology
  - recombination
authors:
- name: Vasily Pavelko
  orcid: 0000-0001-8142-1621
  affiliation: "1"
affiliations: 
- name: International Biotechnology Center Generium, Vladimir region, Volginskiy, Russian Federation.
  index: 1
date: 16 January 2021
bibliography: references.bib
---


# Summary
`Recan` is an R version of the Python library that allows the identification of recombination events through construction and exploration of similarity plots based on genetic distances between nucleotide sequences [@Babin2020]. The R version of `recan` is based on `Bio3D` and `Plotly` libraries. The package includes two functions: `seqSim` and `scanSeqSim`. Both functions require a sequence alignment in fasta format as an input. The output object is a dataframe containing sequence similarities, which is presented as a single `plotly` plot or a list of plots for every pair of sequences with putative recombination sites.

# Statement of need
Genomic changes play an important role in virus evolution through recombination as a major driving force. The identification of a recombination event is more reproducible when performed in an automatic pipeline, as compared to GUI [@Etherington2005; @Lole1999]. R language [@R2018] is the most widely used programming language in bioinformatics. In this regard, we ported the `recan` package from Python to R to make the process easier.
Another R package `seqcombo` provides the same analysis [@Yu2021]. The main difference is that our package data with genetic distances can be further used in downstream applications. Moreover, using a `plotly` library for data visualization allows excluding redundant plot lines from the outoutfigure. A new function `scanSeqSim` was also added for searching areas of recombination . In `scanSeqSim` the user can manually adjust threshold values and recombination lengths of ‘cross-over’ region to reduce the number of false-positive recombination areas. 

# Testing and verification
To validate the package, we selected the same four sets of viral genomes used in the Python version of `recan` [@Babin2020].
The resulting `seqSim` method execution time with the window size of 400 and shift parameters of 200 are shown in Table 1. Both R packages have shorter execution times than the Python version.

Table 1. Comparison of the python and R versions of `recan` package with `seqcombo`.

|           | number of  sequences |   bp   | recan python, ms | recan R, ms | seqcombo R, ms |
|:---------:|:--------------------:|:------:|:----------------:|:-----------:|:--------------:|
|    HIV    |          25          |  3135  |    437 (7.74)    |  152 (13.9) |   116 (16.6)   |
|    HCV    |          23          |  9431  |    579 (58.7)    |  437 (63.5) |  132.8 (16.7)  |
| Norovirus |          19          |  3366  |    648 (44.2)    | 98.5 (6.42) |   94.1 (19.7)  |
|    LSDV   |          14          | 150511 |    3558 (239)    |  2636 (570) |   391 (40.2)   |

Time execution test was performed using a laptop with 2 CPU cores and 8 Gb RAM.
With the `scanSeqSim` function, `recan` is capable of identifying and stripping evidence of recombination from sequence alignments in a large subset of sequences in separated plots. A set of SARS-CoV-2 full genome sequences was used to validate this function [@Paraskevis2020].
The distance plots with recombination events detected by `seqSim` and `scanSeqSim` are shown in Figures 1-2.

![](https://raw.githubusercontent.com/vasily-pavelko/recan/master/plots_paper/SARS_COV_2_all.jpg)
_Figure 1. SARS-CoV-2 sequences._

![](https://raw.githubusercontent.com/vasily-pavelko/recan/master/plots_paper/SARS_COV_2_%2315.jpg)
_Figure 2. Example of scanSecSim function for SARS-CoV-2 sequences._



# Availability and implementation
Recan is run on Linux, macOS and Windows. The package can be installed by `devtools` R package using command line
```
install.packages("devtools")  
library(devtools)  
devtools::install_github("vasily-pavelko/recan")
```
The source code, guide, and datasets are available on the GitHub repository (https://github.com/vasily-pavelko/recan). 


# Acknowledgement 
The author thanks Nataly Pavelko, Yuriy Babin, Vladimir Bobkov and Alexander Sprygin for editing the manuscript.



# References




