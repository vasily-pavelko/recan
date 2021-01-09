---
title: 'Recan: R package for analysis of recombination events in viral genomes'
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
date: 4 January 2021
bibliography: references.bib
---


# Summary
`recan` is an R version of the Python library that allows recombination event analysis through construction and exploration of similarity plots based on genetic distances between nucleotide sequences [@Babin2020]. The R version of `recan` is based on `Bio3D` and `Plotly` libraries. Package included two functions: `seqSim` and `scanSeqSim`. Both functions require a sequence alignment in fasta format as an input. The output is a dataframe contained sequences similarities and illustrated by a single `plotly` plot or a list of plots for every pair of sequences with a highlighted region of putative recombination.

# Statement of need
Nowadays the pandemic of SARS-COV-2 virus attracts attention to viral diseases and the evolution of its agents. Detection of recombinations is more reproducible to perform in an automatic pipeline, rather than use GUI [@Etherington2005; @Lole1999]. As R language [@R2018] is more widely used by biologists I ported the `recan` package from Python to R.
Another package `seqcombo` provides the same analysis [@Yu2021]. In my package data with genetic distances are available for further investigations. Also using a `plotly` library for visualization data allows excluding redundant plot lines from the figure. A new function `scanSeqSim` for searching areas of recombination was added. In `scanSeqSim` you can adjust threshold value and recombination length the minimal number of points in ‘cross-over’ expressed in the number of the sliding windows for reducing the number of false-positively detected recombination areas.  

# Testing and verification
To test the package, we used the same four sets of viral genomes from the article with a description of the Python version of `recan` [@Babin2020].
The resulting `seqSim` method execution time with the window size = 400 and shift parameters = 200 showed in table 1. Both R packages are faster than the Python version.

Table 1. Benchmarking of python and R version of recan package

|           | number of sequences |   bp   |   recan python   |      recan R     |     seqcombo R     |
|:---------:|:-------------------:|:------:|:----------------:|:----------------:|:------------------:|
|    HIV    |          25         |  3135  | 437 ms ± 7.74 ms | 152 ms ± 13.9 ms |  116 ms ± 16.6 ms  |
|    HCV    |          23         |  9431  | 579 ms ± 58.7 ms | 437 ms ± 63.5 ms | 132.8 ms ± 16.7 ms |
| Norovirus |          19         |  3366  | 648 ms ± 44.2 ms |  98 ms ± 6,42 ms |  94.1 ms ± 19.7 ms |
|    LSDV   |          14         | 150511 |  3.55 s ± 239 ms |  2.63 s ± 570 ms | 391.4 ms ± 40.2 ms |

Time execution test was performed using a laptop with 2 CPU cores and 8 Gb RAM.
With the `scanSeqSim` function, `recan` can be used to identify and analyze recombination events in a large subset of sequences in separated plots. A set of SARS-CoV-2 full genome sequences was used to probe this function [@Paraskevis2020].
The distance plots with recombination events detected by `seqSim` and `scanSeqSim` are shown in Figures 1-2.
![](https://raw.githubusercontent.com/vasily-pavelko/recan/master/plots_paper/SARS_COV_2_all.jpg)
_Figure 1. SARS-CoV-2 sequences._

![](https://raw.githubusercontent.com/vasily-pavelko/recan/master/plots_paper/SARS_COV_2_%2315.jpg)
_Figure 2. Example of scanSecSim function for SARS-CoV-2 sequences._



# Availability and implementation
Recan is supported on Linux, macOS and Windows. The package can be installed by `devtools` R package using command
```
install.packages("devtools")   
library(devtools)   
devtools::install_github("vasily-pavelko/recan")
```
The source code, guide and datasets are available on the GitHub repository (https://github.com/vasily-pavelko/recan). 


# Acknowledgement 
The author thanks Nataly Pavelko, Yuriy Babin and Alexander Sprygin for editing the manuscript.



# References



