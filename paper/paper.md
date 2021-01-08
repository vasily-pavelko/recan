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
`recan` is an R version of Python library that allows recombination event analysis through construction and exploration of similarity plots based on genetic distances between nucleotide sequences [@Babin2020]. The R version of `recan` is based on `Bio3D` and `Plotly` libraries. Package included two functions: `seqSim` and `scanSeqSim`. Both functions require a sequence alignment in fasta format as an input. Output is a dataframe with sequences similarities illustrated  by a single `plotly` plot or a list of plots for every pair of sequences with highlighted region of putative recombination.

# Statement of need
Nowadays the pandemic of SARS-COV-2 virus attracts attention to viral diseases and the evolution of its agents. Detection of recombinations is more reproducible to perform in an automatic pipeline, rather than use GUI [@Etherington2005; @Lole1999]. As R language [@R2018] is more widely used by biologists I ported `recan` package from Python to R.

#The user can change parameters for seqSim function: the sequence that is used for similarity calculation, the sliding window size, the window shift, and a region of interest (to plot only an area where breakpoints occur). 

#HCV_all_seq

#But if there is a lot of sequences in alignment it is difficult to find exact recombination. For nice looking picture it is required to exclude unnecessary lines.

#HCV_2_seq
Another package `seqcombo` provides the same analysis [@Yu2021]. It is important that data with genetic distances for further investigations are available in our package. Using `plotly` library for visualization data allows to exclude redundant plot lines from figure. A new function `scanSeqSim` for searching areas of recombination was added. In `scanSeqSim` you can adjust threshold value and recombination length the minimal number of points  in ‘cross-over’ expressed in the number of the sliding windows for reducing a number of false-positively detected recombination areas.  

# Testing and verification
To test the package, we used the same four sets of viral genomes from the article with a description of Python version of `recan` [@Babin2020].
The resulting `seqSim` method execution time with the window size = 400 and shift parameters = 200 showed in table 1. Both R packages are faster than Python version.

Table 1. Benchmarking of python and R version of recan package

|           | number of sequences |   bp   |   recan python   |      recan R     |     seqcombo R     |
|:---------:|:-------------------:|:------:|:----------------:|:----------------:|:------------------:|
|    HIV    |          25         |  3135  | 437 ms ± 7.74 ms | 152 ms ± 13.9 ms |  116 ms ± 16.6 ms  |
|    HCV    |          23         |  9431  | 579 ms ± 58.7 ms | 437 ms ± 63.5 ms | 132.8 ms ± 16.7 ms |
| Norovirus |          19         |  3366  | 648 ms ± 44.2 ms |  98 ms ± 6,42 ms |  94.1 ms ± 19.7 ms |
|    LSDV   |          14         | 150511 |  3.55 s ± 239 ms |  2.63 s ± 570 ms | 391.4 ms ± 40.2 ms |

Time execution test was performed using a laptop with 2 CPU cores and 8 Gb RAM.
With `scanSeqSim` function, `recan` can be used to identify and analyze recombination events in a large subset of sequences in separated plots. A set of SARS-CoV-2 full genome sequences was used to probe this function [@Paraskevis2020].
The distance plots with recombination events detected by `seqSim` and `scanSeqSim` are shown in Figures 1-2.
![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/lsdv_rec_sar.png)
_Figure 1. SARS-CoV-2 sequences._

![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/lsdv_rec_sar.png)
_Figure 2. Example of scanSecSim function for SARS-CoV-2 sequences._



# Availability and implementation
Recan is supported on Linux, MacOS and Windows. The package can be installed by `devtools` R package using command
```
install.packages("devtools")   
library(devtools)   
devtools::install_github("vasily-pavelko/recan")
```
The source code, guide and datasets are available on the GitHub repository (https://github.com/vasily-pavelko/recan). 


# Acknowledgement 
The author thanks Nataly Pavelko, Yuriy Babin and Alexander Sprygin for editing the manuscript.



# References
