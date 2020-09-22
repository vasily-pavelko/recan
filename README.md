
# recan
ported version of recan Python library https://github.com/babinyurii/recan for genetic distance plotting for recombination events analysis
## Requirements
***
R  
Biostrings  
plotly  
dplyr  


## Instalation
***
To install the package
```
install.packages("devtools")
library(devtools)
devtools::install_github("vasily-pavelko/recan")
```

## Usage example
***
Read multiple sequences alignment file
```
seq <- readBStringSet("hbv_C_Bj_Ba.fasta",  "fasta")
```
Inspect alignment   
```
seq
```

```{r}
BStringSet object of length 3:
    width seq                                                            names               
[1]  3215 TTCCACAGCATTCCACCAAGCTCTGCAGGA...GACACTCACCCTCAGGCCATGCAGTGGAA AB048704.1_genoty...
[2]  3215 CTCCACCACGTTCCACCAAACTCTTCAAGA...GACACTCATCCTCAGGCAATGCAGTGGAA AB033555.1_Ba
[3]  3215 CTCCACCACTTTCCACCAAACTCTTCAAGA...GACACTCATCCTCAGGCCGTGCAGTGGAA AB010291.1_Bj
```
We have three sequences in our alignment. BStringSet class is based upon the S4 class of the Biostrings library. Index corresponds to the sequence. 
After you've created the object you can draw the similarity plot. Call the method seqSim() of the BStringSet object to draw the plot. Pass the following parameters to the method:

window: sliding window size. The number of nucleotides the sliding window will span. It has the value of 200 by default.  
shift: this is the step our window slides downstream the alignment. It's value is set to 50 by default.  
region: the index of the potential recombinant. All the other sequences will be plotted as function of distance to that sequence.  
Load recan package and use seqSim() function to the seq object.  
The isolate of Ba genotype is the recombinant between the virus of C genotype and genotype Bj. Let's plot it. We set genotype Ba as the potential recombinant :
```
library(recan)
seqSim(seq)

```
![hbv_C_Bj_Ba](plots/hbv_C_Bj_Ba.png)

## Example datasets
***
To download the datasets use the following link: https://drive.google.com/drive/folders/1v2lg5yUDFw_fgSiulsA1uFeuzoGz0RjH?usp=sharing

## References
***
1. Recombination Analysis Tool (RAT): a program for the high-throughput detection of recombination. Bioinformatics, Volume 21, Issue 3, 1 February 2005, Pages 278–281, https://doi.org/10.1093/bioinformatics/bth500
https://sray.med.som.jhmi.edu/SCRoftware/simplot/

