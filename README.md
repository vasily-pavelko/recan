![](https://zenodo.org/badge/293751662.svg)
# recan
ported version of recan Python library https://github.com/babinyurii/recan for genetic distance plotting for recombination events analysis
## Requirements

R  
bio3d  
combinat   
plotly  
dplyr  


## Instalation

To install the package
```
install.packages("devtools")
library(devtools)
devtools::install_github("vasily-pavelko/recan")
```

## Usage example

Read multiple sequences alignment file
```
seq <- read.fasta("hbv_C_Bj_Ba.fasta")
```
Inspect alignment   
```
seq$id
```

```{r}
[1] "AB048704.1_genotype_C_" "AB033555.1_Ba"          "AB010291.1_Bj"         
```
We have three sequences in our alignment. Fasta class is created with read.fasta function from the bio3d library. Index corresponds to the sequence. After you've created the object you can draw the similarity plot. Call the method `seqSim()` of the fasta object to draw the plot. Pass the following parameters to the method:

`ref`: determine which sequence will be the referent. Similarity of other sequences will be calculated comparing with it.   
`window`: sliding window size. The number of nucleotides the sliding window will span. It has the value of 200 by default.  
`shift`: this is the step our window slides downstream the alignment. It's value is set to 50 by default.  
`region`: the index of the potential recombinant. All the other sequences will be plotted as function of distance to that sequence.  
Load `recan` package and use `seqSim()` function to the seq object.  
The isolate of Ba genotype is the recombinant between the virus of C genotype and genotype Bj. Let's plot it. We set genotype Ba as the potential recombinant :
```
library(recan)
seqSim(seq, ref = 2)

```
![hbv_C_Bj_Ba](plots/HBV_1_rec_C_B.png)


Potential recombinant is not shown in the plot, as the distances are calculated relative to it. The higher is the distance function (i.e. the closer to 1), the closer is the sequence to the recombinant and vice versa.
We can see typical 'crossover' of the distances which is the indicator of the possible recombination event. The distance of one isolate 'drops down' whereas the distance of the other remains the same of even gets closer to the potential recombinant, this abrupt drop shows that recombination could take place.
The figure from the article is shown below. It's just turned upside down relative to our plot, and instead of distance drop we see distance rising. Here Bj 'goes away' from the genotype C, whereas Ba keeps the same distance
![](plots/hbv_C_Bj_Ba.png)

By default `seqSim()` method plots the whole alignment. But after initial exploration, we can take a closer look at a particular region by passing the region parameter to the simgen method. We can slice the alignment by using this parameter. region must be a tuple or a list with two integers: the start and the end position of the alignment slice.
```
seqSim(seq, ref = 2, region = c(1000, 2700))
```
![](plots/hbv_slice_1.png)

To customize the plot or just to export and store the data, use object `seqSim_data`. seqSim_data returns matrix object with sequences as samples, and distances at given points as features.
```
head(seqSim_data)
```
```
AB048704.1_genotype_C_ start_positions AB010291.1_Bj
[1,]                  0.876            1000         0.915
[2,]                  0.891            1050         0.935
[3,]                  0.920            1100         0.955
[4,]                  0.905            1150         0.955
[5,]                  0.896            1200         0.950
[6,]                  0.900            1250         0.955
```

Once you've returned the data, you can easily customize the plot by using your favourite plotting library:
```
fig <- plot_ly()
ref = 2
index <- c(1:ncol(seqSim_data))[c(1:ncol(seqSim_data)) != ref]
for (l in index) {
  fig <- add_trace(fig, 
  y=seqSim_data[,l], 
  x=seqSim_data[,ref], 
  name = colnames(seqSim_data)[l], 
  type = 'scatter', 
  mode = 'lines')
}
fig <- fig %>% layout(xaxis = list(title = "Nucleotide position"),
                      yaxis = list (title = "Sequence identity"))
add_segments(fig, x = 1700, xend = 1700, color = I("red"), y = 0.8, yend = 1, line = list(dash = "dash"), name = "putative recombination break points")%>%
add_segments(x = 2200, xend = 2200, color = I("red"), y = 0.8, yend = 1, 
               line = list(dash = "dash"), showlegend = F)
fig
```
![](plots/hbv_slice_1.png)

to save the distance data in excel or csv format use the method `write.csv`:
```
write.csv(seqSim_data, file = "hbv_distance_data")
```
If there are about 20 or 30 sequences in the input file and their names are long, legend element may hide the plot. So, to be able to analyze many sequences at once, it's better to use short concise sequence names instead of long ones. Like this:

![](plots/short_names.png)

To illustrate how typical breakpoints may look, here are shown some examples of previously described recombinations in the genomes of different viruses. The fasta alignments used are available at dataset folder.

Recombination in HIV genome [5]:
![](plots/hcv_2k_1b_rec.png)

Norovirus recombinant isolate [7]:
![](plots/norovirus_rec.png)

For alignments with huge number of sequences a function `scanSeqSim()` can be applied. This function takes multiple alignment in `bio3d` fasta format. You should define a reference sequence, against which all other sequence similarities will be calculated. Then, for the each pair of sequence similarities function builds the plot and find the possible recombination regions. The start and the end of each region save as a matrix.
You can show results in `scan_data` data.frame

```
scan_data
```
```
     number_seq1 name_seq1                number_seq2 name_seq2      
[1,] 1           "AB048704.1_genotype_C_" 3           "AB010291.1_Bj"
     plots_rec_lines region_recomb_bp
[1,] List,8          Numeric,2 
```
Print plot.
```
data_scan[,5]
```
![](plots/hbv_scanSecSim.png)
And recombination region matrix.
```
data_scan[,6]
```
```
$region_recomb_bp
     [,1]
[1,] 1851
[2,] 2051
```


## Example datasets

To download the datasets use the following link: https://drive.google.com/drive/folders/1v2lg5yUDFw_fgSiulsA1uFeuzoGz0RjH?usp=sharing

## References

1. Recombination Analysis Tool (RAT): a program for the high-throughput detection of recombination. Bioinformatics, Volume 21, Issue 3, 1 February 2005, Pages 278–281, https://doi.org/10.1093/bioinformatics/bth500
https://sray.med.som.jhmi.edu/SCRoftware/simplot/
2. https://sray.med.som.jhmi.edu/SCRoftware/simplot/
3. Sprygin A, Babin Y, Pestova Y, Kononova S, Wallace DB, Van Schalkwyk A, et al. (2018) Analysis and insights into recombination signals in lumpy skin disease virus recovered in the field. PLoS ONE 13(12): e0207480. https://doi.org/ 10.1371/journal.pone.0207480

