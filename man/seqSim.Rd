% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqSim.R
\name{seqSim}
\alias{seqSim}
\title{A seqSim function}
\usage{
seqSim(seq, ref = 1, shift = 50, window = 200, region = c(0, 0))
}
\arguments{
\item{seq}{as input we use fasta object from bio3d package. It is a multiple alignment of three or more sequences.}

\item{ref}{is a sequence which is used for sequences similarity calculation.}

\item{shift}{the step our window slides downstream the alignment. It's value is set to 200 by default}

\item{window}{sliding window size. The number of nucleotides the sliding window will span. It has the value of 50 by default.}

\item{region}{we can slice the alignment by using this parameter. region must be a a vector with two integers:
the start and the end position of the alignment slice. A sequence similarity plot will be built for specified region.}
}
\description{
recan is a R package which allows to construct genetic distance plots to explore and discover recombination events in viral genomes.
This method has been previously implemented in desktop software tools: RAT, Simplot and RDP4.
Single function calculates sequence similarity along the alignment in a narrow window with customized step and save in dataframe.
This dataframe can be used for plots and further analysis.
}
\examples{
seqSim(seq, ref = 2, shift = 100, window = 100, region = c(100, 3000))
seqSim_data
}
