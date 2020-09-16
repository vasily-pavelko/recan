
#'A seqSim function
#'
#'recan is a R package which allows to construct genetic distance plots to explore and discover recombination events in viral genomes.
#'This method has been previously implemented in desktop software tools: RAT, Simplot and RDP4.
#'Single function calculates sequence similarity along the alignment in a narrow window with customized step and save in dataframe.
#'This dataframe can be used for plots and further analysis.

#' @param seq as input we use XString object from BioString package. It is a multiple alignment of three or more sequences.
#' @param shift the step our window slides downstream the alignment. It's value is set to 200 by default
#' @param window sliding window size. The number of nucleotides the sliding window will span. It has the value of 50 by default.
#' @param region we can slice the alignment by using this parameter. region must be a a vector with two integers: the start and the end position of the alignment slice. A sequence similarity plot will be built for specified region.
#' @examples
#' seqSim(seq, shift = 100, window = 100, region = c(100, 3000))
#' seqSim_data
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings pid
#' @importFrom Biostrings subseq
#' @importFrom plotly plot_ly
#' @importFrom dplyr %>%
#' @importFrom plotly layout
#' @importFrom plotly add_trace
#'
#' @export

seqSim <- function(seq, shift = 50, window = 200, region = c(0, 0)) {
  HV <- seq

  if (sum(region) == 0) {
    HV2 <-HV
  } else {
    HV2 <- subseq(HV, start = region[1], end = region[2])
  }

  start_positions <- seq(1, length(HV2[[1]]), by = shift)
  seq_sim_df <- matrix(0, ncol = length(HV2)-1, nrow = length(start_positions))
  seq_sim_df <- data.frame(seq_sim_df)
  for (j in 2:length(HV2))
  {
    seq_sim <- c()
    seq_sim_names <-c()
    for (i in start_positions) {
      if (i+window <= length(HV2[[1]]))
      {
        subseq <- subseq(HV2, start = i, end = (i+window))
        pwA <- pairwiseAlignment(subseq[[1]], subseq[[j]])
        seq_sim[length(seq_sim)+1] <- pid(pwA)/100
      }
      else
      {
        subseq <- subseq(HV2, start = i, end = length(HV2[[1]]))
        seq1 <- pairwiseAlignment(subseq[[1]], subseq[[j]])
        seq_sim_names[length(seq_sim)+1] <- i
        seq_sim[length(seq_sim)+1] <- pid(pwA)/100
      }
    }
    seq_sim_trans <- as.data.frame(seq_sim)
    seq_sim_df[ ,j-1] <- seq_sim_trans
    colnames(seq_sim_df) <- names(HV2)[2:length(HV2)]
  }


seq_sim_t_id <- cbind(start_positions, seq_sim_df)
data <- as.data.frame(seq_sim_t_id)
assign("seqSim_data", data, envir = .GlobalEnv)
fig <- plot_ly(data, x= ~data[,1], name = names(subseq)[2], y= ~data[,2],
               type =  'scatter', mode = 'lines')
fig <- fig %>% layout(xaxis = list(title = "Nucleotide position"),
                      yaxis = list (title = "Sequence identity"))
l <- 3
fig <- fig %>% add_trace(y = ~data[,l], name = names(subseq)[l], mode = 'lines')
fig
}

