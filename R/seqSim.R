
#'A seqSim function
#'
#'recan is a R package which allows to construct genetic distance plots to explore and discover recombination events in viral genomes.
#'This method has been previously implemented in desktop software tools: RAT, Simplot and RDP4.
#'Single function calculates sequence similarity along the alignment in a narrow window with customized step and save in dataframe.
#'This dataframe can be used for plots and further analysis.

#' @param seq as input we use fasta object from bio3d package. It is a multiple alignment of three or more sequences.
#' @param ref is a sequence which is used for sequences similarity calculation.
#' @param shift the step our window slides downstream the alignment. It's value is set to 200 by default
#' @param window sliding window size. The number of nucleotides the sliding window will span. It has the value of 50 by default.
#' @param region we can slice the alignment by using this parameter. region must be a a vector with two integers:
#' the start and the end position of the alignment slice. A sequence similarity plot will be built for specified region.
#' @examples
#' seqSim(seq, ref = 2, shift = 100, window = 100, region = c(100, 3000))
#' seqSim_data
#' @importFrom bio3d seqidentity
#' @importFrom plotly plot_ly
#' @importFrom dplyr %>%
#' @importFrom plotly layout
#' @importFrom plotly add_trace
#'
#' @export

seqSim <- function(seq, ref = 1, shift = 50, window = 200, region = c(0, 0)) {

  if (sum(region) == 0)
  {
    #steps for subsetting
    start_positions <- seq(from = 1, to = length(seq$ali[1,]), by = shift)
    #for full-length sequences
    start_positions_full <- start_positions[start_positions<length(seq$ali[1,])-window]
    #for sequences in the end, they are partial, less then window
    start_positions_part <- start_positions[start_positions>=length(seq$ali[1,])-window]
    start_positions_part <- start_positions_part[start_positions_part != length(seq$ali[1,])]
  } else {
    #steps for subsetting
    start_positions <- seq(from = region[1], to = region[2], by = shift)
    #for full-length sequences
    start_positions_full <- start_positions[start_positions<region[2]-window]
    #for sequences in the end, they are partial, less then window
    start_positions_part <- start_positions[start_positions>=region[2]-window]
    start_positions_part <- start_positions_part[start_positions_part != region[2]]
  }
  start_positions <- append(start_positions_full, start_positions_part)

  seq_sim <- list()
  for (i in start_positions_full) {
    seq_sim[[length(seq_sim)+1]] <- seqidentity(seq$ali[ ,i:(i+window)])
  }
  for (i in start_positions_part) {
    seq_sim[[length(seq_sim)+1]] <- seqidentity(seq$ali[ ,i:(length(seq$ali[1,]))])
  }
  seq_sim_vec <- lapply(seq_sim,"[", ref, )
  seq_sim_mat <- do.call(rbind, seq_sim_vec)
  seq_sim_mat[,ref] <- start_positions
  colnames(seq_sim_mat)[ref] <- "start_positions"
  assign("seqSim_data", seq_sim_mat, envir = .GlobalEnv)

  fig <- plot_ly()
  index <- c(1:ncol(seq_sim_mat))[c(1:ncol(seq_sim_mat)) != ref]
  for (l in index) {
    fig <- add_trace(fig, y=seq_sim_mat[,l], x=seq_sim_mat[,ref], name = colnames(seq_sim_mat)[l],
                     type = 'scatter', mode = 'lines')
  }
  fig <- fig %>% layout(xaxis = list(title = "Nucleotide position"),
                        yaxis = list (title = "Sequence identity"))
  fig
}

