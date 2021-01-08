
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
  #error output
  if  (region[1] > region[2]) {
    stop("The value of the first nucleotide position should be less than the second one")
  } else if (shift <= 0) {
    stop("shift  parameter can't be a negative or zero")
  } else if (window <= 0) {
    stop("window parameter can't be a negative or zero")
  } else if (ref > length(seq$id)) {
    stop("index of reference sequence is out of number of sequences in your alighment")
  } else if (length(seq) < 3) {
    stop("alighment containes less then 3 sequences")
  }
  #if else condition if full length sequence is used for plot
  if (sum(region) == 0)
  {
    #create vector with initial position for subsetting in the length of the input sequence
    start_positions <- seq(from = 1, to = length(seq$ali[1,]), by = shift)
    #choose steps which are produced a full length window
    start_positions_full <- start_positions[start_positions<length(seq$ali[1,])-window]
    #choose steps in the end of the sequences, which are produced only partial subsequence
    start_positions_part <- start_positions[start_positions>=length(seq$ali[1,])-window]
    start_positions_part <- start_positions_part[start_positions_part != length(seq$ali[1,])]
  } else {
    # create vector with initial position for subsetting in the length of the input sequence
    start_positions <- seq(from = region[1], to = region[2], by = shift)
    #choose steps which are produced a full length window
    start_positions_full <- start_positions[start_positions<region[2]-window]
    #choose steps closer to the end of the sequences, which are produced subsequence shorter then window
    start_positions_part <- start_positions[start_positions>=region[2]-window]
    #test, if last element of start_position matches the end of the sequences
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
  #check are arguments good enough to draw plot
  if (length(seqSim_data[1, ]) < 2) {
    stop("too little elements in the output dataframe")
  }
  fig

}

