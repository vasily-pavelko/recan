
#' A ScanSeqSim function
#'
#' This function similar to seqSim function calculates sequence similarity for user-defined reference sequence.
#' For every pair of sequences similarity substraction are performed.
#' Threshold parameter reduces the difference between two sequences similarities.
#' For areas where 'cross-over' of two is occur potential recombination event is detected.
#' Function creates object data_scan with indexes and names of the first and the second sequence, similarity plot and putative regions of recombination.
#'
#' @param seq as input we use fasta object from bio3d package. It is a multiple sequence alignment of three or more sequences.
#' @param ref is a sequence which is used for sequences similarity calculation.
#' @param shift the step our window slides downstream the alignment. Its value is set to 200 by default
#' @param window sliding window size. The number of nucleotides the sliding window will span. It has the value of 50 by default.
#' @param threshold value which reduce number of false-positive detected recombination events.
#' @param rec_length length of recombination. Minimal length of 'cross-over'. Expressed in number of shifts.
#' @param rec_detect logical value for adding red line of recombination to the plots.
#' @examples
#' scanSeqSim(seq, ref = 2, shift = 50, window = 200, threshold = 0.05,
#' rec_length = 3, rec_detect = TRUE)
#' data_scan
#' #Prints all information for the first pair of sequences
#' data_scan[1, ]
#' #Prints all recombination plots
#' data_scan[ ,5]
#' #Print all regions of recombination
#' data_scan[ ,6]
#'
#' @importFrom bio3d seqidentity
#' @importFrom plotly plot_ly
#' @importFrom dplyr %>%
#' @importFrom plotly layout
#' @importFrom plotly add_trace
#' @importFrom combinat combn
#' @export

scanSeqSim <- function(seq, ref = 1, shift = 50, window = 200, threshold = 0.05,
                       rec_length = 3, rec_detect = TRUE) {
  #error output
  if (shift <= 0) {
    stop("shift parameter can't be a negative or zero")
  } else if (window <= 0) {
    stop("window parameter can't be a negative or zero")
  } else if (length(seq) < 3) {
    stop("alignment contains less than 3 sequences")
  }
  #create vector with initial position for subsetting in the length of the input sequence
  start_positions <- seq(from = 1, to = length(seq$ali[1,]), by = shift)
  #choose steps which are produced a full-length window
  start_positions_full <- start_positions[start_positions<length(seq$ali[1,])-window]
  #choose steps in the end of the sequences, which are produced only partial subsequence
  start_positions_part <- start_positions[start_positions >= length(seq$ali[1,])-window]
  start_positions_part <- start_positions_part[start_positions_part != length(seq$ali[1,])]
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
  rownames(seq_sim_mat) <- start_positions
  assign("seqSim_data", seq_sim_mat, envir = .GlobalEnv)

  number_of_seq <- c(1:length(seq$id))
  number_of_seq_without_ref <- number_of_seq[number_of_seq != ref]
  combinations <- as.matrix(combn(number_of_seq_without_ref, 2))

  detectRecombination <- function(pair) {
    seq_sim_1 <- seq_sim_mat[ ,pair[[1]]]
    names(seq_sim_1) <- start_positions
    seq_sim_2 <- seq_sim_mat[ ,pair[[2]]]
    names(seq_sim_2) <- start_positions

    if (sum(seq_sim_1) > sum(seq_sim_2)) {
      find_difference <- (seq_sim_1 - seq_sim_2)+threshold < 0
      rec_region <- with(rle(as.numeric(find_difference)), {
        ok <- values == 1 & lengths > rec_length
        ends <- cumsum(lengths)
        starts <- ends - lengths + 1
        data.frame(starts, ends)[ok, ]
      })
    } else {
      find_difference <- (seq_sim_2 - seq_sim_1)+threshold > 0
      rec_region <- with(rle(as.numeric(find_difference)), {
        ok <- values == 0 & lengths > rec_length
        ends <- cumsum(lengths)
        starts <- ends - lengths + 1
        data.frame(starts, ends)[ok, ]
      })
    }
  }

  region_recomb <- apply(combinations, MARGIN = 2, detectRecombination)
  #add bp position
  change_bp <- function(y) {
    if (dim(y)[1] > 0) {
      apply(y, MARGIN = c(1,2), function(x) {replace(x, 1, start_positions[x])})
    }
    else {0
    }
  }
  region_recomb_bp <- sapply(region_recomb, change_bp)


#build plots
  combine_plots <- function(pair) {
      fig <- plot_ly()
      fig <- add_trace(fig,
              y=seq_sim_mat[ ,pair[[1]]],
              x=start_positions,
              name = colnames(seq_sim_mat)[pair[[1]]],
              type = 'scatter',
              mode = 'lines')
      fig <-  add_trace(fig,
              y=seq_sim_mat[ ,pair[[2]]],
              x=start_positions,
              name = colnames(seq_sim_mat)[pair[[2]]],
              type = 'scatter',
              mode = 'lines')
  }

  list_scan_plots <- apply(combinations, MARGIN = 2, combine_plots)

  #add red lines
  counter <<- 0
  plot_list <- function(matrices) {
    counter <<- counter + 1
    plot_matrix(matrices)
  }
  plot_matrix <- function(matrix) {
    fig_rec_detect <- list_scan_plots[[counter]]
    for (i in 1:length(matrix)) {
      fig_rec_detect <- add_segments(fig_rec_detect,
                                     x = matrix[[i]],
                                     xend = matrix[[i]],
                                     color = I("red"),
                                     y = min(seqSim_data[ ,combinations[1, counter]],
                                             seqSim_data[ ,combinations[2, counter]]),
                                     yend = 1,
                                     line = list(dash = "dash"),
                                     showlegend = FALSE)
    }
    fig_rec_detect
  }
  if (is.null(dim(region_recomb_bp))) {
    plots_rec_lines <- lapply(region_recomb_bp, plot_list)
  } else if  (dim(region_recomb_bp)[2] == 1) {
    counter <<- 1
    plots_rec_lines <- plot_matrix(region_recomb_bp)
  }
  #combine matrix
  number_seq1 <- combinations[1,]
  seq_names <- colnames(seq_sim_mat)
  name_seq1 <- seq_names[number_seq1]
  number_seq2 <- combinations[2,]
  name_seq2 <- seq_names[number_seq2]
  if (length(number_seq1) == 1) {
    plots_rec_lines_list <- list(plots_rec_lines)
    region_recomb_bp_list <- list(region_recomb_bp)
  } else {
    plots_rec_lines_list <- plots_rec_lines
    region_recomb_bp_list <- region_recomb_bp
    }
  data_recomb <- data.frame()
  data_recomb <- cbind(number_seq1, name_seq1, number_seq2, name_seq2, plots_rec_lines_list, region_recomb_bp_list)
  assign("data_scan", data_recomb, envir = .GlobalEnv)
}




