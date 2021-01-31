context("running scanSeqSim")

test_that("check dimentions", {
  load("./data/hcv_C_Bj_Ba_sequences.rda")
  scanSeqSim(seq)
  expect_equal(dim(data_scan), c(1, 6))
})

test_that("check calculations", {
  scanSeqSim(seq, ref = 2, window = 100, shift = 100, threshold = 0 , rec_detect = TRUE)
  expect_equal(data_scan[,6][[1]][1,1], 1801)
  expect_equal(data_scan[,6][[1]][2,1], 2101)
})
