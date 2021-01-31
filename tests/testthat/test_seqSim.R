context("running seqSim")

test_that("check dimentions", {
  load("./data/hcv_C_Bj_Ba_sequences.rda")
  seqSim(seq)
  expect_equal(dim(seqSim_data), c(65, 3))
})

test_that("check calculations", {
  seqSim(seq)
  expect_equal(round(mean(seqSim_data), digits = 3), 534.266)
})


test_that("error if shift less 0", {

  expect_error( seqSim(seq, shift = 0) )

})

test_that("error if region out of range", {

  expect_error( seqSim(seq, region = c(10, 0)) )

})
