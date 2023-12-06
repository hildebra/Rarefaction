library(rtk)
require(testthat)
context("Test rare")

test_that("rare works on columns and rows", {
  data <- matrix(c(10,0,5,2,3,4,0,4,5,2,4,10,3,0,0,2,3,5,10,0,0,5,19,20,0,4,5),3)
  a     <- rep(min(colSums(data)), ncol(data))
  b     <- rep(min(rowSums(data)), nrow(data))

  expect_equal(colSums(rtk(data, depth= min(a), ReturnMatrix = 1, margin = 2, verbose = F)$raremat[[1]]), a)
  expect_equal(rowSums(rtk(data, depth= min(b), ReturnMatrix = 1, margin = 1, verbose = F)$raremat[[1]]), b)
  #expect_equal(length("abc"), 3)
})


test_that("rare has the right names", {
  data <- matrix(seq(from = 1, to = 100), 10)
  cnames <- paste(rep("TestColNames"), 1:ncol(data))
  rnames <- paste(rep("TestRowNames"), 1:nrow(data))
  colnames(data) <- cnames
  rownames(data) <- rnames
  expect_equal(colnames(rtk(data, depth=min(colSums(data)), ReturnMatrix = 1, margin = 2, verbose=F)$raremat[[1]]), cnames)
  expect_equal(rownames(rtk(data, depth=min(colSums(data)), ReturnMatrix = 1, margin = 2, verbose=F)$raremat[[1]]), rnames)

  expect_equal(colnames(rtk(data, depth=min(rowSums(data)), ReturnMatrix = 1, margin = 1, verbose=F)$raremat[[1]]), cnames)
  expect_equal(rownames(rtk(data, depth=min(rowSums(data)), ReturnMatrix = 1, margin = 1, verbose=F)$raremat[[1]]), rnames)

  # Test rare has the right names when run with multiple threads
  # If column / row order is truly random, there is still a small chance this
  # test passes by randomly coming back in the right order. Could consider
  # running multiple times if this is a consistent issue.
  expect_equal(colnames(rtk(data, depth=min(colSums(data)), ReturnMatrix = 1, margin = 2, verbose=F, threads=4)$raremat[[1]]), cnames)
  expect_equal(rownames(rtk(data, depth=min(rowSums(data)), ReturnMatrix = 1, margin = 2, verbose=F, threads=4)$raremat[[1]]), rnames)

  expect_equal(colnames(rtk(data, depth=min(rowSums(data)), ReturnMatrix = 1, margin = 1, verbose=F, threads=4)$raremat[[1]]), cnames)
  expect_equal(rownames(rtk(data, depth=min(rowSums(data)), ReturnMatrix = 1, margin = 1, verbose=F, threads=4)$raremat[[1]]), rnames)

  # Test name order is consistent within rarefaction matrices from the same call
  rare_mats <- rtk(data, depth=min(colSums(data)), ReturnMatrix = 5, margin = 2, verbose=F, threads=4, repeats=5)$raremat
  for (i in seq(1, length(rare_mats) - 1)) {
    for (j in seq(i + 1, length(rare_mats))) {
      expect_equal(colnames(rare_mats[[i]]), colnames(rare_mats[[j]]))
      expect_equal(rownames(rare_mats[[i]]), rownames(rare_mats[[j]]))
    }
  }
  rare_mats <- rtk(data, depth=min(rowSums(data)), ReturnMatrix = 5, margin = 1, verbose=F, threads=4, repeats=5)$raremat
  for (i in seq(1, length(rare_mats) - 1)) {
    for (j in seq(i + 1, length(rare_mats))) {
      expect_equal(colnames(rare_mats[[i]]), colnames(rare_mats[[j]]))
      expect_equal(rownames(rare_mats[[i]]), rownames(rare_mats[[j]]))
    }
  }
})


test_that("divvs has right names", {
  # Single thread
  data <- matrix(seq(from = 1, to = 100), 10)
  cnames <- paste(rep("TestColNames"), 1:ncol(data))
  rnames <- paste(rep("TestRowNames"), 1:nrow(data))
  colnames(data) <- cnames
  rownames(data) <- rnames

  # Single thread cases
  res <- rtk(data, depth = min(colSums(data)), ReturnMatrix = 1, margin = 2,
             verbose = FALSE)
  expect_equal(
    lapply(res$divvs, \(x) {x$samplename}) |> unlist(),
    cnames
  )
  res <- rtk(data, depth = min(colSums(data)), ReturnMatrix = 1, margin = 1,
             verbose = FALSE)
  expect_equal(
    lapply(res$divvs, \(x) {x$samplename}) |> unlist(),
    rnames
  )

  # Multithread cases
  res <- rtk(data, depth = min(colSums(data)), ReturnMatrix = 1, margin = 2,
             verbose = FALSE, threads = 4)
  expect_equal(
    lapply(res$divvs, \(x) {x$samplename}) |> unlist(),
    cnames
  )
  res <- rtk(data, depth = min(colSums(data)), ReturnMatrix = 1, margin = 1,
             verbose = FALSE, threads = 4)
  expect_equal(
    lapply(res$divvs, \(x) {x$samplename}) |> unlist(),
    rnames
  )
})


test_that("zeros are reproduced", {
  data <- matrix(c(0,10,20,30,0, 0,0,10,20,30, 10,20,0,0,30, 10,20,30,0,0, 10,0,20,0,30, 10, 0,0,20,30), 5)
  nullpos <- which(data==0,arr.ind = T)
  samplesize <- min(colSums(data))
  data.r <- rtk(data, depth=samplesize, ReturnMatrix = 1, margin = 2, verbose=F)$raremat[[1]]
  nullpos.after <- which(data==0,arr.ind = T)
  # equal because we rarefy every column to its max, sow e actually do not even lose values!
  # the tables should also be the same
  # so this is an santiy check
  expect_equal(nullpos,nullpos.after)
  expect_equal(data.r,data)
})

test_that("Skipped samples are reported", {
    # columns that contain only zeros or are saller than depth
    # are skipped. this should be mentioned
    data  <- matrix(c(0,10,20,30,0, 0,0,0,2,3, 10,20,0,0,30, 10,20,30,0,0, 10,0,20,0,30, 10, 0,0,20,30), 5)
    data.r <- rtk(data, 1, 10)
    expect_that(length(data.r$skipped),equals(1)) 

})

test_that("Seed produces reproducible runs", {
    data  =  matrix(sample(0:100, 1000, replace = TRUE),
                            20)
    a = rtk(data, seed  = 42, depth = 5)
    b = rtk(data, seed  = 42, depth = 5)
    expect_equal(a,b)
})



