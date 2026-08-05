test_that("run_editr runs on the example data", {
  res <- run_editr(example = TRUE)
  expect_type(res, "list")
  expect_true(!is.null(res$editing.df))
  expect_true(!is.null(res$sangs))
  expect_true(!is.null(res$sangs.filt))
  expect_true(!is.null(res$null.m.params))
  expect_true(!is.null(res$base.info))
  expect_true(!is.null(res$guide.coord))
})

test_that("run_editr returns editing probabilities", {
  res <- run_editr(example = TRUE)
  editing.df <- res$editing.df
  expect_true(all(c("A.pval", "C.pval", "G.pval", "T.pval") %in%
                    colnames(editing.df)))
  expect_true(all(editing.df$A.pval >= 0 & editing.df$A.pval <= 1))
})

test_that("run_editr works on test files", {
  ab1 <- system.file("testfiles", "BE MAFB5.ab1", package = "editR")
  res <- run_editr(
    ab1_file = ab1,
    guide = "AGCCGGCTGGCTGCAGGCGT"
  )
  expect_true(!is.null(res$editing.df))
})

test_that("run_editr requires an ab1 file and guide", {
  expect_error(run_editr(), "Please provide an .ab1")
  expect_error(run_editr(ab1_file = "test.ab1"), "Please provide a guide RNA")
})
