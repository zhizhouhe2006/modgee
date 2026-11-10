test_that("the public fitting API uses M and no starts parameter", {
  expect_identical(names(formals(fit_modgee)),
                   c("M", "groups", "reference_group", "do_k", "bootstrap", "CI"))
  expect_error(fit_modgee(M = example_M, starts = list()), "unused argument")
  expect_error(fit_modgee(detections = example_M), "unused argument")
})

test_that("installed example inputs agree with their CSV sources", {
  expect_type(example_file, "character")
  expect_length(example_file, 1L)
  expect_true(file.exists(example_file))
  rebuilt <- suppressMessages(read_matrix(example_file, timestep = 10,
                                          baseline_min_species = 10))
  expect_equal(example_M, rebuilt)
  csv <- system.file("extdata", "birds_groups.csv", package = "modgee",
                     mustWork = TRUE)
  expect_identical(example_groups,
    utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE))
})

test_that("stored fits are complete and plots build without refitting", {
  expect_identical(example_result$M, example_M)
  expect_identical(example_grouped_result$M, example_M)
  expect_true(example_result$converged)
  expect_true(example_grouped_result$converged)
  expect_null(example_result$groups)
  expect_identical(example_grouped_result$groups$reference_group,
                   "primary + old secondary")
  expect_false(is.null(example_grouped_result$do_k))
  for (result in list(example_result, example_grouped_result)) {
    for (type in c("curve", "probability", "time")) {
      expect_no_warning(p <- plot_modgee(result, type))
      expect_s3_class(p, "ggplot")
      expect_no_error(ggplot2::ggplot_build(p))
    }
  }
})
