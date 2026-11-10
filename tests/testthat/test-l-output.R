test_that("group outputs report actual hazard multipliers", {
  fake_fit <- function(M, H, ...) {
    list(M = M, H = H, ds_est = rep(1, nrow(M)),
         hs_est = rep(1, ncol(M)), ms_est = c(0, 0, 0.1, 0.2),
         ls_est = c(0, log(2), log(0.5)), LL = 1,
         optimizer_runs = list(data.frame(convcode = 0, value = 1,
                                          row.names = "BFGS")))
  }
  local_mocked_bindings(
    .est_modgee_grouped = fake_fit,
    est_detect_extinct_parms_ML = fake_fit,
    compute_P_extinct_i = function(...) list(rep(0, 4), 4, 0)
  )
  M <- matrix(1, 6, 4, dimnames = list(letters[1:6], 2000:2003))
  starts <- list(ds_start = rep(1, 6), hs_start = rep(1, 4),
                 ms_start = c(0, 0, 0.1, 0.2))
  result <- .fit_modgee_with_diagnostics(M, starts = starts, groups = rep(c("A", "B", "C"), each = 2),
                       reference_group = "B")
  expected <- c(B = 1, A = 2, C = 0.5)
  expect_equal(result$parameters$l, expected)
  s <- result$group_results$summary
  expect_true("l" %in% names(s))
  expect_false("hazard_multiplier" %in% names(s))
  expect_equal(result$statistics$group_results$summary$l, unname(expected))
  expect_equal(result$fit$ls_est, log(expected))
  expect_equal(result$starts$ls_start, c(B = 0, A = 0, C = 0))
  for (g in names(expected)) {
    curve <- subset(result$group_results$curve, group == g)
    expect_equal(curve$m, c(0, 0, 0.1, 0.2) * expected[[g]])
    expect_equal(curve$cumulative_extinction_rate, -expm1(-cumsum(curve$m)))
  }
  one <- .fit_modgee_with_diagnostics(M, starts = starts, groups = rep("A", 6))
  expect_equal(one$parameters$l, c(A = 1))
  expect_equal(one$parameters$l, c(A = 1))
  expect_false("baseline_richness" %in% names(s))
  expect_equal(result$summary$original_richness, 2 * sum(exp(0.3 * expected)))
  ungrouped <- .fit_modgee_with_diagnostics(M, starts = starts)
  expect_equal(ungrouped$summary$original_richness, 6 * exp(0.3))
})
