test_that("public results omit internal diagnostic objects", {
  fake_fit <- function(M, H, ...) {
    list(M = M, ds_est = rep(1, nrow(M)), hs_est = rep(1, ncol(M)),
      ms_est = c(0, 0, 0.1, 0), LL = -10,
      optimizer_runs = rep(list(data.frame(convcode = 0, value = 10,
                                          row.names = "BFGS")), 2))
  }
  fake_stats <- function(fit) {
    list(Mus_overall = -expm1(-cumsum(fit$ms_est)),
      p_exts = rep(0, nrow(fit$M)), E_t_exts = rep(NaN, nrow(fit$M)),
      var_t_exts = rep(NaN, nrow(fit$M)),
      X_T_hat = nrow(fit$M) * expm1(0.1))
  }
  local_mocked_bindings(est_detect_extinct_parms_ML = fake_fit,
                        compute_summary_stats = fake_stats)
  M <- matrix(1, 2, 4, dimnames = list(c("a", "b"), seq(1990, 2020, 10)))
  starts <- list(ds_start = c(1, 1), hs_start = rep(1, 4),
                 ms_start = c(0, 0, 0.1, 0))
  set.seed(123)
  result <- fit_modgee(M)
  removed <- c("fit_species_order", "reference_species", "fit", "statistics")
  expect_false(any(removed %in% names(result)))
  expect_true(all(c("summary", "parameters", "species", "convergence") %in% names(result)))
  set.seed(123)
  detailed <- .fit_modgee_with_diagnostics(M)
  detailed["bootstrap"] <- list(NULL)
  expect_identical(result, detailed[setdiff(names(detailed), removed)])
  expect_equal(result$summary$AIC, 32)
})
