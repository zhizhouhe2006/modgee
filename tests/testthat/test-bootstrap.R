# Deterministic test fits exercise persistence without expensive optimization.
bootstrap_test_setup <- function(env = parent.frame()) {
  withr::local_dir(withr::local_tempdir(.local_envir = env), .local_envir = env)
  control <- new.env(parent = emptyenv())
  control$original <- 0L
  control$replicate <- 0L
  control$interrupt <- Inf
  control$fail <- integer()
  control$nonconverged <- integer()
  base <- example_grouped_result
  base$do_k <- NULL
  base$converged <- TRUE
  testthat::local_mocked_bindings(
    .fit_modgee_impl = function(M, groups = NULL, reference_group = NULL,
                               do_k = FALSE, .capture = NULL, ...) {
      control$original <- control$original + 1L
      if (!is.null(.capture)) .capture(list(), TRUE)
      base
    },
    .bootstrap_fit_replicate = function(plan, M, info) {
      control$replicate <- control$replicate + 1L
      if (control$replicate == control$interrupt)
        stop(structure(list(message = "test interruption", call = NULL),
                       class = c("interrupt", "condition")))
      if (control$replicate %in% control$fail) stop("test numerical failure")
      .bootstrap_set_rng(plan$fit_rng)
      x <- base
      value <- stats::runif(1)
      x$summary$dark_extinct <- -10 * value
      x$curve$cumulative_extinction_rate <- seq(0, value, length.out = nrow(x$curve))
      x$converged <- !control$replicate %in% control$nonconverged
      x
    },
    .package = "modgee", .env = env)
  control
}
bootstrap_test_call <- function(B = 3L, CI = .95, ...) {
  suppressWarnings(suppressMessages(fit_modgee(example_M, example_groups,
    reference_group = "primary + old secondary", bootstrap = B, CI = CI, ...)))
}

test_that("bootstrap parameters are strict", {
  for (x in list(0, 1, 2, -1, 1.2, NA_real_, Inf, TRUE, "2", c(1, 2), matrix(2),
                 structure(2, class = "custom")))
    expect_error(fit_modgee(example_M, bootstrap = x), "bootstrap")
  for (x in list(0, 1, -1, NA_real_, Inf, TRUE, "0.95", c(.9, .95), matrix(.95)))
    expect_error(fit_modgee(example_M, CI = x), "CI")
})

test_that("completed runs start fresh and preserve old checkpoints", {
  ctl <- bootstrap_test_setup()
  set.seed(10)
  seed <- .Random.seed
  a <- bootstrap_test_call()
  expect_identical(.Random.seed, seed)
  path <- a$bootstrap$settings$checkpoint_dir
  expect_true(file.exists(file.path(path, "completed.rds")))
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  sums <- tools::md5sum(files)
  b <- bootstrap_test_call(CI = .9)
  expect_equal(ctl$original, 2L)
  expect_equal(ctl$replicate, 6L)
  expect_false(b$bootstrap$status$resumed)
  expect_false(identical(path, b$bootstrap$settings$checkpoint_dir))
  expect_identical(sums, tools::md5sum(files))
  c <- bootstrap_test_call(5)
  expect_equal(ctl$original, 3L)
  expect_equal(ctl$replicate, 11L)
  expect_equal(c$bootstrap$status$completed, 5L)
})

test_that("interruption resumes the exact sample and random starts", {
  ctl <- bootstrap_test_setup()
  set.seed(123)
  continuous <- bootstrap_test_call(5)
  initial_dir <- getwd()
  resume_dir <- file.path(initial_dir, "resume")
  dir.create(resume_dir)
  setwd(resume_dir)
  ctl$replicate <- 0L
  ctl$interrupt <- 4L
  set.seed(123)
  interrupted <- tryCatch(bootstrap_test_call(5), interrupt = function(e) e)
  expect_s3_class(interrupted, "interrupt")
  paths <- list.dirs("modgee_bootstrap", recursive = TRUE, full.names = TRUE)
  path <- paths[grepl("^run_[0-9]{8}$", basename(paths))]
  files <- list.files(file.path(path, "results"), full.names = TRUE)
  expect_length(files, 3L)
  sums <- tools::md5sum(files)
  plan_path <- file.path(path, "plans", "rep_00000004.rds")
  plan_sum <- tools::md5sum(plan_path)
  ctl$interrupt <- Inf
  set.seed(999)
  resumed <- bootstrap_test_call(5)
  expect_identical(sums, tools::md5sum(files))
  expect_identical(plan_sum, tools::md5sum(plan_path))
  expect_equal(resumed$bootstrap$replicates, continuous$bootstrap$replicates)
  expect_equal(resumed$bootstrap$summary, continuous$bootstrap$summary)
  expect_equal(ctl$replicate, 6L)
})

test_that("errors and nonconvergence are terminal and not resampled", {
  ctl <- bootstrap_test_setup()
  ctl$fail <- 2L
  ctl$nonconverged <- 3L
  a <- bootstrap_test_call()
  expect_equal(a$bootstrap$diagnostics$state, c("success", "failed", "nonconverged"))
  expect_equal(a$bootstrap$status$successful, 1L)
  expect_equal(a$bootstrap$status$failed, 1L)
  expect_equal(a$bootstrap$status$nonconverged, 1L)
  expect_true(all(is.na(a$bootstrap$curve$lower)))
  expect_true(all(a$bootstrap$curve$n_valid == 1L))
  path <- a$bootstrap$settings$checkpoint_dir
  expect_true(file.exists(file.path(path, "completed.rds")))
  # Simulate interruption before committing completion.
  unlink(file.path(path, "completed.rds"))
  b <- bootstrap_test_call()
  expect_equal(ctl$replicate, 3L)
  expect_equal(a$bootstrap$diagnostics, b$bootstrap$diagnostics)
})

test_that("corrupt results recover from plans and corrupt plans stop", {
  ctl <- bootstrap_test_setup()
  set.seed(123)
  a <- bootstrap_test_call()
  path <- a$bootstrap$settings$checkpoint_dir
  rp <- file.path(path, "results", "rep_00000002.rds")
  unlink(file.path(path, "completed.rds"))
  writeLines("interrupted write", rp)
  b <- bootstrap_test_call()
  expect_equal(ctl$replicate, 4L)
  expect_equal(a$bootstrap$replicates, b$bootstrap$replicates)
  expect_length(list.files(file.path(path, "results"), pattern = "corrupt"), 1L)
  unlink(file.path(path, "completed.rds"))
  writeLines("corrupt", file.path(path, "plans", "rep_00000002.rds"))
  expect_error(bootstrap_test_call(), "Invalid bootstrap checkpoint")
})

test_that("NULL bootstrap has no disk or RNG side effects beyond the old fit", {
  ctl <- bootstrap_test_setup()
  set.seed(123)
  before <- .Random.seed
  x <- fit_modgee(example_M)
  expect_true("bootstrap" %in% names(x))
  expect_null(x$bootstrap)
  expect_identical(.Random.seed, before)
  expect_false(dir.exists("modgee_bootstrap"))
})

test_that("do_k comparison runs once per new run", {
  ctl <- bootstrap_test_setup()
  comparisons <- 0L
  local_mocked_bindings(.modgee_compare_k = function(...) {
    comparisons <<- comparisons + 1L
    data.frame(model = c("fixed_k", "free_k"), k = c(0, 1),
               AIC = c(20, 21), converged = TRUE)
  })
  a <- bootstrap_test_call(do_k = FALSE)
  b <- bootstrap_test_call(do_k = TRUE)
  c <- bootstrap_test_call(do_k = TRUE)
  d <- bootstrap_test_call(do_k = FALSE)
  expect_null(a$do_k)
  expect_null(d$do_k)
  expect_equal(b$do_k, c$do_k)
  expect_equal(comparisons, 2L)
  expect_equal(ctl$original, 4L)
  expect_equal(ctl$replicate, 12L)
})

test_that("sample copies have unique IDs and a fixed group model", {
  M <- example_M
  info <- .normalize_modgee_groups(example_groups, rownames(M), rowSums(M),
                                   "primary + old secondary")
  seen <- NULL
  local_mocked_bindings(.fit_modgee_impl = function(M, groups, reference_group, do_k) {
    seen <<- list(M = M, groups = groups, reference = reference_group, do_k = do_k)
    example_grouped_result
  })
  idx <- rep(seq_len(nrow(M)), 2)
  plan <- list(indices = idx, fit_rng = .bootstrap_rng())
  .bootstrap_fit_replicate(plan, M, info)
  expect_equal(anyDuplicated(rownames(seen$M)), 0L)
  expect_identical(names(seen$groups), rownames(seen$M))
  expect_equal(as.character(seen$groups), unname(info$species_group[idx]))
  expect_identical(seen$reference, info$reference_group)
  expect_false(seen$do_k)
  plan$indices <- which(info$species_group == info$reference_group)
  expect_error(.bootstrap_fit_replicate(plan, M, info), "lacks two")
})

test_that("ungrouped collection and nonfinite metrics are handled separately", {
  fit <- example_result
  vals <- .bootstrap_values(fit)
  rec <- list(state = "success", error = NA_character_, warnings = character(),
              elapsed = 0, convergence = fit$convergence, values = vals)
  records <- list(rec, rec, rec)
  records[[1]]$values$summary$value[1] <- Inf
  records[[2]]$values$summary$value[2] <- NA_real_
  records[[3]]$values$summary$value[2] <- NA_real_
  x <- suppressWarnings(.bootstrap_collect(fit, records, .95, "unused", FALSE))
  expect_null(x$group_results)
  expect_null(x$replicates$group_curve)
  expect_equal(x$summary$n_valid[1:2], c(2L, 1L))
  expect_true(is.na(x$summary$lower[2]))
  expect_equal(x$summary$estimate, vals$summary$value)
})
