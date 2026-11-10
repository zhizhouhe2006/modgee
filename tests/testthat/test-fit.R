test_that("the standard fit matches the original reference", {

  ref <- readRDS(test_path("fixtures", "birds-reference.rds"))
  colnames(ref$M) <- as.character(ref$times)

  result <- .fit_modgee_with_diagnostics(
    ref$M,
    starts = ref$starts
  )

  actual <- list(
    d = as.numeric(result$parameters$d),
    h = as.numeric(result$parameters$h),
    m = as.numeric(result$parameters$m),
    LL = as.numeric(result$fit$LL),
    curve = as.numeric(result$curve$cumulative_extinction_rate),
    p_extinct = as.numeric(result$species$p_extinct),
    extinction_time =
      as.numeric(result$statistics$E_t_exts[match(result$species$species, result$fit_species_order)]),
    extinction_variance =
      as.numeric(result$statistics$var_t_exts[match(result$species$species, result$fit_species_order)]),
    dark_extinction = as.numeric(result$summary$dark_extinct)
  )

  for (name in names(ref$expected)) {
    expect_equal(
      actual[[name]],
      ref$expected[[name]],
      tolerance = 1e-8,
      info = name
    )
  }

  expect_true(result$status$completed)
  expect_true(result$status$converged)
  expect_identical(result$species$species, trimws(rownames(ref$M)))
})

test_that("all-zero species are retained without changing the fit", {

  ref <- readRDS(test_path("fixtures", "birds-reference.rds"))
  colnames(ref$M) <- as.character(ref$times)

  zero_name <- "test_zero_species"
  stopifnot(!zero_name %in% rownames(ref$M))

  # Put the zero row first to test both filtering and row mapping.
  M <- rbind(rep(0, ncol(ref$M)), ref$M)
  rownames(M)[1] <- zero_name

  starts <- ref$starts
  starts$ds_start <- c(NA_real_, starts$ds_start)

  result <- .fit_modgee_with_diagnostics(M, starts = starts)

  expect_identical(result$species$species, trimws(rownames(M)))
  expect_false((result$species$detected_timesteps > 0)[1])
  expect_true(is.na(result$species$d[1]))
  expect_true(is.na(result$species$p_extinct[1]))
  expect_equal(result$status$excluded_zero_species, 1)
  expect_equal(result$status$fitted_species, nrow(ref$M))

  expect_equal(
    as.numeric(result$species$p_extinct[-1]),
    ref$expected$p_extinct,
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(result$curve$cumulative_extinction_rate),
    ref$expected$curve,
    tolerance = 1e-8
  )
  expect_equal(
    result$summary$dark_extinct,
    ref$expected$dark_extinction,
    tolerance = 1e-8
  )

  expect_true(result$converged)
})

test_that("moving singletons preserves species-result correspondence", {

  ref <- readRDS(test_path("fixtures", "birds-reference.rds"))
  colnames(ref$M) <- as.character(ref$times)
  counts <- rowSums(ref$M)

  index <- c(which(counts == 1), which(counts >= 2))
  M <- ref$M[index, , drop = FALSE]

  starts <- ref$starts
  starts$ds_start <- starts$ds_start[index]

  result <- .fit_modgee_with_diagnostics(M, starts = starts)

  expect_identical(result$species$species, trimws(rownames(M)))

  original_order <- match(
    trimws(rownames(ref$M)),
    result$species$species
  )

  expect_equal(
    as.numeric(result$species$p_extinct[original_order]),
    ref$expected$p_extinct,
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(result$species$d[original_order]),
    ref$expected$d,
    tolerance = 1e-8
  )

  expect_true(result$converged)
})

test_that("automatic starting values are reproducible with a fixed seed", {

  ref <- readRDS(test_path("fixtures", "birds-reference.rds"))
  colnames(ref$M) <- as.character(ref$times)

  set.seed(123)
  first <- .fit_modgee_with_diagnostics(ref$M)

  set.seed(123)
  second <- .fit_modgee_with_diagnostics(ref$M)

  expect_equal(first$starts, second$starts, tolerance = 0)
  expect_equal(first$parameters, second$parameters, tolerance = 1e-8)
  expect_equal(first$curve, second$curve, tolerance = 1e-8)
  expect_equal(first$species, second$species, tolerance = 1e-8)
  expect_equal(
    first$summary$dark_extinct,
    second$summary$dark_extinct,
    tolerance = 1e-8
  )

  expect_true(first$converged)
  expect_true(second$converged)
})
