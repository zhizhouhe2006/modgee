test_that("public summaries preserve model totals and calendar years", {
  M <- rbind(zero = c(0,0,0,0,0), a = c(1,1,0,0,0),
    single = c(0,0,1,0,0), b = c(0,1,0,0,1), removed = c(0,0,0,0,0))
  times <- seq(1980, 2020, 10)
  colnames(M) <- times
  sp <- data.frame(d = c(NA,1,0.2,0.7,NA),
    p_extinct = c(NA,0.8,0.2,0,NA),
    expected_extinction_time = c(NA,2005,2010,NaN,NA))
  fit <- list(hs_est = rep(1,5), ms_est = c(0,0,0.2,0.1,0), LL = -10,
    optimizer_runs = rep(list(data.frame(convcode = 0)), 2))
  stats <- list(X_T_hat = -0.5, Mus_overall = c(0,0,0.1,0.2,0.2))
  x <- .modgee_public_output(M, times, fit, stats, sp, NULL)
  expected <- c("baseline_timestep", "final_timestep", "excluded_species",
    "known_species", "singleton_species", "known_extant", "known_extinct",
    "dark_extinct", "total_extinct", "original_richness",
    "cumulative_extinction_rate", "AIC")
  expect_identical(names(x$summary), expected)
  expect_equal(unlist(x$summary), setNames(c(1980,2020,2,3,1,2,1,-0.5,0.5,2.5,0.2,36), expected))
  expect_identical(names(x$species), c("species", "detected_timesteps", "is_singleton",
    "first_detection_timestep", "last_detection_timestep", "d", "p_extinct",
    "estimated_extinction_year"))
  expect_equal(x$species$first_detection_timestep, c(NA,1980,2000,1990,NA))
  expect_equal(x$species$last_detection_timestep, c(NA,1990,2000,2020,NA))
  expect_identical(x$species$is_singleton, c(FALSE,FALSE,TRUE,FALSE,FALSE))
  expect_true(is.na(x$species$estimated_extinction_year[4]))
  expect_null(x$group_results)
  expect_identical(names(x$curve), c("time", "h", "m", "cumulative_extinction_rate"))
  fit$optimizer_runs[[2]]$convcode <- 1
  expect_true(is.na(.modgee_public_output(M, times, fit, stats, sp, NULL)$summary$AIC))

  g <- list(levels = "A", species_group = setNames(c("A","A","A","A","B"), rownames(M)))
  stats$group_results <- list(
    summary = data.frame(group = "A", dark_extinction = -0.5, cumulative_extinction = 0.2),
    curve = data.frame(group = "A", time = times, m = fit$ms_est,
                       cumulative_extinction = stats$Mus_overall))
  fit$ls_est <- c(A = 0)
  y <- .modgee_public_output(M, times, fit, stats, sp, g)
  expect_identical(names(y$group_results$summary), c("group", "l", expected[-length(expected)]))
  expect_identical(names(y$species)[2], "belonged_group")
  s <- y$group_results$summary
  expect_equal(s$l, c(1, NA))
  expect_equal(s$excluded_species, c(1,1))
  expect_equal(s$known_species, c(3,0))
  expect_equal(s$known_extinct, c(1,0))
  expect_equal(s$original_richness, c(2.5,NA))
  expect_equal(sum(s$excluded_species), y$summary$excluded_species)
  expect_false("cumulative_extinction" %in% names(y$group_results$curve))
})
