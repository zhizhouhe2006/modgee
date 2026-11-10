test_that("group representations align by species and preserve explicit levels", {
  sp <- c("a", "b", "c", "d")
  cnt <- rep(2, 4)
  tab <- data.frame(species = c("d", "b", "c", "a"), group = c("Z", "Z", "A", "A"))
  a <- .normalize_modgee_groups(tab, sp, cnt)
  b <- .normalize_modgee_groups(c(a = "A", b = "Z", c = "A", d = "Z"), sp, cnt)
  expect_identical(a$H, b$H)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  write.csv(tab, path, row.names = FALSE)
  expect_identical(.normalize_modgee_groups(path, sp, cnt)$H, a$H)
  expect_identical(.normalize_modgee_groups(a$H[c(4, 2, 1, 3), ], sp, cnt)$H, a$H)
  f <- factor(c("A", "Z", "A", "Z"), levels = c("Z", "A"))
  expect_identical(.normalize_modgee_groups(f, sp, cnt)$reference_group, "Z")
  expect_error(.normalize_modgee_groups(tab[-1, ], sp, cnt), "Missing groups")
  expect_error(.normalize_modgee_groups(rbind(tab, tab[1, ]), sp, cnt), "unique")
  expect_error(.normalize_modgee_groups(matrix(1, 4, 2), sp, cnt), "exactly one")
  expect_error(.normalize_modgee_groups(c("A", "A", "A", "Z"), sp, cnt), "at least two")
  tab$species[4] <- " a "
  expect_identical(.normalize_modgee_groups(tab, sp, cnt)$H, a$H)
})

test_that("empty groups and group starts are handled explicitly", {
  sp <- c("a", "b", "c")
  f <- factor(c("A", "A", "B"), levels = c("B", "A", "C"))
  expect_warning(g <- .normalize_modgee_groups(f, sp, c(2, 2, 0)), "removed")
  expect_identical(g$levels, "A")
  expect_setequal(g$removed_groups, c("B", "C"))
  expect_identical(unname(g$H[, 1]), c(1, 1, 0))
  expect_equal(.modgee_group_starts(c(B = 2, A = 0), c("A", "B")), c(A = 0, B = 2))
  expect_error(.modgee_group_starts(c(1, 0), c("A", "B")), "zero")
})

test_that("group summaries obey richness and extinction identities", {
  M <- rbind(a = c(1,1,0,0,0), b = c(1,0,1,0,0),
             c = c(1,1,0,1,1), d = c(0,1,1,0,1))
  g <- .normalize_modgee_groups(c("A","A","B","B"), rownames(M), rowSums(M))
  fit <- list(M = M, H = g$H, ds_est = c(1,0.7,1.2,0.8),
              hs_est = rep(1,5), ms_est = c(0,0,0.1,0.2,0), ls_est = c(0,log(0.5)))
  x <- .summary_modgee_grouped(fit, 1:5, g)
  z <- x$group_results$summary
  expect_equal(x$X_T_hat, sum(z$dark_extinction))
  expect_false("baseline_richness" %in% names(z))
  richness <- z$surviving_species / exp(-c(0.3, 0.15))
  expect_equal(x$p_T_hat, weighted.mean(z$cumulative_extinction, richness))
  expect_equal(sum(richness), nrow(M) + x$X_T_hat)
  expect_equal(z$cumulative_extinction, -expm1(-c(0.3,0.15)))
  # Changing the reference preserves every species' hazards and predictions.
  fit$ms_est <- fit$ms_est * 0.5
  fit$ls_est <- c(log(2),0)
  y <- .summary_modgee_grouped(fit, 1:5, g)
  expect_equal(y$p_exts, x$p_exts)
  expect_equal(y$Mus_overall, x$Mus_overall)
})
