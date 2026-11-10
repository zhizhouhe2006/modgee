test_that("explicit reference overrides input ordering consistently", {
  sp <- letters[1:6]
  g <- c("A", "A", "B", "B", "C", "C")
  cnt <- rep(2L, 6)
  normal <- .normalize_modgee_groups(g, sp, cnt)
  expect_identical(.normalize_modgee_groups(g, sp, cnt, NULL), normal)
  ref <- .normalize_modgee_groups(g, sp, cnt, " B ")
  expect_identical(ref$levels, c("B", "A", "C"))
  expect_identical(ref$reference_group, "B")
  expect_identical(ref$H, normal$H[, c("B", "A", "C")])
  tab <- data.frame(species = rev(sp), group = rev(g))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::write.csv(tab, path, row.names = FALSE)
  forms <- list(tab, path, stats::setNames(g, sp),
                factor(g, levels = c("A", "B", "C")), normal$H)
  for (x in forms) {
    z <- .normalize_modgee_groups(x, sp, cnt, "B")
    expect_identical(z$H, ref$H)
    expect_identical(z$reference_group, "B")
  }
  expect_identical(.normalize_modgee_groups(c("A", "A"), sp[1:2], cnt[1:2], "A")$levels, "A")
})

test_that("invalid or inactive reference groups are rejected", {
  sp <- letters[1:4]; cnt <- rep(2L, 4); g <- c("A", "A", "B", "B")
  for (x in list(NA_character_, "", "  ", character(), c("A", "B"),
                 1, TRUE, factor("A"), matrix("A")))
    expect_error(.normalize_modgee_groups(g, sp, cnt, x), "single nonempty")
  expect_error(.normalize_modgee_groups(NULL, sp, cnt, "A"), "requires groups")
  expect_error(.normalize_modgee_groups(g, sp, cnt, "X"), "Available active groups: A, B")
  expect_error(.normalize_modgee_groups(g, sp, c(2, 2, 0, 0), "B"), "no recorded species")
  expect_error(.normalize_modgee_groups(g, sp, c(2, 2, 1, 1), "B"), "at least two")
})

test_that("group starts follow the final reference order", {
  lev <- c("B", "A", "C")
  expect_identical(.modgee_group_starts(NULL, lev), c(B = 0, A = 0, C = 0))
  expect_equal(.modgee_group_starts(c(A = -1, B = 0, C = 2), lev), c(B = 0, A = -1, C = 2))
  expect_equal(.modgee_group_starts(c(0, -1, 2), lev), c(B = 0, A = -1, C = 2))
  expect_error(.modgee_group_starts(c(A = 0, B = 1, C = 2), lev), "zero")
})
