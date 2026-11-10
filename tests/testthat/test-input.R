test_that("invalid detection matrices are rejected", {
  M <- matrix(c(1, 1, 0, 1, 1, 0, 1, 1), nrow = 2, byrow = TRUE,
              dimnames = list(c("species_a", "species_b"),
                              c("1990", "2000", "2010", "2020")))
  bad <- M
  bad[1, 1] <- NA
  expect_error(fit_modgee(bad), "Missing values are not allowed", fixed = TRUE)
  bad <- M
  bad[1, 1] <- 2
  expect_error(fit_modgee(bad), "Detection must be 0 or 1", fixed = TRUE)
  bad <- M
  storage.mode(bad) <- "character"
  expect_error(fit_modgee(bad), "numeric or logical matrix", fixed = TRUE)
  expect_error(fit_modgee(M[, 1:3, drop = FALSE]), "4 to 215", fixed = TRUE)
  expect_error(fit_modgee(matrix(1, 2, 216)), "4 to 215", fixed = TRUE)
  singletons <- diag(4)
  dimnames(singletons) <- list(paste0("sp", 1:4), colnames(M))
  expect_error(fit_modgee(singletons), "At least two species", fixed = TRUE)
  expect_error(fit_modgee(M * 0), "At least two species", fixed = TRUE)
})

test_that("invalid time labels are rejected", {
  M <- matrix(1, 2, 4, dimnames = list(c("species_a", "species_b"), NULL))
  expect_error(fit_modgee(M), "Time column names are required", fixed = TRUE)
  colnames(M) <- c("2020", "2010", "2000", "1990")
  expect_error(fit_modgee(M), "strictly increasing", fixed = TRUE)
  colnames(M) <- c("1990", "2000", NA, "2020")
  expect_error(fit_modgee(M), "present and non-empty", fixed = TRUE)
  colnames(M) <- c("1990", "2000", "2000", "2020")
  expect_error(fit_modgee(M), "must be unique", fixed = TRUE)
  colnames(M) <- c("1990", "2001", "2010", "2020")
  expect_error(fit_modgee(M), "equally spaced", fixed = TRUE)
})

test_that("duplicate species names are rejected", {
  M <- matrix(1, 2, 4,
              dimnames = list(c("species_a", "species_a"),
                              c("1990", "2000", "2010", "2020")))
  expect_error(fit_modgee(M), "unique after trimming whitespace", fixed = TRUE)
})
