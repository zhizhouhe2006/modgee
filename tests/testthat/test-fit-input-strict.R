strict_fit_matrix <- function(times = seq(1900, 1930, 10)) {
  matrix(1, 2, length(times),
         dimnames = list(c("A", "B"), as.character(times)))
}

strict_fit_error <- function(expr) {
  tryCatch({ force(expr); NA_character_ },
           error = function(e) conditionMessage(e))
}

test_that("fit_modgee no longer accepts a times argument", {
  expect_identical(
    names(formals(fit_modgee)),
    c("M", "groups", "reference_group", "do_k", "bootstrap", "CI")
  )
  x <- strict_fit_matrix()
  # Argument matching fails before model execution; message text is localized.
  expect_error(fit_modgee(x, times = seq(1900, 1930, 10)))
  expect_error(fit_modgee(), "Provide a detection matrix")
})

test_that("row names are mandatory and checked after whitespace trimming", {
  x <- strict_fit_matrix()
  rownames(x) <- NULL
  expect_error(fit_modgee(x), "Species row names are required")
  rownames(x) <- c(" ", "B")
  expect_error(fit_modgee(x), "Invalid rows: 1", fixed = TRUE)
  rownames(x) <- c("A", NA)
  expect_error(fit_modgee(x), "Invalid rows: 2", fixed = TRUE)
  rownames(x) <- c(" A ", "A")
  expect_error(fit_modgee(x), "unique after trimming")
  rownames(x) <- c(" B ", "A")
  before <- x
  validated <- .validate_modgee_matrix(x)
  expect_identical(rownames(validated$M), c("B", "A"))
  expect_identical(x, before)
})

test_that("missing and invalid time labels are located", {
  x <- strict_fit_matrix()
  colnames(x) <- NULL
  expect_error(fit_modgee(x), "Time column names are required")
  colnames(x) <- c("1900", "", NA, "1930")
  error <- strict_fit_error(fit_modgee(x))
  expect_match(error, 'time column 2 ("")', fixed = TRUE)
  expect_match(error, 'time column 3 (<NA>)', fixed = TRUE)
  colnames(x) <- c("1900", "X1910", "Inf", "1930")
  error <- strict_fit_error(fit_modgee(x))
  expect_match(error, 'time column 2 ("X1910")', fixed = TRUE)
  expect_match(error, 'time column 3 ("Inf")', fixed = TRUE)
  colnames(x) <- c("1900", "1900.0", "1920", "1930")
  error <- strict_fit_error(fit_modgee(x))
  expect_match(error, 'year 1900: time column 1 ("1900"); time column 2 ("1900.0")',
               fixed = TRUE)
})

test_that("fit input rejects unsorted and unequal time intervals", {
  x <- strict_fit_matrix(c(1900, 1920, 1910, 1930))
  expect_error(fit_modgee(x), "strictly increasing")
  x <- strict_fit_matrix(c(1900, 1910, 1921, 1930))
  expect_error(fit_modgee(x), "equally spaced")
})

test_that("real binary matrices and logical matrices are accepted", {
  x <- strict_fit_matrix()
  storage.mode(x) <- "integer"
  numeric_input <- .validate_modgee_matrix(x)
  logical_input <- .validate_modgee_matrix(x == 1L)
  expect_identical(numeric_input, logical_input)
  expect_identical(typeof(numeric_input$M), "double")
  expect_identical(numeric_input$times, c(1900, 1910, 1920, 1930))
  expect_null(attr(x, "preparation"))
  expect_error(fit_modgee(x + 0i), "real numeric or logical")
  character_input <- x
  storage.mode(character_input) <- "character"
  expect_error(fit_modgee(character_input), "real numeric or logical")
  expect_error(fit_modgee(as.data.frame(x)), "real numeric or logical")
})

test_that("invalid detections identify the species and time", {
  x <- strict_fit_matrix()
  for (value in c(2, -1, 0.5, NA, NaN, Inf)) {
    bad <- x
    bad[2, 3] <- value
    expect_error(fit_modgee(bad), "species 'B', time '1920'", fixed = TRUE)
  }
})

test_that("read_matrix and fit use the same absolute time-spacing tolerance", {
  within <- strict_fit_matrix(c(0, 10, 20 + 5e-9, 30))
  outside <- strict_fit_matrix(c(0, 10, 20 + 2e-8, 30))
  # Avoid calendar-alignment constraints here: compare the shared spacing rule.
  expect_true(all(.matrix_near(diff(as.numeric(colnames(within))), 10)))
  expect_identical(.validate_modgee_matrix(within)$times,
                   as.numeric(colnames(within)))
  expect_false(all(.matrix_near(diff(as.numeric(colnames(outside))), 10)))
  expect_error(.validate_modgee_matrix(outside), "equally spaced")
})

test_that("read_matrix output passes the stricter fitting validation", {
  x <- matrix(1, 2, 47, dimnames = list(c("A", "B"), 1980:2026))
  M <- suppressMessages(read_matrix(data = x, timestep = 10))
  validated <- .validate_modgee_matrix(M)
  expect_identical(validated$M, M)
  expect_equal(validated$times, seq(1980, 2020, 10))
  expect_false(tail(attr(validated$M, "preparation")$coverage$complete, 1))
  attr(M, "preparation") <- NULL
  expect_identical(.validate_modgee_matrix(M)$M, M)
  # The removed public argument must be rejected before optimization.
  expect_error(fit_modgee(M, starts = list()), "unused argument")
})

test_that("all-zero rows remain available and fit data limits are unchanged", {
  x <- strict_fit_matrix()
  x <- rbind(Zero = c(0, 0, 0, 0), x)
  expect_identical(.validate_modgee_matrix(x)$M, x)
  expect_error(fit_modgee(x, starts = list()), "unused argument")
  expect_error(fit_modgee(x[, 1:3, drop = FALSE]), "4 to 215")
  too_long <- strict_fit_matrix(seq_len(216))
  expect_error(fit_modgee(too_long), "4 to 215")
  x[2:3, ] <- 0
  expect_error(fit_modgee(x), "At least two species")
})

test_that("starting values must be real numeric vectors", {
  x <- strict_fit_matrix()
  starts <- list(ds_start = c(1, 1), hs_start = rep(1, 4),
                 ms_start = c(0, 0, 0.1, 0))
  starts$ds_start <- c(1 + 0i, 1 + 0i)
  expect_error(.fit_modgee_impl(x, starts = starts), "Invalid starting-value vector")
  starts$ds_start <- matrix(c(1, 1), 2, 1)
  expect_error(.fit_modgee_impl(x, starts = starts), "Invalid starting-value vector")
})
