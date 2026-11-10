diagnostic_matrix <- function(times = seq(1900, 1930, 10)) {
  matrix(1, nrow = 2, ncol = length(times),
         dimnames = list(c("A", "B"), as.character(times)))
}

diagnostic_error <- function(expr) {
  tryCatch({ force(expr); NA_character_ },
           error = function(e) conditionMessage(e))
}

test_that("all invalid time labels include their names and original positions", {
  x <- diagnostic_matrix(c("1900", "X1910", "Inf", "year"))
  error <- diagnostic_error(read_matrix(data = x, timestep = 10))
  expect_match(error, "must be finite numbers", fixed = TRUE)
  expect_match(error, 'time column 2 ("X1910")', fixed = TRUE)
  expect_match(error, 'time column 3 ("Inf")', fixed = TRUE)
  expect_match(error, 'time column 4 ("year")', fixed = TRUE)
})

test_that("empty and missing labels are located", {
  x <- diagnostic_matrix(c("1900", "", NA, " "))
  error <- diagnostic_error(read_matrix(data = x, timestep = 10))
  expect_match(error, "present and non-empty", fixed = TRUE)
  expect_match(error, 'time column 2 ("")', fixed = TRUE)
  expect_match(error, 'time column 3 (<NA>)', fixed = TRUE)
  expect_match(error, 'time column 4 (" ")', fixed = TRUE)
  colnames(x) <- NULL
  expect_error(read_matrix(data = x), "all column names are missing")
})

test_that("duplicate years list every original label and position", {
  x <- diagnostic_matrix(c("1910", "1900", "1900.0", "1910.0"))
  error <- diagnostic_error(read_matrix(data = x, timestep = 10))
  expect_match(error, "must be unique", fixed = TRUE)
  expect_match(error, 'year 1900: time column 2 ("1900"); time column 3 ("1900.0")',
               fixed = TRUE)
  expect_match(error, 'year 1910: time column 1 ("1910"); time column 4 ("1910.0")',
               fixed = TRUE)
})

test_that("CSV diagnostics use time-column positions excluding species", {
  x <- data.frame(species = c("A", "B"), a = c(1, 1), b = c(1, 1),
                  c = c(1, 1), check.names = FALSE)
  names(x) <- c("species", "1900", "X1910", "Inf")
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::write.csv(x, path, row.names = FALSE)
  error <- diagnostic_error(read_matrix(data = path, timestep = 10))
  expect_match(error, 'time column 2 ("X1910")', fixed = TRUE)
  expect_match(error, 'time column 3 ("Inf")', fixed = TRUE)
  names(x) <- c("species", "1900", "1900.0", "1910")
  utils::write.csv(x, path, row.names = FALSE)
  error <- diagnostic_error(read_matrix(data = path, timestep = 10))
  expect_match(error, 'year 1900: time column 1 ("1900"); time column 2 ("1900.0")',
               fixed = TRUE)
})

test_that("partly out-of-range requests retain the overlap with a message", {
  x <- diagnostic_matrix(seq(1900, 2020, 10))
  old_options <- options(warn = 2)
  on.exit(options(old_options))
  expect_message(
    M <- read_matrix(data = x, timestep = 10, time_range = c(1800, 2020)),
    "Retained interval labels before baseline selection: [1900, 2020]",
    fixed = TRUE
  )
  expect_identical(colnames(M), colnames(x))
  expect_equal(attr(M, "preparation")$time_range, c(1800, 2020))
  expect_message(
    M <- read_matrix(data = x, timestep = 10, time_range = c(1900, 2050)),
    "extends beyond available interval labels", fixed = TRUE
  )
  expect_identical(colnames(M), colnames(x))
  expect_message(
    M <- read_matrix(data = x, timestep = 10, time_range = c(1800, 2050)),
    "Retained interval labels before baseline selection: [1900, 2020]",
    fixed = TRUE
  )
  expect_identical(colnames(M), colnames(x))
})

test_that("range messages report selected labels before baseline selection", {
  x <- diagnostic_matrix(seq(1900, 2020, 10))
  expect_message(
    M <- read_matrix(data = x, timestep = 10,
                     time_range = c(1800, 2005), baseline = 1950),
    "Retained interval labels before baseline selection: [1900, 2000]",
    fixed = TRUE
  )
  expect_identical(colnames(M), as.character(seq(1950, 2000, 10)))
  expect_equal(attr(M, "preparation")$baseline, 1950)
})

test_that("in-range requests do not produce a clipping message", {
  x <- diagnostic_matrix(seq(1900, 2020, 10))
  messages <- character()
  M <- withCallingHandlers(
    read_matrix(data = x, timestep = 10, time_range = c(1910, 2000)),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_length(messages, 0L)
  expect_identical(colnames(M), as.character(seq(1910, 2000, 10)))
  expect_error(read_matrix(data = x, timestep = 10,
                           time_range = c(1800, 1890)), "No intervals fall")
})

test_that("CSV range clipping matches matrix input", {
  x <- diagnostic_matrix(seq(1900, 2020, 10))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  d <- data.frame(species = rownames(x), x, check.names = FALSE)
  utils::write.csv(d, path, row.names = FALSE)
  expect_message(
    M <- read_matrix(data = path, timestep = 10, time_range = c(1800, 2020)),
    "Retained interval labels before baseline selection: [1900, 2020]",
    fixed = TRUE
  )
  expect_identical(colnames(M), colnames(x))
  expect_identical(as.numeric(M), as.numeric(x))
})
