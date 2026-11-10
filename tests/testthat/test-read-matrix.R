make_detection_csv <- function(times, values, species = rownames(values)) {
  if (is.null(species)) species <- paste0("sp", seq_len(nrow(values)))
  x <- data.frame(species = species, values, check.names = FALSE)
  names(x) <- c("species", as.character(times))
  path <- tempfile(fileext = ".csv")
  utils::write.csv(x, path, row.names = FALSE)
  path
}

matrix_values <- function(x) {
  attr(x, "preparation") <- NULL
  x
}

test_that("the public input argument is named data", {
  expect_identical(names(formals(read_matrix))[1], "data")
  expect_error(read_matrix(), "'data'")
})

test_that("CSV and matrix inputs produce the same detection histories", {
  x <- matrix(0, 3, 20, dimnames = list(c("A", "B", "C"), 1900:1919))
  x[1, c(1, 2, 11)] <- 1
  x[2, 11:20] <- 1
  path <- make_detection_csv(1900:1919, x)
  on.exit(unlink(path))
  csv <- suppressMessages(read_matrix(data = path, timestep = 10))
  mat <- suppressMessages(read_matrix(data = x, timestep = 10))
  expect_identical(matrix_values(csv), matrix_values(mat))
  expect_identical(colnames(mat), c("1900", "1910"))
  expect_identical(rownames(mat), c("A", "B", "C"))
  expect_equal(as.numeric(mat[, 1]), c(1, 0, 0))
  expect_equal(as.numeric(mat[, 2]), c(1, 1, 0))
  expect_identical(attr(csv, "preparation")$input_type, "csv")
  expect_identical(attr(mat, "preparation")$input_type, "matrix")
})

test_that("integer and logical matrices are accepted without modifying input", {
  x <- matrix(c(1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L), 2,
              dimnames = list(c("A", "B"), 1900:1903))
  attr(x, "preparation") <- list(baseline = -999)
  attr(x, "custom") <- "original"
  original <- x
  M <- read_matrix(data = x)
  L <- read_matrix(data = x == 1L)
  expect_identical(x, original)
  expect_identical(typeof(M), "double")
  expect_identical(matrix_values(M), matrix_values(L))
  expect_null(attr(M, "custom"))
  expect_equal(attr(M, "preparation")$baseline, 1900)
})

test_that("grouped input can be preserved or coarsened but not split", {
  x <- rbind(A = c(1, 0, 0, 1), B = c(0, 1, 1, 0))
  colnames(x) <- c("1900", "1910", "1920", "1930")
  path <- make_detection_csv(colnames(x), x)
  on.exit(unlink(path))
  for (input in list(x, path)) {
    M <- read_matrix(data = input, timestep = 10)
    expect_identical(matrix_values(M), x)
    M20 <- suppressMessages(read_matrix(data = input, timestep = 20))
    expect_identical(colnames(M20), c("1900", "1920"))
    expect_true(all(M20 == 1))
    expect_error(read_matrix(data = input), "integer multiple")
    expect_error(read_matrix(data = input, timestep = 15), "integer multiple")
  }
})

test_that("range filtering precedes threshold baseline selection for both inputs", {
  x <- rbind(A = c(1, 1, 1, 1), B = c(1, 0, 1, 0), C = c(0, 0, 0, 0))
  colnames(x) <- c("1900", "1910", "1920", "1930")
  path <- make_detection_csv(colnames(x), x)
  on.exit(unlink(path))
  for (input in list(x, path)) {
    M <- suppressMessages(read_matrix(
      data = input, timestep = 10, time_range = c(1910, 1930),
      baseline_min_species = 2
    ))
    expect_identical(colnames(M), c("1920", "1930"))
    expect_identical(rownames(M), c("A", "B", "C"))
    expect_true(all(M["C", ] == 0))
    expect_equal(attr(M, "preparation")$baseline, 1920)
    M2 <- suppressMessages(read_matrix(
      data = input, timestep = 10, time_range = c(1910, 1930), baseline = 1930
    ))
    expect_identical(colnames(M2), "1930")
    expect_error(read_matrix(data = input, timestep = 10, baseline = 1910,
                             baseline_min_species = 2), "not both")
    expect_error(read_matrix(data = input, timestep = 10,
                             baseline_min_species = 4), "No interval meets")
    expect_error(read_matrix(data = input, timestep = 10,
                             time_range = c(1910, 1930), baseline = 1900),
                 "must match")
    expect_error(read_matrix(data = input, timestep = 10,
                             time_range = c(2000, 2010)), "No intervals")
  }
})

test_that("time sorting preserves species order for both input types", {
  x <- rbind(B = c(0, 1, 1, 0), A = c(1, 0, 0, 1))
  colnames(x) <- c("1902", "1900", "1903", "1901")
  path <- make_detection_csv(colnames(x), x)
  on.exit(unlink(path))
  for (input in list(x, path)) {
    expect_message(M <- read_matrix(data = input), "reordered")
    expect_identical(colnames(M), as.character(1900:1903))
    expect_equal(as.numeric(M["B", ]), c(1, 0, 0, 1))
    expect_identical(rownames(M), c("B", "A"))
  }
})

test_that("invalid matrix names and values are rejected", {
  x <- matrix(1, 2, 4, dimnames = list(c("A", "B"), 1900:1903))
  bad <- x; rownames(bad) <- NULL
  expect_error(read_matrix(data = bad), "row names")
  bad <- x; rownames(bad) <- c("A", "A")
  expect_error(read_matrix(data = bad), "Duplicate species")
  bad <- x; rownames(bad) <- c(" ", "B")
  expect_error(read_matrix(data = bad), "missing or empty")
  bad <- x; colnames(bad) <- NULL
  expect_error(read_matrix(data = bad), "Time column names")
  bad <- x; colnames(bad)[2] <- "1900"
  expect_error(read_matrix(data = bad), "must be unique")
  bad <- x; colnames(bad)[2] <- "year"
  expect_error(read_matrix(data = bad), "finite numbers")
  for (value in c(2, -1, NA, Inf, NaN)) {
    bad <- x; bad[1, 2] <- value
    expect_error(read_matrix(data = bad), "species 'A', time '1901'")
  }
  bad <- x; storage.mode(bad) <- "character"
  expect_error(read_matrix(data = bad), "numeric or logical")
  expect_error(read_matrix(data = x + 0i), "numeric or logical")
  expect_error(read_matrix(data = as.data.frame(x)), "CSV file path or")
  expect_error(read_matrix(data = x[, 1, drop = FALSE]), "at least two")
  expect_error(read_matrix(data = tempfile()), "existing CSV")
})

test_that("CSV validation rejects duplicate species and nonbinary detections", {
  x <- matrix(1, 2, 4)
  duplicate <- make_detection_csv(1900:1903, x, c("A", "A"))
  x[1, 2] <- 2
  invalid <- make_detection_csv(1900:1903, x, c("A", "B"))
  on.exit(unlink(c(duplicate, invalid)))
  expect_error(read_matrix(data = duplicate), "Duplicate species")
  expect_error(read_matrix(data = invalid), "species 'A', time '1901'")
})

test_that("invalid timing and control arguments are rejected", {
  x <- matrix(1, 2, 4, dimnames = list(c("A", "B"), 1900:1903))
  bad <- x; colnames(bad) <- c("1900", "1901", "1903", "1904")
  expect_error(read_matrix(data = bad), "equally spaced")
  bad <- x; colnames(bad) <- c("1905", "1915", "1925", "1935")
  expect_error(read_matrix(data = bad, timestep = 20), "not aligned")
  for (width in list(0, -1, NA_real_, Inf, "10", c(1, 2))) {
    expect_error(read_matrix(data = x, timestep = width), "positive number")
  }
  expect_error(read_matrix(data = x, time_range = c(1910, 1900)), "ordered")
  expect_error(read_matrix(data = x, baseline = NA_real_), "finite number")
  expect_error(read_matrix(data = x, baseline_min_species = 1.5),
               "positive integer")
})

test_that("annual input ending in 2026 retains the partial final decade", {
  x <- matrix(0, 2, 47, dimnames = list(c("A", "B"), 1980:2026))
  x[, c("1980", "1990", "2000", "2010")] <- 1
  x[1, "2026"] <- 1
  path <- make_detection_csv(colnames(x), x)
  on.exit(unlink(path))
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)
  for (input in list(x, path)) {
    expect_message(M <- read_matrix(data = input, timestep = 10),
                   "Incomplete boundary intervals retained: 2020")
    expect_identical(tail(colnames(M), 1), "2020")
    expect_equal(as.numeric(M[, "2020"]), c(1, 0))
    coverage <- tail(attr(M, "preparation")$coverage, 1)
    expect_equal(coverage$observed_start, 2020)
    expect_equal(coverage$observed_end_exclusive, 2027)
    expect_identical(coverage$complete, FALSE)
  }
})

test_that("both incomplete boundaries are retained and one species stays a matrix", {
  x <- matrix(1, 1, 10, dimnames = list("A", 1905:1914))
  M <- suppressMessages(read_matrix(data = x, timestep = 10))
  expect_true(is.matrix(M))
  expect_identical(dim(M), c(1L, 2L))
  coverage <- attr(M, "preparation")$coverage
  expect_equal(coverage$observed_start, c(1905, 1910))
  expect_equal(coverage$observed_end_exclusive, c(1910, 1915))
  expect_identical(coverage$complete, c(FALSE, FALSE))
})
