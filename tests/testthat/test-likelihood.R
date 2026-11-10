test_that("no-extinction likelihood agrees with a hand calculation", {

  # Detection probability is 1 - exp(-log(2)) = 0.5.
  observed <- c(1, 0)
  detection <- rep(log(2), 2)

  actual <- compute_L_i_no_extinction_Cpp(
    observed, detection
  )

  # Detected first, missed second: 0.5 * 0.5.
  expect_equal(actual, 0.25, tolerance = 1e-12)
})

test_that("extinction likelihood agrees with a hand calculation", {

  observed <- c(1, 0)
  detection <- rep(log(2), 2)
  extinction <- c(0, log(2))

  actual <- compute_L_i_Cpp(
    observed, detection, extinction
  )

  # First detection, then extinction: 0.5 * 0.5.
  # First detection, then survival and non-detection: 0.5 * 0.5 * 0.5.
  expect_equal(actual, 0.375, tolerance = 1e-12)
})

test_that("zero extinction reduces to the no-extinction likelihood", {

  observed <- c(1, 0, 1, 0)
  detection <- c(0.2, 0.6, 0.9, 0.4)

  expected <- prod(ifelse(
    observed == 1,
    1 - exp(-detection),
    exp(-detection)
  ))

  expect_equal(
    compute_L_i_Cpp(observed, detection, rep(0, 4)),
    expected,
    tolerance = 1e-12
  )

  expect_equal(
    compute_L_i_no_extinction_Cpp(observed, detection),
    expected,
    tolerance = 1e-12
  )
})

test_that("probabilities of all detection histories sum to one", {

  histories <- as.matrix(expand.grid(
    t1 = 0:1,
    t2 = 0:1,
    t3 = 0:1,
    t4 = 0:1
  ))

  probabilities <- apply(histories, 1, function(observed) {
    compute_L_i_Cpp(
      observed,
      c(0.2, 0.6, 0.9, 0.4),
      c(0, 0.1, 0.3, 0.2)
    )
  })

  expect_true(all(is.finite(probabilities)))
  expect_true(all(probabilities >= 0 & probabilities <= 1))
  expect_equal(sum(probabilities), 1, tolerance = 1e-12)
})
