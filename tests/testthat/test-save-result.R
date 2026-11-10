# Run after source('R/save_result.R') or devtools::load_all().
verify_save_result <- function() {
  root <- tempfile('save-result-test-')
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  t <- c(1830, 1840, 1850, 1860)
  base <- list(
    summary = list(baseline_timestep = 1830, final_timestep = 1860,
      excluded_species = 1L, known_species = 2L, singleton_species = 1L,
      known_extant = 1.5, known_extinct = .5, dark_extinct = -.1,
      total_extinct = .4, original_richness = 1.9,
      cumulative_extinction_rate = .4 / 1.9, AIC = 20),
    curve = data.frame(time = t, h = c(1, 2, 3, 4), m = c(0, 0, .2, 0),
      cumulative_extinction_rate = c(0, 0, .4 / 1.9, .4 / 1.9)),
    species = data.frame(species = c('A', 'B', 'C'), detected_timesteps = c(2, 1, 0),
      is_singleton = c(FALSE, TRUE, FALSE), d = c(1, .2, NA),
      p_extinct = c(.1, .4, NA), estimated_extinction_year = c(1850, 1850, NA)),
    group_results = NULL, do_k = NULL, bootstrap = NULL)
  read <- function(dir, name) utils::read.csv(file.path(dir, name), check.names = FALSE)
  fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
  a <- file.path(root, 'base')
  vis <- withVisible(save_result(base, a, prefix = ''))
  stopifnot(!vis$visible, nrow(vis$value) == 3L,
    identical(sort(list.files(a)), sort(c('summary.csv', 'time_result.csv', 'species_result.csv'))),
    read(a, 'summary.csv')$dark_extinct == -.1,
    identical(read(a, 'species_result.csv')$species, base$species$species),
    is.na(read(a, 'species_result.csv')$d[3]))
  before <- tools::md5sum(vis$value$file)
  fails(save_result(base, a, prefix = ''))
  stopifnot(identical(before, tools::md5sum(vis$value$file)))
  save_result(base, a, prefix = '', overwrite = TRUE)
  stopifnot(identical(before, tools::md5sum(vis$value$file)))

  grouped <- base
  gs <- as.data.frame(base$summary)
  gs$AIC <- NULL
  gs <- gs[rep(1, 3), , drop = FALSE]
  gs$group <- c('overall', '\u6e7f\u5730, "A"', 'excluded')
  gs$l <- c(1, 2, NA)
  gs$cumulative_extinction_rate[3] <- NA
  gc <- do.call(rbind, lapply(gs$group[1:2], function(g) {
    data.frame(group = g, time = rev(t), m = c(0, .2, 0, 0),
      cumulative_extinction_rate = c(.2, .2, 0, 0))
  }))
  grouped$group_results <- list(summary = gs, curve = gc)
  grouped$do_k <- data.frame(model = c('fixed', 'free'), k = c(0, -.2),
                             AIC = c(20, NA))
  bc <- data.frame(time = t, estimate = base$curve$cumulative_extinction_rate,
                   lower = c(0, 0, .1, NA), upper = c(0, 0, .3, NA),
                   n_valid = c(200L, 200L, 198L, 1L))
  bgc <- do.call(rbind, lapply(gs$group[1:2], function(g) {
    data.frame(group = g, bc, check.names = FALSE)
  }))
  grouped$bootstrap <- list(curve = bc, group_results = list(curve = bgc))
  b <- file.path(root, 'grouped')
  files <- save_result(grouped, b, prefix = 'hk_')
  stopifnot(nrow(files) == 5L, 'CI' %in% files$table,
    file.exists(file.path(b, 'hk_CI.csv')),
    !file.exists(file.path(b, 'hk_bootstrap_summary.csv')))
  ss <- read(b, 'hk_summary.csv')
  tt <- read(b, 'hk_time_result.csv')
  bb <- read(b, 'hk_CI.csv')
  stopifnot(nrow(ss) == 4L, ss$scope[1] == 'overall', is.na(ss$group[1]),
    identical(ss$group[-1], gs$group), all(is.na(ss$AIC[-1])),
    is.na(ss$cumulative_extinction_rate[4]), nrow(tt) == 12L,
    isTRUE(all.equal(tt$time[5:8], rev(t))), all(is.na(tt$h[-(1:4)])),
    nrow(bb) == 12L, identical(bb$n_valid[1:4], bc$n_valid),
    is.na(bb$lower[4]), isTRUE(all.equal(read(b, 'hk_k_comparison.csv'), grouped$do_k)))
  single_boot <- base
  single_boot$bootstrap <- list(curve = bc, group_results = NULL)
  stopifnot(nrow(save_result(single_boot, file.path(root, 'boot'))) == 4L)
  grouped$bootstrap <- NULL
  stopifnot(nrow(save_result(grouped, file.path(root, 'k'))) == 4L)
  grouped$do_k <- NULL
  stopifnot(nrow(save_result(grouped, file.path(root, 'groups'))) == 3L)
  one <- grouped
  one$group_results$summary <- gs[1, ]
  one$group_results$curve <- gc[gc$group == gs$group[1], ]
  stopifnot(nrow(save_result(one, file.path(root, 'one'))) == 3L)
  bad <- grouped
  bad$group_results$curve$time[1] <- 9999
  fails(save_result(bad, file.path(root, 'bad')))
  stopifnot(!dir.exists(file.path(root, 'bad')))
  fails(save_result(base, root, prefix = '../bad'))
  fails(save_result(base, root, overwrite = NA))
  fails(save_result(list(), root))
  # A conflict on a later file must prevent earlier files from being created.
  conflict <- file.path(root, 'conflict')
  dir.create(conflict)
  writeLines('unchanged', file.path(conflict, 'modgee_species_result.csv'))
  fails(save_result(base, conflict))
  stopifnot(identical(list.files(conflict), 'modgee_species_result.csv'),
    identical(readLines(file.path(conflict, 'modgee_species_result.csv')), 'unchanged'))
  message('save_result verification passed.')
  invisible(TRUE)
}

testthat::test_that("CSV exports preserve results and protect existing files", {
  testthat::expect_true(verify_save_result())
})
