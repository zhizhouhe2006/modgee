# Run once from the package root; interrupted runs reuse the same checkpoints.
(function() {
  package_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
  pkgload::load_all(package_dir, quiet = TRUE)
  work_dir <- file.path(dirname(package_dir), "tmp", "example_bootstrap")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)
  set.seed(123)
  example_bootstrap_result <- fit_modgee(
    M = example_M,
    groups = example_groups,
    reference_group = "primary + old secondary",
    bootstrap = 10
  )
  b <- example_bootstrap_result$bootstrap
  stopifnot(b$status$completed == 10L, b$status$pending == 0L,
            any(is.finite(b$curve$lower) & is.finite(b$curve$upper)))
  for (g in unique(b$group_results$curve$group)) {
    x <- b$group_results$curve[b$group_results$curve$group == g, ]
    stopifnot(any(is.finite(x$lower) & is.finite(x$upper)))
  }
  example_bootstrap_result$bootstrap$settings["checkpoint_dir"] <- list(NULL)
  example_bootstrap_result$bootstrap$settings$provenance <- list(
    seed = 123L, package_version = "0.1.1", R = R.version.string,
    dependencies = vapply(c("optimx", "Rcpp"), function(p)
      as.character(utils::packageVersion(p)), character(1)),
    generated = as.character(Sys.Date()))
  dir.create(file.path(package_dir, "data"), showWarnings = FALSE)
  save(example_bootstrap_result,
       file = file.path(package_dir, "data", "example_bootstrap_result.rda"),
       compress = "xz", version = 2)
  print(b$status)
})()
