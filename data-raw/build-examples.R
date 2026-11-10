# Run from the package root. Fits are generated once, not at package load time.
pkgload::load_all(".", quiet = TRUE)
example_M <- read_matrix(
  data = example_file,
  timestep = 10,
  baseline_min_species = 10
)
example_groups <- utils::read.csv(
  file.path("inst", "extdata", "birds_groups.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
set.seed(123)
example_result <- fit_modgee(M = example_M)
set.seed(123)
example_grouped_result <- fit_modgee(
  M = example_M,
  groups = example_groups,
  reference_group = "primary + old secondary",
  do_k = TRUE
)
stopifnot(isTRUE(example_result$converged),
          isTRUE(example_grouped_result$converged),
          !is.null(example_grouped_result$do_k))
dir.create("data", showWarnings = FALSE)
for (name in c("example_M", "example_groups", "example_result",
               "example_grouped_result")) {
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "xz", version = 2)
}
