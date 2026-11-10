# modgee 0.1.1

* Renamed the save_result bootstrap curve export to CI.csv and clarified its output documentation and help links.

* Added save_result() to export combined overall/group results, k comparisons and bootstrap curve intervals to CSV.

* Added bootstrap and CI arguments, full-pipeline row resampling, immutable
  per-replicate checkpoints, deterministic resume and cached k comparison.
* Added interval tables, compact replicate estimates and diagnostics.
* Added a developer script to generate example_bootstrap_result locally.
  The dataset is installed only after that script has completed.

# modgee 0.1.0

* Rename the fitting input to M and remove user-supplied starts from the public API.
* Supply example_file, example_M, example_groups, example_result and example_grouped_result.
* Resolve the example CSV path for the installed package and store precomputed fits.
* Update all function examples and cross-references to the new workflow.
