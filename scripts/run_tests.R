## Runs the test suite. From the repository root:  Rscript scripts/run_tests.R
res <- testthat::test_dir("tests/testthat", reporter = testthat::ProgressReporter$new(show_praise = FALSE),
                          stop_on_failure = TRUE)
