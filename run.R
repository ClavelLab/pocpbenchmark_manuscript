#!/usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# Run the first part of the workflow
# https://books.ropensci.org/targets/projects.html#run-each-project
Sys.setenv(TAR_PROJECT = "prepare_pocpbenchmark_data")
targets::tar_make()
Sys.setenv(TAR_PROJECT = "analyze_pocpbenchmark_data")
targets::tar_make()
