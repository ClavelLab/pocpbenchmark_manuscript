## Analyses workflow for POCP benchmark manuscript

`popcbenchmark_manuscript` contains the workflow for analyzing data produced by our benchmark [`ClavelLab/pocpbenchmark`](https://github.com/ClavelLab/pocpbenchmark). We set out to compare proteins alignment tools for improved genus delineation using the Percentage Of Conserved Proteins (POCP).


### Setup the environment for the workflow

These analyses were conducted in R 4.3.1 and in Rstudio. We recommend setting up R and specific versions using [`rig`](https://github.com/r-lib/rig#id-features), and getting Rstudio from [Posit](https://posit.co/download/rstudio-desktop/). We also use [`renv`](https://rstudio.github.io/renv) for reproducible environment, which can be installed in R with `install.packages("renv")`.


1. Open Rstudio and create a new project via "File > New Project..."
2. Select "Version Control" and then "Git"
	1. Type `https://github.com/ClavelLab/pocpbenchmark_manuscript` in Repository URL.
	2. Make sure the project is going to be created in the correct subdirectory on your computer, or else edit accordingly
	3. Click on "Create project"

If you comfortable with the command line and git, clone the repository either with SSH or HTTPS in a suitable location.

3. Rstudio warns you that `One or more packages recorded in the lockfile are not installed` because a couple of R packages and dependencies are needed.
	1. Install the dependencies by typing `renv::restore()` in the Console and agree to the installation of the packages.
	2. Check that all dependencies are set by typing `renv::status()` in the Console where you should have `No issues found`


Our analysis workflow is orchestrated by [`targets`](https://docs.ropensci.org/targets/) and is composed of two subworkflows.

### Prepare the data for the analysis

> [!NOTE] I want to run the workflow!
> You can skip to the next part if you want to start the workflow from already prepared files!


1. Download the raw output files from the workflow using the "Download all" button: <https://doi.org/10.5281/zenodo.14974869>
2. Uncompress the zip archive within your project
3. Create a `data_benchmark` folder within your project.
4. Move all the zip files downloaded from zenodo (`benchmark-gtdb-f__*.zip`) to `data_benchmark`.
5. Ensure the two csv files are at the root of your project.
6. Run the workflow with the following command:

```r
Sys.setenv(TAR_PROJECT = "prepare_pocpbenchmark_data")
targets::tar_make()
```

### Analyze the data and build the manuscript

If you skipped the first workflow, you need to download the cleaned and formatted POCP/POCPu values and metadata tables for analysis from <https://doi.org/10.5281/zenodo.14975029>. These are the files you would have generated with the previous section.

1. Run the workflow with the following command:

```r
Sys.setenv(TAR_PROJECT = "analyze_pocpbenchmark_data")
targets::tar_make()
```

The manuscript is then available in the `_manuscript` folder, both as a HTML document (`index.html`) and a docx document. The figures are generated in the `figures` folder.
