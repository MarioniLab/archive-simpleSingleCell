# Workflow for low-level scRNA-seq data analyses

To run this workflow:

1. Run `add_figure.sh` to generate figure numbers.
2. Download the files containing the count data for each study from public repositories.
(This can be done by running the code block containing `all.urls <-` from `raw_workflow.Rmd`.
Alternatively, set `on.bioc <- TRUE` in `workflow.Rmd`, though this will also delete the files at the end of the workflow.)
3. Run `knitr::knit("workflow.Rmd")` in R.
4. Run `rmarkdown::render("workflow.md")` in R.
