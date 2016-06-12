# Workflow for low-level scRNA-seq data analyses

To run this code:

1. Run `add_figure.sh` to generate figure numbers.
2. Download the data files from public repositories.
(This can be done by running the code block containing `all.urls <-` from `raw_workflow.Rmd`.)
3. Run `knitr::knit("workflow.Rmd")` in R.
4. Run `rmarkdown::render("workflow.md")` in R.
