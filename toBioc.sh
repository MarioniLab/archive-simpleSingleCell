cat workflow.Rmd | sed "s/on.bioc <- FALSE/on.bioc <- TRUE/" > ../bioc/vignettes/workflow.Rmd
cp ref.bib ../bioc/vignettes
