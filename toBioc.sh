cat workflow.Rmd | sed "s/on.bioc <- FALSE/on.bioc <- TRUE/" > ../bioc/workflow.Rmd
cp ref.bib ../bioc
