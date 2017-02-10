cat ../local/workflow.Rmd | sed "s/on.bioc <- FALSE/on.bioc <- TRUE/" > workflow.Rmd
cp ../local/ref.bib .
