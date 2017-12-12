for wf in $(ls *Rmd)
do
    cat $wf | sed "s/on.bioc <- FALSE/on.bioc <- TRUE/" > ../bioc/vignettes/$wf
done    
cp ref.bib ../bioc/vignettes
