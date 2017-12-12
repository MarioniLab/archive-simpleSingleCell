---
title: "Analyzing single-cell RNA-seq data with Bioconductor (read counts)"
author: 
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Davis J. McCarthy
  affiliation: 
  - &EMBL EMBL European Bioinformatics Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SD, United Kingdom
  - St Vincent's Institute of Medical Research, 41 Victoria Parade, Fitzroy, Victoria 3065, Australia
- name: John C. Marioni
  affiliation: 
  - *CRUK
  - *EMBL
  - Wellcome Trust Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SA, United Kingdom
date: 10 December 2017
vignette: >
  %\VignetteIndexEntry{A worfklow for low-level analyses of single-cell RNA-seq data
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
    BiocStyle::html_document
bibliography: ref.bib
---







# Overview

We use a relatively simple dataset [@lun2017assessing] to introduce most of the concepts of scRNA-seq data analysis.
This dataset contains two plates of 416B cells (an immortalized mouse myeloid progenitor cell line), processed using the Smart-seq2 protocol [@picelli2014fulllength].
A constant amount of spike-in RNA from the External RNA Controls Consortium (ERCC) was also added to each cell's lysate prior to library preparation.
High-throughput sequencing was performed and the expression of each gene was quantified by counting the total number of reads mapped to its exonic regions.
Similarly, the quantity of each spike-in transcript was measured by counting the number of reads mapped to the spike-in reference sequences.
Counts for all genes/transcripts in each cell were obtained from ArrayExpress using the accession number [E-MTAB-5522](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5522).

# Setting up the data

## Loading in the count matrix

Our first task is to load the count matrices into memory.
One matrix was generated for each plate of cells used in the study.
In each matrix, each row represents an endogenous gene or a spike-in transcript, and each column represents a cell.
Subsequently, the count in each entry of the matrix represents the number of reads mapped to a particular gene/transcript in a particular cell.


```r
batch1 <- read.delim("counts_Calero_20160113.tsv", 
    header=TRUE, row.names=1, check.names=FALSE)
batch2 <- read.delim("counts_Calero_20160325.tsv", 
    header=TRUE, row.names=1, check.names=FALSE)
gene.lengths <- batch1$Length # First column is the gene length.
batch1 <- as.matrix(batch1[,-1]) # Discarding gene length (as it is not a cell).
batch2 <- as.matrix(batch2[,-1])
rbind(Batch1=dim(batch1), Batch2=dim(batch2))
```

```
##         [,1] [,2]
## Batch1 46703   96
## Batch2 46703   96
```

We combine the two matrices into a single object for further processing.
This is done after verifying that the genes are in the same order between the two matrices.


```r
stopifnot(identical(rownames(batch1), rownames(batch2)))
all.counts <- cbind(batch1, batch2)
```

For convenience, the count matrix is stored in a `SingleCellExperiment` object from the *[SingleCellExperiment](http://bioconductor.org/packages/SingleCellExperiment)* package.
This allows different types of row- and column-level metadata to be stored alongside the counts for synchronized manipulation throughout the workflow.


```r
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce
```

```
## class: SingleCellExperiment 
## dim: 46703 192 
## metadata(0):
## assays(1): counts
## rownames(46703): ENSMUSG00000102693 ENSMUSG00000064842 ... SIRV7 CBFB-MYH11-mcherry
## rowData names(1): GeneLength
## colnames(192): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1 SLX-9555.N701_S503.C89V9ANXX.s_1.r_1
##   ... SLX-11312.N712_S508.H5H5YBBXX.s_8.r_1 SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
## colData names(0):
## reducedDimNames(0):
## spikeNames(0):
```

We identify the rows corresponding to ERCC spike-in transcripts from the row names.
We store this information in the `SingleCellExperiment` object for future use.
This is necessary as spike-ins require special treatment in downstream steps such as normalization.


```r
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))
```

```
##    Mode   FALSE    TRUE 
## logical   46611      92
```

This dataset is slightly unusual in that it contains information from another set of spike-in transcripts, the Spike-In RNA Variants (SIRV) set.
For simplicity, we will only use the ERCC spike-ins in this analysis.
Thus, we must remove the rows corresponding to the SIRV transcripts prior to further analysis, which can be done simply by subsetting the `SingleCellExperiment` object.


```r
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)
```

```
##    Mode   FALSE    TRUE 
## logical   46696       7
```

**Comments from Aaron:**

- Some feature-counting tools will report mapping statistics in the count matrix (e.g., the number of unaligned or unassigned reads).
While these values can be useful for quality control, they would be misleading if treated as gene expression values.
Thus, they should be removed (or at least moved to the `rowData`) prior to further analyses.
- Be aware of using the `^ERCC` regular expression for human data where the row names of the count matrix are gene symbols.
An ERCC gene family actually exists in human annotation, so this would result in incorrect identification of genes as spike-in transcripts.
This problem can be avoided by publishing count matrices with standard identifiers (e.g., Ensembl, Entrez).

## Incorporating cell-based annotation

We load in the metadata for each library/cell from the `sdrf.txt` file.
It is important to check that the rows of the metadata table are in the same order as the columns of the count matrix.
Otherwise, incorrect metadata will be assigned to each cell.


```r
metadata <- read.delim("E-MTAB-5522.sdrf.txt", check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))
```

```
## [1] "Source Name"                "Comment[ENA_SAMPLE]"        "Comment[BioSD_SAMPLE]"     
## [4] "Characteristics[organism]"  "Characteristics[cell line]" "Characteristics[cell type]"
```

We only retain relevant metadata fields to avoid storing unnecessary information in the `colData` of the `SingleCellExperiment` object.
In particular, we keep the plate (i.e., `block`) in which each cell was processed, and the phenotype of each cell.
The second field is relevant as all of the cells contain a _CBFB-MYH11_ oncogene, but the expression of this oncogene is only induced in a subset of the cells.


```r
colData(sce)$Batch <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "WT")
colData(sce)$Phenotype <- pheno
table(colData(sce)$Phenotype, colData(sce)$Batch)
```

```
##          
##           20160113 20160325
##   induced       48       48
##   WT            48       48
```

## Incorporating gene-based annotation

Feature-counting tools typically report genes in terms of standard identifiers from Ensembl or Entrez.
These identifiers are used as they are unambiguous and highly stable.
However, they are difficult to interpret compared to the gene symbols which are more commonly used in the literature.
We can easily convert from one to the other using annotation packages like *[org.Mm.eg.db](http://bioconductor.org/packages/org.Mm.eg.db)*.


```r
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))
```

```
## DataFrame with 6 rows and 3 columns
##   GeneLength            ENSEMBL      SYMBOL
##    <integer>        <character> <character>
## 1       1070 ENSMUSG00000102693          NA
## 2        110 ENSMUSG00000064842          NA
## 3       6094 ENSMUSG00000051951        Xkr4
## 4        480 ENSMUSG00000102851          NA
## 5       2819 ENSMUSG00000103377          NA
## 6       2233 ENSMUSG00000104017          NA
```

It is often desirable to rename the row names of `sce` to the gene symbols, as these are easier to interpret.
However, this requires some work to account for missing and duplicate symbols.
The code below will replace missing symbols with the Ensembl identifier and concatenate duplicated symbols with the (unique) Ensembl identifiers.


```r
new.names <- rowData(sce)$SYMBOL
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ENSEMBL[missing.name]
dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
head(rownames(sce))
```

```
## [1] "ENSMUSG00000102693" "ENSMUSG00000064842" "Xkr4"               "ENSMUSG00000102851"
## [5] "ENSMUSG00000103377" "ENSMUSG00000104017"
```

We also determine the chromosomal location for each gene using the *[TxDb.Mmusculus.UCSC.mm10.ensGene](http://bioconductor.org/packages/TxDb.Mmusculus.UCSC.mm10.ensGene)* package.
This will be useful later as several quality control metrics will be computed from rows corresponding to mitochondrial genes.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   22428      13   24255
```

# Quality control on the cells 

## Defining the quality control metrics

Low-quality cells need to be removed to ensure that technical effects do not distort downstream analysis results.
We use several quality control (QC) metrics:

- The library size is defined as the total sum of counts across all features, i.e., genes and spike-in transcripts.
Cells with small library sizes are of low quality as the RNA has not been efficiently captured (i.e., converted into cDNA and amplified) during library preparation.
- The number of expressed features in each cell is defined as the number of features with non-zero counts for that cell.
Any cell with very few expressed genes is likely to be of poor quality as the diverse transcript population has not been successfully captured.
- The proportion of reads mapped to spike-in transcripts is calculated relative to the library size for each cell.
High proportions are indicative of poor-quality cells, where endogenous RNA has been lost during processing (e.g., due to cell lysis or RNA degradation).
The same amount of spike-in RNA to each cell, so an enrichment in spike-in counts is symptomatic of loss of endogenous RNA.
- In the absence of spike-in transcripts, the proportion of reads mapped to genes in the mitochondrial genome can also be used.
High proportions are indicative of poor-quality cells [@islam2014quantitative;@ilicic2016classification], possibly because of loss of cytoplasmic RNA from perforated cells.
The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.

For each cell, we calculate these quality control metrics using the `calculateQCMetrics` function from the *[scater](http://bioconductor.org/packages/scater)* package [@mccarthy2016scater].
These are stored in the row- and column-wise metadata of the `SingleCellExperiment` for future reference.


```r
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=isSpike(sce),
    Mt=which(rowData(sce)$CHR=="chrM")))
head(colnames(colData(sce)), 10)
```

```
##  [1] "Batch"                       "Phenotype"                   "total_features"             
##  [4] "log10_total_features"        "total_counts"                "log10_total_counts"         
##  [7] "pct_counts_top_50_features"  "pct_counts_top_100_features" "pct_counts_top_200_features"
## [10] "pct_counts_top_500_features"
```

The distributions of these metrics are shown in Figure \@ref(fig:qcplot416b).
The aim is to remove putative low-quality cells that have low library sizes, low numbers of expressed features, and high spike-in (or mitochondrial) proportions.


```r
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e6, xlab="Library size (millions)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
```

![Histograms of various QC metrics for all cells in the 416B data set. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes.](figure/qcplot416b-1.png)

It is also valuable to examine how the QC metrics behave with respect to each other (Figure \@ref(fig:qcbiplot416b)).
Generally, they will be in rough agreement, i.e., cells with low total counts will also have low numbers of expressed features and high ERCC/mitochondrial proportions.
Clear discrepancies may correspond to technical differences between batches of cells (see below) or genuine biological differences in RNA content.


```r
par(mfrow=c(1,3))
plot(sce$total_features, sce$total_counts/1e6, xlab="Number of expressed genes",
    ylab="Library size (millions)")
plot(sce$total_features, sce$pct_counts_ERCC, xlab="Number of expressed genes",
    ylab="ERCC proportion (%)")
plot(sce$total_features, sce$pct_counts_Mt, xlab="Number of expressed genes",
    ylab="Mitochondrial proportion (%)")
```

![Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the 416B dataset.](figure/qcbiplot416b-1.png)

## Removing low-quality cells based on outliers

### Identifying outliers for each metric 

Picking a threshold for these metrics is not straightforward as their absolute values depend on the experimental protocol.
For example, sequencing to greater depth will lead to more reads and more expressed features, regardless of the quality of the cells.
Similarly, using more spike-in RNA in the protocol will result in higher spike-in proportions.
To obtain an adaptive threshold, we assume that most of the dataset consists of high-quality cells, and identify cells that are outliers for the various QC metrics.

Outliers are defined based on the median absolute deviation (MADs) from the median value of each metric across all cells.
We remove cells with log-library sizes that are more than 3 MADs below the median log-library size.
A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.
We also remove cells where the log-transformed number of expressed genes is 3 MADs below the median value.


```r
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
```

We identify outliers for the proportion-based metrics in a similar manner.
Here, no transformation is required as we are identifying large outliers, for which the distinction should be fairly clear on the raw scale.
We do not use the mitochondrial proportions as we already have the spike-in proportions for this data set.


```r
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher")
```

Subsetting by column will retain only the high-quality cells that pass each filter described above.
We examine the number of cells removed by each filter as well as the total number of retained cells.
Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality.


```r
keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=sum(keep))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         4         0       1       188
```

At this point, we could apply `keep` to subset `sce` to only retain high-quality cells.
However, the above approach ignores prior information about the cells, which can be incorporated to improve the filtering procedure.
We describe why this is important and how it can be achieved in the following sections.

### Assumptions of outlier identification 

We have already mentioned the assumption that most cells are of high quality.
This is usually reasonable, and can be experimentally supported in some situations by visually checking that the cells are intact (e.g., on the microwell plate or in the microfluidics system).
Another assumption is that the QC metrics are independent on the biological state of each cell.
This ensures that any outlier values for these metrics are driven by technical factors rather than biological processes.
Thus, removing cells based on the metrics will not misrepresent the biology in downstream analyses.

The second assumption is most likely to be violated in highly heterogeneous cell populations.
For example, some cell types may naturally have less RNA or express fewer genes than other cell types.
Such cell types are more likely to be considered outliers and removed, even if they are of high quality.
The use of the MAD mitigates this problem by accounting for biological variability in the QC metrics.
A heterogeneous population should have higher variability in the metrics among high-quality cells, increasing the MAD and reducing the chance of incorrectly removing particular cell types (at the cost of reducing power to remove low-quality cells).
Nonetheless, filtering based on outliers may not be appropriate in extreme cases where one cell type is very different from the others.

We can explore the effect of known biological factors on the QC metrics using the `plotPhenoData` function.
Figure \@ref(qcplotbatch416b) demonstrates that induction of the CBFB-MYH11 oncogene results in some modest changes to the QC metric distributions.
This suggests that we could improve our QC step by considering the condition of the cell during outlier identification. 
Analyzing all conditions together would unnecessarily inflate the MAD and compromise the removal of low-quality cells. 


```r
multiplot(
    plotPhenoData(sce, aes_string(y="total_counts", x="Phenotype")),
    plotPhenoData(sce, aes_string(y="total_features", x="Phenotype")),
    plotPhenoData(sce, aes_string(y="pct_counts_ERCC", x="Phenotype")),
    plotPhenoData(sce, aes_string(y="pct_counts_Mt", x="Phenotype")),
    cols=2)
```

![Distribution of each QC metric for wild-type (WT) and oncogene-induced cells in the 416B data set.](figure/qcplotbatch416b-1.png)

### Blocking on known conditions

Systematic differences in the QC metrics can be handled to some extent using the `batch` argument in the `isOutlier` function.
Setting `batch` to the plate of origin will identify outliers within each level of `batch`, using plate-specific median and MAD estimates.
This is obviously useful for batch effects caused by known differences in experimental processing, e.g., sequencing at different depth or had different amounts of spike-in added.
We can also include known biological factors in `batch`, if those factors could result in systematically fewer expressed genes or lower RNA content.


```r
blocking <- paste0(sce$Batch, sce$Phenotype)
libsize.drop2 <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE, batch=blocking)
feature.drop2 <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE, batch=blocking)
spike.drop2 <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher", batch=blocking)
keep2 <- !(libsize.drop2 | feature.drop2 | spike.drop2)
data.frame(ByLibSize=sum(libsize.drop2), ByFeature=sum(feature.drop2),
    BySpike=sum(spike.drop2), Remaining=sum(keep2))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         5         4       6       183
```

The use of this blocking approach in `isOutlier` results in a small increase in the number of discarded cells.
This is expected given that the variability within each level of `batch` is lower, resulting in more power to detect outliers.
We then subset the `SingleCellExperiment` object to retain only the putative high-quality cells.


```r
sce <- sce[,keep2]
dim(sce)
```

```
## [1] 46696   183
```

## Alternative approaches to quality control

An alternative approach to quality control is to set pre-defined thresholds on each QC metric.
For example, we might remove all cells with library sizes below 100000 and numbers of expressed genes below 4000.
This generally requires a great deal of experience to determine appropriate thresholds for each experimental protocol and biological system.
Indeed, even with the same protocol and system, the appropriate threshold can vary from run to run due to the vagaries of RNA capture and sequencing.

Another strategy is to perform a principal components analysis (PCA) based on the quality metrics for each cell, e.g., the total number of reads, the total number of features and the proportion of mitochondrial or spike-in reads.
Outliers on a PCA plot may be indicative of low-quality cells that have aberrant technical properties compared to the (presumed) majority of high-quality cells.
In Figure \@ref(fig:pcaqualplothsc), no obvious outliers are present, which is consistent with the removal of suspect cells in the preceding quality control steps.


```r
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize
```

![PCA plot for cells in the 416Bdataset, constructed using quality metrics. The first and second components are shown on each axis, along with the percentage of total variance explained by each component. Bars represent the coordinates of the cells on each axis.](figure/pcaqualplot416b-1.png)

Methods like PCA-based outlier detection and support vector machines can provide more power to distinguish low-quality cells from high-quality counterparts [@ilicic2016classification].
This is because they are able to detect subtle patterns across many quality metrics simultaneously. 
However, this comes at some cost to interpretability, as the reason for removing a given cell may not always be obvious.
Thus, for this workflow, we will use the simple approach whereby each quality metric is considered separately.
Users interested in the more sophisticated approaches are referred to the *[scater](http://bioconductor.org/packages/scater)* and *[cellity](http://bioconductor.org/packages/cellity)* packages.

For completeness, we note that outliers can also be identified based on the gene expression profiles, rather than QC metrics.
However, we consider this to be a risky strategy as it can remove high-quality cells in rare populations.

# Classification of cell cycle phase 

We use the prediction method described by @scialdone2015computational to classify cells into cell cycle phases based on the gene expression data.
Using a training dataset, the sign of the difference in expression between two genes was computed for each pair of genes.
Pairs with changes in the sign across cell cycle phases were chosen as markers.
Cells in a test dataset can then be classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent with one phase or another.
This approach is implemented in the `cyclone` function from the *[scran](http://bioconductor.org/packages/scran)* package, using a pre-trained set of marker pairs for mouse data.


```r
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
```

The `cyclone` result for each cell in the HSC dataset is shown in Figure \@ref(fig:phaseplot416b).
Each cell is assigned a score for each phase, with a higher score corresponding to a higher probability that the cell is in that phase.
We focus on the G1 and G2/M scores as these are the most informative for classification.


```r
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)
```

![Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases.](figure/phaseplot416b-1.png)

Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M score; 
    in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score; 
    and in S phase if neither score is above 0.5.
Here, the vast majority of cells are classified as being in G1 phase.
We save these assignments into the `SingleCellExperiment` object for later use.


```r
sce$phases <- assignments$phases
table(sce$phases)
```

```
## 
##  G1 G2M   S 
##  99  62  22
```

Pre-trained classifiers are available in *[scran](http://bioconductor.org/packages/scran)* for human and mouse data. 
While the mouse classifier used here was trained on data from embryonic stem cells, it is still accurate for other cell types [@scialdone2015computational].
This may be due to the conservation of the transcriptional program associated with the cell cycle [@bertoli2013control;@conboy2007cell].
The pair-based method is also a non-parametric procedure that is robust to most technical differences between datasets.

__Comments from Aaron:__

- To remove confounding effects due to cell cycle phase, we can filter the cells to only retain those in a particular phase (usually G1) for downstream analysis.
Alternatively, if a non-negligible number of cells are in other phases, we can use the assigned phase as a blocking factor.
This protects against cell cycle effects without discarding information, and will be discussed later in more detail.
- The classifier may not be accurate for data that are substantially different from those used in the training set, e.g., due to the use of a different protocol.
In such cases, users can construct a custom classifier from their own training data using the `sandbag` function.
This will also be necessary for other model organisms where pre-trained classifiers are not available.
- Do not filter out low-abundance genes before applying `cyclone`.
Even if a gene is not expressed in *any* cell, it may still be useful for classification if it is phase-specific.
Its lack of expression relative to other genes will still yield informative pairs, and filtering them out would reduce power.

# Examining gene-level expression metrics

## Inspecting the most highly expressed genes

We also look at the identities of the most highly expressed genes (Figure \@ref(fig:topgene416b)).
This should generally be dominated by constitutively expressed transcripts, such as those for ribosomal or mitochondrial proteins.
The presence of other classes of features may be cause for concern if they are not consistent with expected biology.
For example, a top set containing many spike-in transcripts suggests that too much spike-in RNA was added during library preparation, while the absence of ribosomal proteins and/or the presence of their pseudogenes are indicative of suboptimal alignment.


```r
plotQC(sce, type = "highest-expression", n=50) + fontsize
```

![Percentage of total counts assigned to the top 50 most highly-abundant features in the 416B dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature.](figure/topgene416b-1.png)

## Filtering out low-abundance genes

Low-abundance genes are problematic as zero or near-zero counts do not contain much information for reliable statistical inference [@bourgon2010independent].
These genes typically do not provide enough evidence to reject the null hypothesis during testing, yet they still increase the severity of the multiple testing correction.
In addition, the discreteness of the counts may interfere with statistical procedures, e.g., by compromising the accuracy of continuous approximations.
Thus, low-abundance genes are often removed in many RNA-seq analysis pipelines before the application of downstream methods.

The "optimal" choice of filtering strategy depends on the downstream application.
A more aggressive filter is usually required to remove discreteness (e.g., for normalization) compared to that required for removing underpowered tests.
For hypothesis testing, the filter statistic should also be independent of the test statistic under the null hypothesis.
Thus, we (or the relevant function) will filter at each step as needed, rather than applying a single filter for the entire analysis.

Several metrics can be used to define low-abundance genes.
The most obvious is the average count for each gene, computed across all cells in the data set.
We calculate this using the `calcAverage` function, which also performs some adjustment for library size differences between cells 
We typically observe a peak of moderately expressed genes following a plateau of lowly expressed genes (Figure \@ref(fig:abhist416b)).


```r
ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))
```

![Histogram of log-average counts for all genes in the 416B dataset.](figure/abhist416b-1.png)

A minimum threshold can be applied to this value to filter out genes that are lowly expressed.
The example below demonstrates how we could remove genes with average counts less than 1.
The number of `TRUE` values in `demo.keep` corresponds to the number of retained rows/genes after filtering.


```r
demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
```

```
##    Mode   FALSE    TRUE 
## logical   33490   13206
```

We also examine the number of cells that express each gene.
This is closely related to the average count for most genes, as expression in many cells will result in a higher average (Figure \@ref(fig:nexprshist416b)).
Genes expressed in very few cells are often uninteresting as they are driven by amplification artifacts (though they may also also arise from rare populations).
We could then remove genes that are expressed in fewer than _n_ cells.


```r
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
    xlab=expression(Log[10]~"average count"))
```

![The number of cells expressing each gene in the 416B data set, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location.](figure/nexprshist416b-1.png)

As mentioned above, these filters will be applied at each step (automatically, in most cases, within the relevant function) rather than applied globally by subsetting `sce`.
This ensures that the most appropriate filter is used in each application.
Nonetheless, we remove genes that are not expressed in any cell to reduce computational work in downstream steps. 
Such genes provide no information and would be removed by any filtering strategy.


```r
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)
```

```
##    Mode   FALSE    TRUE 
## logical   22833   23863
```

# Normalization of cell-specific biases

## Using the deconvolution method to deal with zero counts

Read counts are subject to differences in capture efficiency and sequencing depth between cells [@stegle2015computational].
Normalization is required to eliminate these cell-specific biases prior to downstream quantitative analyses.
This is often done by assuming that most genes are not differentially expressed (DE) between cells.
Any systematic difference in count size across the non-DE majority of genes between two cells is assumed to represent bias and is removed by scaling.
More specifically, "size factors" are calculated that represent the extent to which counts should be scaled in each library.

Size factors can be computed with several different approaches, e.g., using the `estimateSizeFactorsFromMatrix` function in the *[DESeq2](http://bioconductor.org/packages/DESeq2)* package [@anders2010differential;@love2014moderated], or with the `calcNormFactors` function [@robinson2010scaling] in the *[edgeR](http://bioconductor.org/packages/edgeR)* package.
However, single-cell data can be problematic for these bulk data-based methods due to the dominance of low and zero counts.
To overcome this, we pool counts from many cells to increase the count size for accurate size factor estimation [@lun2016pooling].
Pool-based size factors are then "deconvolved" into cell-based factors for cell-specific normalization.


```r
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3432  0.7287  0.9224  1.0000  1.1485  3.5348
```

The size factors are well-correlated with the library sizes for all cells (Figure \@ref(fig:normplothsc)).
This suggests that most of the systematic differences between cells are driven by differences in capture efficiency or sequencing depth.
Any DE between cells would yield a non-linear trend between the total count and size factor, and/or increased scatter around the trend.
We observe some evidence of this after oncogene induction, where the size factors after induction are systematically lower.
This is consistent with composition biases introduced by upregulation of genes after induction.


```r
plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
    xlab="Library size (millions)", ylab="Size factor",
    col=c("red", "black")[sce$Phenotype], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
    legend=levels(sce$Phenotype))
```

![Size factors from deconvolution, plotted against library sizes for all cells in the 416B dataset. Axes are shown on a log-scale. Wild-type cells are shown in black and oncogene-induced cells are shown in red.](figure/normplot416b-1.png)

__Comments from Aaron:__

- While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero.
This manifests as negative size factors, which are obviously nonsensical.
To avoid this, the `computeSumFactors` function will automatically remove low-abundance genes prior to the calculation of size factors.
Genes with an average count below a specified threshold (`min.mean`) are ignored.
For read count data, the default value of 1 is usually satisfactory.
For UMI data, counts are lower so a threshold of 0.1 is recommended. 
- Cell-based QC should always be performed prior to normalization, to remove cells with very low numbers of expressed genes.
If this is not done, the `computeSumFactors` function may yield negative size factors for low-quality cells.
- The `sizes` argument can be used to specify the number of pool sizes to use to compute the size factors.
More `sizes` yields more precise estimates at the cost of some computational time and memory.
In general, `sizes` should not be below 20 cells, to ensure that there are sufficient non-zero expression values in each pool.
We also recommend that the total number of cells should be at least 100 for effective pooling.
- For highly heterogeneous data sets, it is advisable to perform a rough clustering of the cells.
This can be done with the `quickCluster` function and the results passed to `computeSumFactors` via the `cluster` argument.
Cells in each cluster are normalized separately, and the size factors are rescaled to be comparable across clusters.
This avoids the need to assume that most genes are non-DE across the entire population - only a non-DE majority is required between pairs of clusters.
We demonstrate this approach later with a larger dataset, as there are not enough cells in the 416B dataset.

## Computing separate size factors for spike-in transcripts

Size factors computed from the counts for endogenous genes are usually not appropriate for normalizing the counts for spike-in transcripts.
Consider an experiment without library quantification, i.e., the amount of cDNA from each library is _not_ equalized prior to pooling and multiplexed sequencing.
Here, cells containing more RNA have greater counts for endogenous genes and thus larger size factors to scale down those counts.
However, the same amount of spike-in RNA is added to each cell during library preparation.
This means that the counts for spike-in transcripts are not subject to the effects of RNA content.
Attempting to normalize the spike-in counts with the gene-based size factors will lead to over-normalization and incorrect quantification of expression.
Similar reasoning applies in cases where library quantification is performed. 
For a constant total amount of cDNA, any increases in endogenous RNA content will suppress the coverage of spike-in transcripts.
As a result, the bias in the spike-in counts will be opposite to that captured by the gene-based size factor.

To ensure normalization is performed correctly, we compute a separate set of size factors for the spike-in set.
For each cell, the spike-in-specific size factor is defined as the total count across all transcripts in the spike-in set.
This assumes that none of the spike-in transcripts are differentially expressed, which is reasonable given that the same amount and composition of spike-in RNA should have been added to each cell.
(See below for a more detailed discussion on spike-in normalization.)
These size factors are stored in a separate field of the `SingleCellExperiment` object by setting `general.use=FALSE` in `computeSpikeFactors`.
This ensures that they will only be used with the spike-in transcripts but not the endogenous genes.


```r
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
```

## Applying the size factors to normalize gene expression

The count data are used to compute normalized log-expression values for use in downstream analyses.
Each value is defined as the log~2~-ratio of each count to the size factor for the corresponding cell, after adding a prior count of 1 to avoid undefined values at zero counts.
Division by the size factor ensures that any cell-specific biases are removed.
If spike-in-specific size factors are present in `sce`, they will be automatically applied to normalize the spike-in transcripts separately from the endogenous genes. 


```r
sce <- normalize(sce)
```

The log-transformation is useful as it means that any differences in the values directly represent log~2~-fold changes in expression between cells.
This is usually more relevant than the absolute differences in coverage, which need to be interpreted in the context of the overall abundance.
The log-transformation also provides some measure of variance stabilization [@law2014voom], so that high-abundance genes with large variances do not dominate downstream analyses.
The computed values are stored as an `"logcounts"` matrix in addition to the other assay elements.



# Modelling the technical noise in gene expression

## Fitting a trend to the spike-in variances

Variability in the observed expression values across genes can be driven by genuine biological heterogeneity or uninteresting technical noise. 
To distinguish between these two possibiltiies, we need to model the technical component of the variance of the expression values for each gene.
We do so using the set of spike-in transcripts, which were added in the same quantity to each cell.
Thus, the spike-in transcripts should exhibit no biological variability, i.e., any variance in their counts should be technical in origin.

We use the `trendVar` function to fit a mean-dependent trend to the variances of the log-expression values for the spike-in transcripts.
We set `design` to block on the plate of origin for each cell, to ensure that technical differences between plates do not inflate the variances.
Given the mean abundance of a gene, the fitted value of the trend is then used as an estimate of the technical component for that gene.
The biological component of the variance is finally calculated by subtracting the technical component from the total variance of each gene with the `decomposeVar` function.


```r
design <- model.matrix(~sce$Batch)
var.fit <- trendVar(sce, parametric=TRUE, design=design, span=0.3)
var.out <- decomposeVar(sce, var.fit)
head(var.out)
```

```
##                           mean      total         bio       tech   p.value       FDR
## ENSMUSG00000103377 0.008029604 0.01179812 -0.02378567 0.03558380 1.0000000 1.0000000
## ENSMUSG00000103147 0.034571996 0.07202590 -0.08118252 0.15320842 1.0000000 1.0000000
## ENSMUSG00000103161 0.005210840 0.00496927 -0.01812296 0.02309223 1.0000000 1.0000000
## ENSMUSG00000102331 0.018575647 0.03262176 -0.04969763 0.08231939 1.0000000 1.0000000
## ENSMUSG00000102948 0.059116545 0.08826678 -0.17371266 0.26197944 1.0000000 1.0000000
## Rp1                0.097464502 0.45637197  0.02445072 0.43192124 0.2865203 0.8056285
```

We visually inspect the trend to confirm that it corresponds to the spike-in variances (Figure \@ref(fig:hvgplot416b))). 
The wave-like shape is typical of the mean-variance trend for log-expression values.
A linear increase in the variance is observed as the mean increases from zero, as larger variances are possible when the counts increase.
At very high abundances, the effect of sampling noise decreases due to the law of large numbers, resulting in a decrease in the variance.


```r
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
```

![Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red).](figure/hvgplot416b-1.png)

We check the distribution of expression values for the genes with the largest biological components.
This ensures that the variance estimate is not driven by one or two outlier cells (Figure \@ref(fig:hvgvioplot416b)).


```r
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize
```

![Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the 416B dataset. Each point represents the log-expression value in a single cell.](figure/hvgvioplot416b-1.png)

__Comments from Aaron:__

- In practice, trend fitting is complicated by the small number of spike-in transcripts and the uneven distribution of their abundances.
For low numbers of cells, these issues are exacerbated by the low precision of the variance estimates.
Some tuning of trend parameters such as `span` may be required to achieve a suitable fit - see `?trendVar` for more details.
Setting `parametric=TRUE` is especially useful for modelling the expected wave-like shape of the mean-variance relationship.
(This is not the default setting as it is not robust for arbitrary trend shapes.)
- The `trendVar` function will automatically filter out low-abundance genes prior to trend fitting.
This ensures that low-abundance genes do not interfere with the fit - either due to discreteness, which biases the estimate of variability of the variances around the trend;
or due to the frequency of low-abundance genes, which reduces the sensitivity of span-based smoothing algorithms at higher abundances.
    - Filtering uses the average of log-expression values rather than the (library size-adjusted) average count.
    The mean log-expression is independent of the variance estimate in a linear modelling framework [@bourgon2010independent], 
    which ensures that the filter does not introduce spurious trends in the variances at the filter boundary.
    - The filter threshold is specified with the `min.mean` argument in `trendVar`.
    We use the default threshold of 0.1 (`min.mean`) based on the appearance of discrete patterns in the variance estimates for simulated Poisson-distributed counts.
    Lower thresholds of 0.001-0.01 may be more suitable for very sparse data, e.g., from droplet-based protocols.
    - The filter used in `trendVar` is _not_ applied in `decomposeVar` by default.
    Retention of all genes ensures that weak biological signal from rare subpopulations is not discarded.
    To apply the filter in `decomposeVar`, users should set `subset.row=rowMeans(logcounts(sce)) > 0.1` in the function call.
- Negative biological components are often obtained from `decomposeVar`. 
These are intuitively meaningless as it is impossible for a gene to have total variance below technical noise.
Nonetheless, such values occur due to imprecise estimation of the total variance, especially for low numbers of cells.
- `decomposeVar` also yields _p_-values that can be used to define HVGs at a specific threshold for the false discovery rate (FDR).
We will discuss this in more detail later, as formal detection of HVGs is not necessary for feature selection during data exploration.

## Trend fitting when spike-ins are unavailable

If spike-in RNA has not been added in appropriate quantities (or at all), an alternative approach is to fit the trend to the variance estimates of the endogenous genes.
This is done using the `use.spikes=FALSE` setting in `trendVar`, as shown below.
The resulting trend is higher than most of the spike-in variances (Figure \@ref(fig:hvgplot416b2)), consistent with the presence of non-zero biological variation for most genes.


```r
var.fit.nospike <- trendVar(sce, parametric=TRUE, design=design,
    use.spikes=FALSE, span=0.2)
var.out.nospike <- decomposeVar(sce, var.fit.nospike)
plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)
```

![Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes (black), with spike-in transcripts shown in red.](figure/hvgplot416b2-1.png)

The output of `decomposeVar` is more difficult to interpret when the trend is fitted to the variances of endogenous genes.
The fitted value of the trend can no longer be generally interpreted as the technical component, as it contains some biological variation as well.
Instead, recall that `decomposeVar` quantifies changes in variances for each gene over the majority of genes with the same abundance.
One could assume that the variabilities of most genes are driven by constitutive "house-keeping" processes, which are generally uninteresting.
Any gene with an increase in its variance is _relatively_ highly variable and can be prioritized for further study.

# Denoising expression values using PCA

Once the technical noise is modelled, we can use principal components analysis (PCA) to remove random technical noise.
Consider that each cell represents a point in the high-dimensional expression space, where the spread of points represents the total variance.
PCA identifies axes in this space that capture as much of this variance as possible.
Each axis is a principal component (PC), where any early PC will explain more of the variance than a later PC.

We assume that biological processes involving co-regulated groups of genes will account for the most variance in the data.
If this is the case, this process should be represented by one or more of the earlier PCs.
In contrast, random technical noise affects each gene independently and will be represented by later PCs.
The `denoisePCA` function removes later PCs until the total discarded variance is equal to the sum of technical components for all genes used in the PCA.


```r
sce <- denoisePCA(sce, technical=var.fit$trend, design=design) 
dim(reducedDim(sce, "PCA")) 
```

```
## [1] 183  23
```

The function returns a `SingleCellExperiment` object containing the PC scores for each cell in the `reducedDims` slot.
The aim is to eliminate technical noise and enrich for biological signal in the retained PCs.
This improves resolution of the underlying biology during downstream procedures such as clustering.

__Comments from Aaron:__

- `denoisePCA` will internally filter to only genes that have positive biological components in `denoisePCA`.
This guarantees that the total technical variance to be discarded will not be greater than the total variance in the data.
- No filtering is performed on abundance here, which ensures that PCs corresponding to rare subpopulations can still be detected. 
Discreteness is less of an issue as low-abundance genes also have lower variance, thus reducing their contribution to the PCA.
- It is also possible to obtain a low-rank approximation of the original expression matrix, capturing the variance equivalent to the retained PCs.
This is useful for denoising prior to downstream procedures that require gene-wise expression values.


```r
sce2 <- denoisePCA(sce, technical=var.fit$trend, design=design, value="lowrank") 
assayNames(sce2)
```

```
## [1] "counts"    "logcounts" "lowrank"
```



# Data exploration with dimensionality reduction 

We visualize the relationships between cells by constructing pairwise PCA plots for the first three components (Figure \@ref(fig:pcaplothsc)).
Cells with similar expression profiles should be located close together in the plot, while dissimilar cells should be far apart.
In this case, we observe a clear separation of cells based on the oncogene induction status, consistent with the expected effects on the transcriptome.


```r
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, colour_by="Phenotype") + fontsize
```

![Pairwise PCA plots of the first three PCs in the 416B data set, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to the oncogene induction status. Bars represent the coordinates of the cells on each axis.](figure/pcaplot416b-1.png)

Another widely used approach is the _t_-stochastic neighbour embedding (_t_-SNE) method [@van2008visualizing].
_t_-SNE tends to work better than PCA for separating cells in more diverse populations.
This is because the former can directly capture non-linear relationships in high-dimensional space, whereas the latter must represent them on linear axes.
However, this improvement comes at the cost of more computational effort and requires the user to consider parameters such as the random seed and perplexity (see comments).
We demonstrate the generation of _t_-SNE plots in Figure \@ref(fig:tsneplot416b), using the low-rank approximation of the data to take advantage of the denoising step.


```r
out5 <- plotTSNE(sce, use_dimred="PCA", perplexity=5, colour_by="Phenotype", 
    rand_seed=100) + fontsize + ggtitle("Perplexity = 5")
out10 <- plotTSNE(sce, use_dimred="PCA", perplexity=10, colour_by="Phenotype",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 10")
out20 <- plotTSNE(sce, use_dimred="PCA", perplexity=20, colour_by="Phenotype",
    rand_seed=100) + fontsize + ggtitle("Perplexity = 20")
multiplot(out5, out10, out20, cols=3)
```

![_t_-SNE plots constructed from the denoised PCs in the 416B data set, using a range of perplexity values. Each point represents a cell, coloured according to its oncogene induction status. Bars represent the coordinates of the cells on each axis.](figure/tsneplot416b-1.png)

There are many other dimensionality reduction techniques that we do not consider here but could also be used, e.g., multidimensional scaling, diffusion maps.
These have their own advantages and disadvantages -- for example, diffusion maps (see `plotDiffusionMap`) place cells along a continuous trajectory and are suited for visualizing graduated processes like differentiation [@angerer2016destiny].

__Comments from Aaron:__

- For each visualization method, additional cell-specific information can be incorporated into the colour, size or shape of each point.
Here, cells are coloured by the total number of expressed features, which is strongly correlated with the first PC.
We will discuss this in more detail in the next section.
- For PCA, more components can be shown but these are usually less informative (and more difficult to interpret) as they explain less of the variance. 
- _t_-SNE is a stochastic method, so users should run the algorithm several times to ensure that the results are representative.
Scripts should set a seed (via the `rand_seed` argument) to ensure that the chosen results are reproducible.
It is also advisable to test different settings of the "perplexity" parameter as this will affect the distribution of points in the low-dimensional space.
A good guide on how to interpret _t_-SNE plots can be found at http://distill.pub/2016/misread-tsne/.

# Clustering cells into putative subpopulations

## Defining cell clusters from expression data

The denoised log-expression values are used to cluster cells into putative subpopulations.
Specifically, we perform hierarchical clustering on the Euclidean distances between cells, using Ward's criterion to minimize the total variance within each cluster.
This yields a dendrogram that groups together cells with similar expression patterns across the chosen genes.
An alternative approach is to cluster on a matrix of distances derived from correlations (e.g., as in `quickCluster`).
This is more robust to noise and normalization errors, but is also less sensitive to subtle changes in the expression profiles. 


```r
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")
```

Clusters are explicitly defined by applying a dynamic tree cut [@langfelder2008defining] to the dendrogram.
This exploits the shape of the branches in the dendrogram to refine the cluster definitions, and is more appropriate than `cutree` for complex dendrograms.
Greater control of the empirical clusters can be obtained by manually specifying `cutHeight` in `cutreeDynamic`.
Here, we also set `minClusterSize` to a lower value than the default of 20, to avoid spurious aggregation of distant small clusters.


```r
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
    minClusterSize=10, verbose=0))
```

We examine the distribution of cells in each cluster with respect to known factors.
Each cluster is comprised of cells from both batches, indicating that the clustering is not driven by a batch effect.
Differences in the composition of each cluster are observed with respect to `Phenotype`, consistent with a biological effect of oncogene induction.


```r
table(my.clusters, sce$Batch)
```

```
##            
## my.clusters 20160113 20160325
##           1       40       37
##           2       16       19
##           3       18       12
##           4       11       16
##           5        6        8
```

```r
table(my.clusters, sce$Phenotype)
```

```
##            
## my.clusters induced WT
##           1      77  0
##           2       0 35
##           3       0 30
##           4       2 25
##           5      14  0
```

We visualize the cluster assignments for all cells on the _t_-SNE plot in Figure \@(fig:tsnecluster416b).
Adjacent cells are generally assigned to the same cluster, indicating that the clustering procedure was applied correctly.


```r
sce$cluster <- factor(my.clusters)
plotTSNE(sce, use_dimred="PCA", colour_by="cluster", 
    perplexity=20, rand_seed=200) + fontsize
```

![_t_-SNE plot of the denoised PCs of the 416B data set. Each point represents a cell and is coloured according to the cluster identity to which it was assigned.](figure/tsnecluster416b-1.png)

We check the separatedness of the clusters using the silhouette width (Figure \@(fig:silhouette416b)).
Cells with large positive silhouette widths are closer to other cells in the _same_ cluster than to cells in _different_ clusters.
Conversely, cells with negative widths are closer to other clusters than to other cells in the cluster to which it was assigned.
Each cluster would ideally contain many cells with large positive widths, indicating that it is well-separated from other clusters.
This can be used to gauge the optimal parameter values (e.g., cut height, number of clusters) that maximize the separation between clusters.
For example, we could vary the cut height in `cutreeDynamic` to maximize the average silhouette width across all cells.


```r
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
    border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
```

![Barplot of silhouette widths for cells in each cluster. Each cluster is assigned a colour and cells with positive widths are coloured according to the colour of its assigned cluster. Any cell with a negative width is coloured according to the colour of the cluster that it is closest to. The average width for all cells in each cluster is shown, along with the average width for all cells in the data set.](figure/silhouette416b-1.png)



## Detecting marker genes between clusters

Once putative subpopulations are identified by clustering, we can identify marker genes for each cluster using the `findMarkers` function.
This fits a linear model to the log-expression values for each gene using *[limma](http://bioconductor.org/packages/limma)* [@ritchie2015limma].
The aim is to test for DE in each cluster compared to the others while blocking on uninteresting factors in `design`.
The top DE genes are likely to be good candidate markers as they can effectively distinguish between cells in different clusters.


```r
markers <- findMarkers(sce, my.clusters, design=design)
```

For each cluster, the DE results of the relevant comparisons are consolidated into a single output table.
This allows a set of marker genes to be easily defined by taking the top DE genes from each pairwise comparison between clusters.
For example, to construct a marker set for cluster 1 from the top 10 genes of each comparison, one would filter `marker.set` to retain rows with `Top` less than or equal to 10.
Other statistics are also reported for each gene, including the adjusted p-values (see below) and the log-fold changes relative to every other cluster.




```r
marker.set <- markers[["1"]]
head(marker.set, 10)
```

```
##    Top               Gene      FDR logFC.2 logFC.3 logFC.4 logFC.5
## 1    1             Pimreg 4.63e-47   -7.36   -6.18  -1.130  -7.046
## 2    1              Ccna2 1.10e-52   -7.40   -7.30  -2.495  -7.378
## 3    1              Myh11 8.01e-40    4.45    4.34   4.119   0.905
## 4    2              Pclaf 7.30e-21   -5.72   -7.70  -2.983  -5.471
## 5    2               Mcm5 2.66e-33   -6.16   -7.42  -6.582  -4.933
## 6    3              Aurkb 6.50e-30   -7.36   -6.88  -1.912  -6.466
## 7    3 ENSMUSG00000046057 4.63e-47   -3.88   -5.67  -1.797  -3.806
## 8    3              Kif11 8.15e-42   -6.19   -5.58  -0.627  -7.286
## 9    3               Mcm2 5.33e-30   -4.76   -6.35  -6.387  -3.905
## 10   4               Prc1 9.74e-40   -6.96   -6.01  -0.552  -6.914
```



We save the list of candidate marker genes for further examination.
The `overlapExprs` function may also be useful here, to prioritize candidates where there is clear separation between the distributions of expression values of different clusters.


```r
write.table(marker.set, file="marker_1.tsv", sep="\t", quote=FALSE, col.names=NA)
```

We visualize the expression profiles of the top candidates to verify that the DE signature is robust.
Figure \@ref(fig:heatmapmarker416b) indicates that most of the top markers have strong and consistent up- or downregulation in cells of cluster 1 compared to some or all of the other clusters.
Thus, cells from the subpopulation of interest can be identified as those that express the upregulated markers and do not express the downregulated markers.


```r
top.markers <- marker.set$Gene[marker.set$Top <= 10]
top.exprs <- logcounts(sce)[top.markers,,drop=FALSE]
heat.vals <- top.exprs - rowMeans(top.exprs)

library(pheatmap)
unique.clusters <- sort(unique(my.clusters))
chosen.cols <- clust.col[seq_along(unique.clusters)]
pheatmap(heat.vals, cluster_cols=my.tree,
    annotation_col=data.frame(Cluster=factor(my.clusters), Batch=factor(sce$Batch),
        Phenotype=sce$Phenotype, row.names=colnames(sce)),
    annotation_colors=list(Cluster=setNames(chosen.cols, unique.clusters),
        Batch=setNames(c("grey50", "grey80"), levels(sce$Batch)),
        Phenotype=setNames(c("red", "black"), levels(sce$Phenotype))))
```

![Heatmap of mean-centred normalized and corrected log-expression values for the top set of markers for cluster 1 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.](figure/heatmapmarker416b-1.png)

Many of the markers in Figure \@ref(fig:heatmapmarker416b) are not uniquely up- or downregulated in the chosen cluster.
Testing for unique DE tends to be too stringent as it overlooks important genes that are expressed in two or more clusters.
For example, in a mixed population of CD4^+^-only, CD8^+^-only, double-positive and double-negative T cells, neither _Cd4_ or _Cd8_ would be detected as subpopulation-specific markers because each gene is expressed in two subpopulations.
With our approach, both of these genes will be picked up as candidate markers as they will be DE between at least one pair of subpopulations.
A combination of markers can then be chosen to characterize a subpopulation, which is more flexible than trying to find uniquely DE genes.

__Comments from Aaron:__

- To avoid problems with discreteness when modelling the mean-variance relationship, `findMarkers` will automatically partition the data into low- and high-abundance genes.
Empirical Bayes shrinkage is performed in each partition separately, prior to calculation of _p_-values using the shrunk variance estimates.
This ensures that discreteness does not affect the inferences for high-abundance genes, without needing to entirely discard the low-abundance genes.
- `findMarkers` can also be directed to find genes that are DE between the chosen cluster and _all_ other clusters.
This should be done by setting `pval.type="all"`, which defines the p-value for each gene as the maximum value across all pairwise comparisons involving the chosen cluster.
Combined with `direction="up"`, this can be used to identify unique markers for each cluster.
However, this is sensitive to overclustering, as unique marker genes will no longer exist if a cluster is split into two smaller subclusters.
- It must be stressed that the (adjusted) _p_-values computed here cannot be properly interpreted as measures of significance.
This is because the clusters have been empirically identified from the data.
*[limma](http://bioconductor.org/packages/limma)* does not account for the uncertainty of clustering, which means that the _p_-values are much lower than they should be. 
This is not a concern in other analyses where the groups are pre-defined.

# Software availability

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R version 3.4.2 Patched (2017-10-30 r73642)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/R-3-4-branch_release/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/R-3-4-branch_release/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
##  [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.8                         cluster_2.0.6                         
##  [3] dynamicTreeCut_1.63-1                  scran_1.6.6                           
##  [5] scater_1.6.1                           ggplot2_2.2.1                         
##  [7] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0 GenomicFeatures_1.30.0                
##  [9] org.Mm.eg.db_3.5.0                     AnnotationDbi_1.40.0                  
## [11] SingleCellExperiment_1.0.0             SummarizedExperiment_1.8.0            
## [13] DelayedArray_0.4.1                     matrixStats_0.52.2                    
## [15] Biobase_2.38.0                         GenomicRanges_1.30.0                  
## [17] GenomeInfoDb_1.14.0                    IRanges_2.12.0                        
## [19] S4Vectors_0.16.0                       BiocGenerics_0.24.0                   
## [21] mvoutlier_2.0.8                        sgeostat_1.0-27                       
## [23] Rtsne_0.13                             BiocParallel_1.12.0                   
## [25] knitr_1.17                             BiocStyle_2.6.1                       
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.1          igraph_1.1.2             plyr_1.8.4              
##   [4] lazyeval_0.2.1           sp_1.2-5                 shinydashboard_0.6.1    
##   [7] splines_3.4.2            digest_0.6.12            htmltools_0.3.6         
##  [10] viridis_0.4.0            magrittr_1.5             memoise_1.1.0           
##  [13] limma_3.34.3             Biostrings_2.46.0        prettyunits_1.0.2       
##  [16] colorspace_1.3-2         blob_1.1.0               rrcov_1.4-3             
##  [19] dplyr_0.7.4              RCurl_1.95-4.8           tximport_1.6.0          
##  [22] lme4_1.1-14              bindr_0.1                zoo_1.8-0               
##  [25] glue_1.2.0               gtable_0.2.0             zlibbioc_1.24.0         
##  [28] XVector_0.18.0           MatrixModels_0.4-1       car_2.1-6               
##  [31] kernlab_0.9-25           prabclus_2.2-6           DEoptimR_1.0-8          
##  [34] SparseM_1.77             VIM_4.7.0                scales_0.5.0            
##  [37] mvtnorm_1.0-6            DBI_0.7                  GGally_1.3.2            
##  [40] edgeR_3.20.1             Rcpp_0.12.14.1           sROC_0.1-2              
##  [43] viridisLite_0.2.0        xtable_1.8-2             progress_1.1.2          
##  [46] laeken_0.4.6             bit_1.1-12               mclust_5.4              
##  [49] DT_0.2                   vcd_1.4-4                htmlwidgets_0.9         
##  [52] FNN_1.1                  RColorBrewer_1.1-2       fpc_2.1-10              
##  [55] modeltools_0.2-21        pkgconfig_2.0.1          reshape_0.8.7           
##  [58] XML_3.98-1.9             flexmix_2.3-14           nnet_7.3-12             
##  [61] locfit_1.5-9.1           labeling_0.3             rlang_0.1.4             
##  [64] reshape2_1.4.2           munsell_0.4.3            tools_3.4.2             
##  [67] RSQLite_2.0              pls_2.6-0                evaluate_0.10.1         
##  [70] stringr_1.2.0            cvTools_0.3.2            yaml_2.1.15             
##  [73] bit64_0.9-7              robustbase_0.92-8        bindrcpp_0.2            
##  [76] nlme_3.1-131             mime_0.5                 quantreg_5.34           
##  [79] biomaRt_2.34.0           compiler_3.4.2           pbkrtest_0.4-7          
##  [82] beeswarm_0.2.3           e1071_1.6-8              statmod_1.4.30          
##  [85] tibble_1.3.4             robCompositions_2.0.6    pcaPP_1.9-72            
##  [88] stringi_1.1.6            highr_0.6                lattice_0.20-35         
##  [91] trimcluster_0.1-2        Matrix_1.2-12            nloptr_1.0.4            
##  [94] lmtest_0.9-35            cowplot_0.9.1            data.table_1.10.4-3     
##  [97] bitops_1.0-6             httpuv_1.3.5             rtracklayer_1.38.1      
## [100] R6_2.2.2                 RMySQL_0.10.13           KernSmooth_2.23-15      
## [103] gridExtra_2.3            vipor_0.4.5              boot_1.3-20             
## [106] MASS_7.3-47              assertthat_0.2.0         rhdf5_2.22.0            
## [109] rprojroot_1.2            rjson_0.2.15             GenomicAlignments_1.14.1
## [112] Rsamtools_1.30.0         GenomeInfoDbData_0.99.1  diptest_0.75-7          
## [115] mgcv_1.8-22              grid_3.4.2               class_7.3-14            
## [118] minqa_1.2.4              rmarkdown_1.8            shiny_1.0.5             
## [121] ggbeeswarm_0.6.0
```

# Author contributions

A.T.L.L. developed and tested the workflow on all datasets.
A.T.L.L. and D.J.M. implemented improvements to the software packages required by the workflow.
J.C.M. provided direction to the software and workflow development.
All authors wrote and approved the final manuscript.

# Competing interests

No competing interests were disclosed.

# Grant information

A.T.L.L. and J.C.M. were supported by core funding from Cancer Research UK (award no. A17197).
D.J.M. was supported by a CJ Martin Fellowship from the National Health and Medical Research Council of Australia.
D.J.M and J.C.M. were also supported by core funding from EMBL.

# Acknowledgements

We would like to thank Antonio Scialdone for helpful discussions, as well as Michael Epstein, James R. Smith and John Wilson-Kanamori for testing the workflow on other datasets.

# References

