Gottgens lab single cell RNAseq toolkit.
================
Wajid Jawaid
2017-02-05

<!-- README.md is generated from README.Rmd. Please edit that file -->
Installation
============

**bglab** can be installed from Github.

Github
------

``` r
library(devtools)
install_github("wjawaid/bglab")
```

Sample code
===========

This example code uses data from our gastrulation paper. Load the library, download counts and metdata from <http://gastrulation.stemcells.cam.ac.uk>.

``` r
library(bglab)

system("wget http://gastrulation.stemcells.cam.ac.uk/data/counts.tar.gz")
untar("counts.tar.gz")
counts <- read.table("counts.txt", header = TRUE, row.names = 1, sep = " ", stringsAsFactors = FALSE)
counts <- as.matrix(counts)
meta <- read.table("http://gastrulation.stemcells.cam.ac.uk/data/metadataAll.txt", header = TRUE, stringsAsFactors = FALSE)
```

Download the appropriate gene and ensembl IDs from ensembl.

``` r
library(biomaRt)
host <- "oct2014.archive.ensembl.org"   # Ensembl 77
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host=host)
ensembl <- useDataset("mmusculus_gene_ensembl", ensembl)
geneTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = ensembl)
colnames(geneTable) <- c("Ensembl_Gene_ID", "Associated_Gene_Name")
geneTable <- geneTable[!duplicated(geneTable$Ensembl_Gene_ID),]
```

Extract QC produced by HTSeq and the ERCCs from the previously downloaded counts table.

``` r
qcInd <- grep("^__", rownames(counts))
erccInd <- grep("^ERCC", rownames(counts))

ercc <- counts[erccInd,]
qc <- counts[qcInd,]
counts <- counts[-c(qcInd, erccInd),]
```

Now start loading the data into the single cell dataset (SCD) object which is required for the *bglab* package.

``` r
scd <- newSCD("RNAseq", counts = counts, genoData = geneTable, spike = ercc, qc = qc, phenoData = meta)
```

Perform QC, this will output a separate pdf for each lane of the flowCell (If that is encoded in the metadata) otherwise the user will need to supply a different column in the metadata. For convenience this function will also perform normalisation - as per the DESeq size factor. In future versions this may be separated.

``` r
scd <- performQC(scd, pdf = "qc")
```

Select technically variable genes as per Brennecke et al.

``` r
scd <- techVar(scd, useERCC = FALSE, meanForFit = 10)
plotTechVar(scd@technicalNoise)
```

Visualise the datasets using dimensionality reduction of choice.

``` r
scd <- runTSNE(scd, ndims = 2, seed = 0, verbose = TRUE, use_dist = TRUE)
plot(scd, reduceMethod = "tsne", colorBy = "cluster", plotCols = sort(unique(pData(scd)[,"cluster"])), legPos = "bottomleft", outline = "black")

scd <- runPCA(scd, scale. = TRUE)
plot(scd, reduceMethod = "pca", colorBy = "cluster", plotCols = sort(unique(pData(scd)[,"cluster"])), legPos = "topright", outline = "black")

dmap <- diffuseMat2(exprs(scd))
plot(dmap$vectors, bg = pData(scd)[,"cluster"], col = "black", pch = 21)
plot(dmap$vectors, bg = ggCol(pData(scd)[,"embryoStage"]), col = "black", pch = 21)
```

This si a large package and many other functions are available. Look at the help documents if things are not clear, feel free to contact me by email <wj241@cam.ac.uk>
