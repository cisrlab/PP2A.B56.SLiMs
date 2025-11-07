---
output:
  word_document: default
  html_document: default
---
itcpredictr readme
================

- [1 itcpredictr](#1-itcpredictr)
  - [1.1 INSTALLATION](#11-installation)
    - [1.1.1 R/RStudio](#111-rrstudio)
    - [1.1.2 MEME Suite](#112-meme-suite)
  - [1.2 Usage](#12-usage)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 itcpredictr

## 1.1 INSTALLATION

### 1.1.1 R/RStudio

First we will install and set up R/RStudio

– Install [R](https://cran.r-project.org)

– Optional, Install
[RStudio](https://posit.co/download/rstudio-desktop/)

– Within R/RStudio:

– Install [bioconductor](https://www.bioconductor.org/)

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")                  
BiocManager::install(version = "3.19") 
```

– Install Biostrings

``` r
BiocManager::install(“Biostrings”) 
```

– Install GenomicAlignments

``` r
BiocManager::install(“GenomicAlignments”) 
```

– Install devtools

``` r
install.packages(“devtools”) 
```

– install itcpredictr by first expanding the tar.gz file. Then use
devtools to install

``` r
devtools::install(“itcpredictr-main/itcpredictr-main/”) 
```

### 1.1.2 MEME Suite

The *itcpredictr* package makes calls to the [MEME
Suite](https://meme-suite.org/meme/) on the command-line to provide the
most accurate predictions. You can look at the
[installation](https://meme-suite.org/meme/doc/install.html?man_type=web)
instructions here or try the following depending upon the operating
system.

#### 1.1.2.1 MACOS

Use macports or install from source. Make sure that the meme and fimo
commands can be accessed from the command line by adjusting your *PATH*
variable

Using MacPorts: port install meme

From Source: review <https://meme-suite.org/meme/doc/install.html>

#### 1.1.2.2 LINUX

— Download meme source:

    wget <https://meme-suite.org/meme/meme-software/5.5.6/meme-5.5.6.tar.gz>

— Expand:

    tar zxvf meme-5.5.6.tar.gz

— Read requirements in MEME Suite’s INSTALL document

— Build meme and install

    cd meme-5.5.6

    ./configure --prefix=\$HOME --with-url=<http://meme-suite.org> --enable-build-libxml2 --enable-build-libxslt

    make –j4

    make install

#### 1.1.2.3 WINDOWS

Currently itcpredictr can use the meme suite built using
[Cygwin](https://www.cygwin.com/)

Install [Cygwin](https://www.cygwin.com/) with packages: *gcc*
*automake* *wget* *python-3.9* *perl* *make* *zlib* *zlib-devel*
*imagemagick* *ghostscript*

Within your Cygwin instance:

— Download meme source

    wget <https://meme-suite.org/meme/meme-software/5.5.6/meme-5.5.6.tar.gz>

— Build meme and install

    cd meme-5.5.6

    ./configure --prefix=\$HOME --with-url=<http://meme-suite.org> --enable-build-libxml2 --enable-build-libxslt

    make –j4

    make install

## 1.2 Usage

The *itcpredictr* package contains the *getITC_2024_0_27_cv()* call that
makes the prediction of itc given an extended amino acid sequence.
Within the amino acid sequence, dots (\\) can separate the motif
sequence from the full sequence. Phosporylation sites can be indicates
with a (pS) or (sP) to indicate whether the phosporylation site was
shown to be reflective of an D or E within the amino acid sequence.

``` r
suppressPackageStartupMessages(library(itcpredictr)) 
suppressPackageStartupMessages(getITC_2024_06_27_cv("LDTLRETQE")) 
```

    ## [1] 5.967888

``` r
getITC_2024_06_27_cv("K.L(pS)PIIED(pS)") 
```

    ## [1] 2.22189

You can also pass in a vector of sequences.

``` r
seq <- c(
        "WTSFFSG.CSPIEEEAH",
        "K.L(pS)PIIED(pS)",
        "LDTLRETQE",
        "LSIKK.LSPIIEDSREA",
        "WTSFFSG.LSPIE(Sp)EAH",
        "L(pS)PIIEDDREADH"
    )
seq_df <- data.frame(
    seq = seq,
    itc = getITC_2024_06_27_cv(seq),
    stringsAsFactors = FALSE
)

knitr::kable(seq_df)
```

| seq                  |       itc |
|:---------------------|----------:|
| WTSFFSG.CSPIEEEAH    | 49.074395 |
| K.L(pS)PIIED(pS)     |  2.221890 |
| LDTLRETQE            |  5.967887 |
| LSIKK.LSPIIEDSREA    |  6.747849 |
| WTSFFSG.LSPIE(Sp)EAH | 57.664617 |
| L(pS)PIIEDDREADH     |  1.411930 |
