# nonlinear iBag

<!-- badges: start -->
<!-- badges: end -->

The goal of iBag is to Analyzing data from multi-platform genomics experiments combined with patients’ clinical outcomes helps us understand the complex biological processes that characterize a disease, as well as how these processes relate to the development of the disease. 


## Installation

You can install the released version of iBag from github with:

``` r
devtools::install_github("schi006/nonlineariBAG")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(iBag)
#fit mechanistic model 
M<-mechmodel(meth,mrna,cnv,dsurv)

NRUNS=1000
delta=0.05
#bart clinical model
c_bart<-clinicalmodel_survbart(M,nruns,delta)
c_bart$
#linear clinical model
c_linear<-clinicalmodel_linear(M,nruns,delta)

```

