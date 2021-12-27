\name{CoCoCat}
\alias{CoCoCat}
\title{A concept-generative fine-mapping method from causal diagram perspective}
\usage{
CGWAS()
}
\description{
A Concept-generative multi-stage GWAS fine mapping method with each stage performing Conditional regression by accounting for all the possible Causal diagram structure (CoCoCat).
}
\examples{
#### import data and run the CoCoCat function
source('R/Fun_cgwas.R')

### read the example dataset
#dt_fm <- read.csv('sampledata/SampleData.csv')
causalSNP <- CGWAS(data0=dt_fm, Pvalue=0.05, family0='gaussian')
# causalSNP: "rs2450390"
}
