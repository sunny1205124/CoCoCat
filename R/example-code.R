#### import data and run the CoCoCat function
source('R/Fun_cgwas.R')

### read the example dataset
dt_fm <- read.csv('sampledata/SampleData.csv')
causalSNP <- CGWAS(data0=dt_fm, Pvalue=0.05, family0='gaussian')
# causalSNP: "rs2450390"
