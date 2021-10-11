
setwd("D:/Ddrive/LMM-MQM-TS/")
source("LMMMQMts.R")

setwd("D:/Ddrive/LMM-MQM-TS/S1xB6N")

genotypes <- read.table("genotypes.txt", sep = "\t", check.names=FALSE, colClasses="character")
phenotypes <- read.table("phenotypes.txt", sep = "\t")
covariates <- read.table("covariates.txt", sep = "\t")
map <- read.table("map.txt", sep = "\t")

covariates <- covariates[,-c(1,3)]

results <- LMMMQMts(genotypes, phenotypes, covariates, map, markers = c("gUNC5046545", "gUNC10595065"), mar.code = c("-1", "0", "1"))

