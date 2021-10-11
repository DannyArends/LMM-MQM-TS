
setwd("D:/Ddrive/LMM-MQM-TS/")
source("LMMMQMts.R")

setwd("D:/Ddrive/LMM-MQM-TS/BFMIxB6N")

genotypes <- read.table("genotypes.txt", sep = "\t", na.strings=c("", "NA", "-"), check.names=FALSE)
phenotypes <- read.table("phenotypes.txt", sep = "\t")
covariates <- read.table("covariates.txt", sep = "\t")
map <- read.table("map.txt", sep = "\t")

# Marker we want to compensate for
markers <- c("UNC5048297")
verbose <- TRUE
mar.code <- c("N", "H", "B")
window.size <- 5
cov.type <- c("factor", "numeric", "factor", "factor")
control <- lmeControl(opt="optim", msMaxIter = 150) # Allow some more iterations (default is 50)

results <- LMMMQMts(genotypes, phenotypes, covariates, map, markers = markers, mar.code = mar.code, cov.type = cov.type, control = control,verbose = TRUE)

