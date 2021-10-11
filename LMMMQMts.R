#
# LMM-MQM timeseries analysis
#

library(nlme)

# genotypes     MxN matrix (M = Markers, N = individuals)
# phenotypes    NxT matrix (N = individuals, T = Timepoints), timepoints are dXX, where XX is the day
# covariates    NxC matrix (N = individuals, C = covariates)
# map:          Mx2 matrix (M = Markers, columns: chr, pos (Mb)
# markers:      Character vector with marker names (M1, M2, M3, .. Mx) to compensate for in multiple QTL mapping
# window.size:  Window.size defines when a marker is dropped from the MQM model
# mar.code:     Character vector of length 3 specifying how markers in the genotypes are encoded (e.g. A, H, B, or AA, AB, BB)
# cov.type:     Character vector specifying how each covariate should be interpreted (allowed: factor, numeric)
LMMMQMts <- function(genotypes, phenotypes, covariates, map, markers = NULL, window.size = 5, 
                     mar.code = c("A", "H", "B"), cov.type = rep("factor", ncol(covariates)),
                     control = lmeControl(opt="optim", msMaxIter = 150), verbose = TRUE)
{
  cat("Phase 1a - Checking input data\n")
  if(dim(genotypes)[2] != dim(phenotypes)[1]) stop("Number of individuals doesn't match between genotypes and phenotypes")
  if(dim(genotypes)[1] != dim(map)[1]) stop("Number of markers doesn't match between genotypes and map")
  if(!all(rownames(phenotypes) == colnames(genotypes))) stop("Individual names do NOT match between genotypes and phenotypes")
  if(!all(rownames(genotypes) == rownames(map))) stop("Marker names do NOT match between genotypes and map")
  if(dim(phenotypes)[1] != dim(covariates)[1]) stop("Number of individuals doesn't match between phenotypes and covariates")
  if(dim(map)[2] != 2) stop("Map needs to have only 2 columns: \"chr\" and \"pos\"")
  if(length(mar.code) != 3) stop("The mar.code parameter should be of length 3")
  if(length(cov.type) != dim(covariates)[2]) stop("The cov.type parameter should specify a type for each cofactor (allowed: factor or numeric)")

  gtscodes <- unique(unlist(genotypes))
  if(length(na.omit(gtscodes)) > 3) stop("genotypes are encoded by more than 4 groups. Please use: AA, AB, BB or A,H,B or similar coding")
  if(! all(na.omit(gtscodes) %in% mar.code)) stop("genotype encoding in genotypes does not match parameter: mar.code")
  
  # Check the column names of the map
  # Check that the positions on the map are smaller than 1 mio (prevent users from giving basepair positions)
  
  n.mar <- dim(genotypes)[1]
  n.ind <- dim(genotypes)[2]
  n.cov <- dim(covariates)[2]
  n.cov.num <- sum(cov.type == "numeric")
  n.cov.fac <- sum(cov.type == "factor")
  n.tp <- dim(phenotypes)[2]

  cat("- Markers: ", n.mar, ", Individuals: ", n.ind, ", Covariates (f/n): ", n.cov, " (",n.cov.fac, "/",n.cov.num,")", ", Timepoints: ", n.tp, "\n", sep="")

  ind.names <- rownames(phenotypes)
  cov.names <- colnames(covariates)
  cov.num.names <- colnames(covariates)[which(cov.type == "numeric")]
  cov.fac.names <- colnames(covariates)[which(cov.type == "factor")]
  mar.names <- rownames(map)
  timepoints <- as.numeric(gsub("d", "", colnames(phenotypes)))

  cat("Phase 1b - Restructure into long format\n")
  lnames <- paste0(rep(ind.names, n.tp), ":", unlist(lapply(timepoints,rep,n.ind)))
  numdata <- matrix(NA, n.ind * n.tp, 2 + length(cov.num.names), dimnames=list(lnames, c("response", "time", cov.num.names)))
  covdata <- matrix(NA, n.ind * n.tp, 1 + length(cov.fac.names), dimnames=list(lnames, c("individual", cov.fac.names)))
  for(tp in timepoints) {
    numdata[paste0(ind.names, ":", tp), "response"] <- phenotypes[, paste0("d",tp)]
    numdata[paste0(ind.names, ":", tp), "time"] <- rep(tp, n.ind)
    covdata[paste0(ind.names, ":", tp), "individual"] <- ind.names
    for(cn in cov.num.names) { numdata[paste0(ind.names, ":", tp), cn] <- covariates[, cn]; }
    for(cn in cov.fac.names) { covdata[paste0(ind.names, ":", tp), cn] <- covariates[, cn]; }
  }

  cat("Numeric data: ", dim(numdata)[1], "x", dim(numdata)[2], ", Factor data: ", dim(covdata)[1], "x", dim(covdata)[2], "\n", sep="")
  
  data.com <- data.frame(numdata, covdata)
  ismissing <- which(apply(apply(data.com,1,is.na),2,any))
  if(length(ismissing) > 0){
    data.com <-  data.com[-ismissing,]
    warning("WARN: Removing ",length(ismissing)," observations from model due to missing data\n", sep = "")
  }

  cat("Phase 2a - Model selection (Fixed effect covariates)\n")
  covToTest <- c(cov.names, "time")
  randomE <- "~ 1 | individual"
  currentNull <- "1"

  repeat{
    null <- paste0("response ~ ", currentNull)
    AICs <- c()
    for(cn in covToTest){
      alt <- paste0("response ~ ",currentNull, " + ", cn)
      m0 <- lme(fixed = as.formula(null), data = data.com, random = as.formula(randomE), control=control, method = "ML")
      m1 <- lme(fixed = as.formula(alt), data = data.com, random = as.formula(randomE), control=control, method = "ML")
      if(verbose) cat(null, "versus", alt, ", AIC = ", diff(AIC(m0,m1)[,2]), "\n")
      AICs <- c(AICs, diff(AIC(m0,m1)[,2]))
    }
    if(all(AICs > -10)) break;
    cat("minimum AIC for: ", covToTest[which.min(AICs)], " at ", AICs[which.min(AICs)], "\n",sep = "")
    currentNull <- paste0(currentNull, " + ", covToTest[which.min(AICs)])
    covToTest <- covToTest[-which.min(AICs)]
  }

  cat("Phase 2b - Model selection (Random time effect?)\n")
  currentNull <- paste0("response ~ ", currentNull)
  randomET <- "~ time | individual"
  m0 <- lme(fixed = as.formula(currentNull), data = data.com, random = as.formula(randomE), control=control, method = "ML")
  m1 <- lme(fixed = as.formula(currentNull), data = data.com, random = as.formula(randomET), control=control, method = "ML")
  if(verbose) cat(randomE,"versus", randomET, diff(AIC(m0,m1)[,2]), "\n")
  if(diff(AIC(m0,m1)[,2]) < -10) randomE <- randomET

  cat("Phase 2c - Model selection (Add time inflection points)\n")
  pow <- 2
  repeat{
    alt <- paste0(currentNull, " + I(time^",pow, ")")
    m0 <- lme(fixed = as.formula(currentNull), data = data.com, random = as.formula(randomE), control=control, method = "ML")
    m1 <- lme(fixed = as.formula(alt), data = data.com, random = as.formula(randomE), control=control, method = "ML")
    if(verbose) cat(currentNull,"versus", alt, diff(AIC(m0,m1)[,2]), "\n")
    if(diff(AIC(m0,m1)[,2]) > -10) break;
    pow <- pow + 1
    currentNull <- alt
  }
  cat("Best model found: ", currentNull, "\n", sep = "")

  cat("Phase 2d - Checking markers for MQM mapping\n")
  markerData <- matrix(NA, length(lnames), length(markers), dimnames = list(lnames, markers))
  for(m in markers){
    mEff <- factor(genotypes[m, ], levels = mar.code)
    mEff <- rep(mEff, length(timepoints))
    mData <- matrix(mEff, length(mEff),1, dimnames = list(lnames, m))
    markerData[,m] <- as.character(mEff)
    data.com <- data.frame(numdata, covdata, mData)

    ismissing <- which(apply(apply(data.com,1,is.na),2,any))
    if(length(ismissing) > 0){
      data.com <-  data.com[-ismissing,]
      cat("WARN: Marker: ",m,", removing ", length(ismissing), " observations from model (missing data)\n", sep = "")
    }

    alt <- paste0(currentNull, " + ", m, " + ", m, ":time")
    m0 <- lme(fixed = as.formula(currentNull), data = data.com, random = as.formula(randomE), control=control, method = "ML")
    m1 <- lme(fixed = as.formula(alt), data = data.com, random = as.formula(randomE), control=control, method = "ML")
    if(verbose) cat(currentNull,"versus", alt, diff(AIC(m0,m1)[,2]), "\n")
    if(diff(AIC(m0,m1)[,2]) > -10) {
      cat("WARN: Marker: ", m, " is included. However, no evidence for inclusion (AIC=", round(diff(AIC(m0,m1)[,2]),2), ")\n", sep = "")
    }
  }

  ismissing.md <- which(apply(markerData, 1, function(x){any(is.na(x))}))
  if(length(ismissing.md) > 0){
    cat("WARN: markers selected for MQM leads to removing ", length(ismissing.md), " observations from model (missing data)\n", sep = "")
  }

  cat("Phase 3 - QTL mapping\n")
  pvalues <<- c()
  for(m in rownames(genotypes)){
    mEff <- factor(genotypes[m, ], levels = mar.code)
    mEff <- rep(mEff, length(timepoints))
    mData <- matrix(mEff, length(mEff),1, dimnames = list(lnames, m))
    data.com <- data.frame(numdata, covdata, markerData, mData)
    # Remove the MQM markers selected if they're close to the marker
    chr <- map[m, "chr"]
    p.start <- map[m, "pos"] - window.size
    p.end <- map[m, "pos"] + window.size
    isClose <- rownames(map)[which(map[,"chr"] == chr & map[, "pos"] > p.start & map[, "pos"] < p.end)]
    if(any(markers %in% isClose)){
      idx <- which(markers %in% isClose)
      markersC <- markers[-idx]
      if(length(markersC) > 0){
        cat("Mapping marker ",m," using ", length(markersC), " MQM markers\n")
        markerNull <<- paste0(currentNull, " + ", paste0(markersC, " + ", paste0(markersC, ":time"), collapse=" + "))
      }else{
        cat("Mapping marker ",m," using 0 MQM marker\ns")
        markerNull <<- currentNull
      }
    }else{
      if(verbose) cat("Mapping marker ",m," using ", length(markers), " MQM markers\n")
      markerNull <<- paste0(currentNull, " + ", paste0(markers, " + ", paste0(markers, ":time"), collapse=" + "))
    }
    ismissing <- which(apply(apply(data.com, 1, is.na),2,any))
    if(length(ismissing) > 0){
      data.com <- data.com[-ismissing,]
    }
    if(length(ismissing) > length(ismissing.md)){
      cat("WARN: Marker: ",m,", removing ", length(ismissing) - length(ismissing.md), " observations from model (missing data)\n", sep = "")
    }
    if(verbose) cat("Mapping marker: ", m, " using ", nrow(data.com), " observations on ",length(unique(data.com[, "individual"]))," individuals\n",sep="")
    tryCatch(
      {
        mfstr <<- paste0(markerNull, " + ", m, " + ", m, ":time")
        mmstr <<- paste0(markerNull, " + ", m)

        mFull <- lme(fixed = as.formula(mfstr), data = data.com, random = as.formula(randomE), control=control, method = "ML")
        mMain <- lme(fixed = as.formula(mmstr), data = data.com, random = as.formula(randomE), control=control, method = "ML")
        mNull <- lme(fixed = as.formula(markerNull), data = data.com, random = as.formula(randomE), control=control, method = "ML")

        pTime <- anova(mMain, mFull)["mFull", "p-value"]
        pMain <- anova(mNull, mMain)["mMain", "p-value"]
        pBoth <- anova(mNull, mFull)["mFull", "p-value"]
        pvalues <<- rbind(pvalues, c(pBoth, pMain, pTime))
      }
    , error = function(e){ 
        pvalues <<- rbind(pvalues, c(NA, NA, NA))
        warning(paste0("Error occured at marker: ", m, ":", e))
      })
  }
  rownames(pvalues) <- rownames(genotypes)
  colnames(pvalues) <- c("All", "Main", "Time")
  return(-log10(pvalues))
}

# results:      The object returned by the LMMMQMts function
# map:          Mx2 matrix (M = Markers, columns: chr, pos (Mb)
# what:         What should we plot: All = full model results, Main = Main effects, Time = maker x time interaction
# gap:          Gap size between chromosomes
# pch:          Plot symbol
# cex:          Plot symbol magnification
# chr.col       Chromosome colors
plotEffects <- function(results, map, what = c("All", "Main", "Time"), gap = 10, pch = 20, cex=1, chr.col = c("coral", "blue")){
  chrs <- unique(map[,"chr"])
  map.sorted <- NULL
  chr.lengths <- c()
  chr.starts <- c(0)
  chrmids <- c()
  i <- 1
  for (chr in chrs) {
    onChr <- which(map[,"chr"] == chr)
    map.sorted <- rbind(map.sorted, map[onChr,])
    chr.lengths <- c(chr.lengths, max(map[onChr, "pos"]))
    chr.starts <- c(chr.starts, chr.starts[i] + max(map[onChr, "pos"]) + gap)
    i <- i + 1
  }

  chr.start <- chr.starts[-1]
  chr.ends <- chr.start + chr.lengths
  names(chr.starts) <- chrs
  names(chr.lengths) <- chrs

  for (x in chrs) {
    chrmid <- as.numeric(chr.lengths[x]/2) + as.numeric(chr.starts[x])
    chrmids <- c(chrmids, chrmid)
  }  
  plot(x = c(-gap, tail(chr.starts,1)), y = c(0, max(results,na.rm=TRUE) * 1.2), t = 'n', xlab="Chromosome", ylab="-log10[P]",xaxt='n', xaxs="i", yaxs="i", las=2, main=paste0("Effect profile"))

  i <- 1
  for (chr in chrs) {
    onChr <- rownames(map[map[,"chr"] == chr,])
    points(x=chr.starts[chr] + map[onChr,"pos"], y = results[onChr, what[1]], t ='p', pch = pch, cex = cex, col = chr.col[(i %% 2) + 1])
    i <- i + 1
  }
  axis(1, chrs, at = chrmids)
  abline(h = -log10(0.01/nrow(results)), col="green", lty=3)
  abline(h = -log10(0.05/nrow(results)), col="orange", lty=3)
}
