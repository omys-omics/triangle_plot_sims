setwd("C:/Users/bwien/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/Triangle_Plots_best_practices")

library(readr)
library(vcfR)
library(adegenet)
library(StAMPP)
library(ggplot2)
library(hierfstat)
library(dartR)
library(triangulaR)
library(ggpubr)
library(reshape)
library(ggbreak)
library(bgchm)

######################################
#####      Define Functions     ######
######################################

# MODIFY min_mac from SNPfiltR to accomodate pipe character in vcf files
min.mac <- function (vcfR, min.mac = NULL) 
{
  if (!inherits(vcfR, what = "vcfR")) {
    stop("specified vcfR object must be of class 'vcfR'")
  }
  if (max(nchar(gsub(",", "", vcfR@fix[, "ALT"])), na.rm = T) > 
      1) {
    stop("Input vcf contains SNPs with > 2 alleles. MAC is calculated under a strict assumption that a single SNP can only possess two alleles. Please use 'filter_biallelic(vcfR)' to remove multi-allelic sites before implementing a MAC filter.")
  }
  if (is.null(min.mac)) {
    message("no filtering cutoff provided, vcf will be returned unfiltered")
    gt.matrix <- vcfR::extract.gt(vcfR)
    missingness.og <- sum(is.na(gt.matrix))
    gt.matrix[gt.matrix == "0|0"] <- 0
    gt.matrix[gt.matrix == "0|1"] <- 1
    gt.matrix[gt.matrix == "1|0"] <- 1
    gt.matrix[gt.matrix == "1|1"] <- 2
    class(gt.matrix) <- "numeric"
    missingness.new <- sum(is.na(gt.matrix))
    if (missingness.og != missingness.new) {
      stop("Unrecognized genotype values in input vcf. Only allowed non-missing genotype inputs are '0/0','0/1','1/0','1/1'.")
    }
    sfs <- rowSums(gt.matrix, na.rm = TRUE)
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i, ]))) {
      }
      else {
        sfs[i] <- (sum(!is.na(gt.matrix[i, ])) * 2 - 
                     sfs[i])
      }
    }
    graphics::hist(sfs, main = "folded SFS", xlab = "MAC")
    return(vcfR)
  }
  else {
    if (!inherits(min.mac, what = "numeric")) {
      stop("specified min.mac must be numeric")
    }
    gt.matrix <- vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0|0"] <- 0
    gt.matrix[gt.matrix == "0|1"] <- 1
    gt.matrix[gt.matrix == "1|1"] <- 2
    class(gt.matrix) <- "numeric"
    sfs <- rowSums(gt.matrix, na.rm = TRUE)
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i, ]))) {
      }
      else {
        sfs[i] <- (sum(!is.na(gt.matrix[i, ])) * 2 - 
                     sfs[i])
      }
    }
    graphics::hist(sfs, main = "folded SFS", xlab = "MAC")
    graphics::abline(v = min.mac - 1, col = "red")
    p <- round((sum(sfs < min.mac)/length(sfs)) * 100, 2)
    message(p, "% of SNPs fell below a minor allele count of ", 
            min.mac, " and were removed from the VCF")
    vcfR <- vcfR[sfs >= min.mac, ]
    return(vcfR)
  }
}




cross <- function(m = NULL, parents1 = NULL, parents2 = NULL) {
  n <- data.frame(matrix(ncol = 0, nrow = nrow(m)))
  for (i in 1:length(parents1)) {
    p1 <- parents1[i]
    p2 <- parents2[i]
    off <- c()

    for (g in 1:nrow(m)) {
      
      if (m[g,p1] == 0 & m[g,p2] == 0) {
        off <- c(off, 0)
        next
      }
      
      if (m[g,p1] == 0 & m[g,p2] == 1) {
        r <- runif(1, min = 0, max = 1)
        if (r < 0.5) {
          off <- c(off, 0)
          next
        }
        if (r > 0.5) {
          off <- c(off, 1)
          next
        }
      }
      
      if (m[g,p1] == 0 & m[g,p2] == 2) {
        off <- c(off, 1)
        next
      }
      
      if (m[g,p1] == 1 & m[g,p2] == 0) {
        r <- runif(1, min = 0, max = 1)
        if (r < 0.5) {
          off <- c(off, 0)
          next
        }
        if (r > 0.5) {
          off <- c(off, 1)
          next
        }
        next
      }
      
      if (m[g,p1] == 1 & m[g,p2] == 1) {
        r <- runif(1, min = 0, max = 1)
        if (r < 0.25) {
          off <- c(off, 0)
          next
        }
        if (r > 0.25 & r < 0.75) {
          off <- c(off, 1)
          next
        }
        if (r > 0.75) {
          off <- c(off, 2)
          next
        }
        next
      }
      
      if (m[g,p1] == 1 & m[g,p2] == 2) {
        r <- runif(1, min = 0, max = 1)
        if (r < 0.5) {
          off <- c(off, 1)
          next
        }
        if (r > 0.5) {
          off <- c(off, 2)
          next
        }
        next
      }
      
      if (m[g,p1] == 2 & m[g,p2] == 0) {
        off <- c(off, 1)
        next
      }
      
      if (m[g,p1] == 2 & m[g,p2] == 1) {
        r <- runif(1, min = 0, max = 1)
        if (r < 0.5) {
          off <- c(off, 1)
          next
        }
        if (r > 0.5) {
          off <- c(off, 2)
          next
        }
        next
      }
      
      if (m[g,p1] == 2 & m[g,p2] == 2) {
        off <- c(off, 2)
        next
      }
    }
    n <- cbind(n, off)
  }
  return(n)  
}





makeHybrids <- function(vcf = NULL, seed = 36428709, samples.p1 = NULL, samples.p2 = NULL) {
  vcfR <- read.vcfR(vcf)
  vcfR <- min.mac(vcfR, min.mac = 1)
  
  m <- extract.gt(vcfR)
  # recode missing data
  m[is.na(m)] <- -9
  # recode to allele counts
  m[m=="0|0"] <- 0
  m[m=="0|1"] <- 1
  m[m=="1|0"] <- 1
  m[m=="1|1"] <- 2
  m[m=="0/0"] <- 0
  m[m=="0/1"] <- 1
  m[m=="1/0"] <- 1
  m[m=="1/1"] <- 2
  
  # set seed
  set.seed(seed)
  
  # make f1s
  p0.parents <- sample(colnames(m)[1:samples.p1],20)
  p20.parents <- sample(colnames(m)[(samples.p1+1):(samples.p1+samples.p2)],20)
  f1 <- cross(m = m, parents1 = p0.parents, parents2 = p20.parents)
  colnames(f1) <- paste0(rep("f1", 20), ".", 1:20)
  m <- cbind(m, f1)
  
  # make f2s
  f1.parents.1 <- sample(colnames(f1))
  f1.parents.2 <- rev(f1.parents.1)
  f2 <- cross(m = m, parents1 = f1.parents.1, parents2 = f1.parents.2)
  colnames(f2) <- paste0(rep("f2", 20), ".", 1:20)
  m <- cbind(m, f2)
  
  # make p1 backcrosses
  p0.parents <- sample(colnames(m)[1:samples.p1],20)
  f1.parents.1 <- sample(colnames(f1))
  bc1 <- cross(m = m, parents1 = p0.parents, parents2 = f1.parents.1)
  colnames(bc1) <- paste0(rep("bc1", 20), ".", 1:20)
  m <- cbind(m, bc1)
  
  # make p2 backcrosses
  p20.parents <- sample(colnames(m)[(samples.p1+1):(samples.p1+samples.p2)],20)
  f1.parents.1 <- sample(colnames(f1))
  bc2 <- cross(m = m, parents1 = p20.parents, parents2 = f1.parents.1)
  colnames(bc2) <- paste0(rep("bc2", 20), ".", 1:20)
  m <- cbind(m, bc2)
  
  return(m)
}




AIMnames <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- matrix[, pm[pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  loci <- rownames(af[af$diff >= difference,])
  
  return(loci)
}



alleleFreqDiff.matrix <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL, differences = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }

  # Filter and subset the genotypes for the two populations
  p1.gts <- matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- matrix[, pm[pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  # creat list to store matrices in
  matrix.diffs <- list()
  for (d in differences) {
    # get names of loci with allele frequency difference above threshold
    loci <- rownames(af[af$diff >= d,])
    
    # print statement
    s <- length(loci)
    print(paste0(s, " sites passed allele frequency difference threshold: ", d))
    
    # subset matrix by loci with allele frequencies difference in parental pops above threshold
    matrix.diff <- matrix[loci,]
    matrix.diffs[as.character(d)] <- list(matrix.diff)
  }

  return(matrix.diffs)
}



specFreqDiff <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- matrix[, pm[pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  # Filter by the difference threshold
  af <- af[af$diff >= difference, ]
  
  # Print statement
  s <- nrow(af)
  cat(s, "sites passed allele frequency difference threshold\n")
  
  return(af)
}


aimFreqDist.matrix <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL, difference = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- matrix[, pm[pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  # only keep loci with allele frequency difference above threshold
  af <- af[af$diff >= difference,]
  
  # polarize allele freqs
  af[af$p1>0.5,"p1"] <- 1 - af[af$p1>0.5,"p1"]
  af[af$p2<0.5,"p2"] <- 1 - af[af$p2<0.5,"p2"]
  
  # print statement
  s <- nrow(af)
  print(paste0(s, " sites passed allele frequency difference threshold"))
  return(af)
}




trueFreqDist.matrix <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL, 
                                difference = NULL, num.pars = NULL, total.samples.per.par = NULL,
                                seed = 103498) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }
  
  if(seed != "rep") {
    set.seed(seed)
  }

  
  ids <- c(sample(pm[pm$pop=="p0","id"], num.pars, replace = F), sample(pm[pm$pop=="p20","id"], num.pars, replace = F))
  sub.pm <- pm[pm$id %in% ids,]
  sub.matrix <- matrix[,sub.pm$id]
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- sub.matrix[, sub.pm[sub.pm$pop == p1, "id"]]
  p2.gts <- sub.matrix[, sub.pm[sub.pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  # get names of loci with allele frequency difference above threshold
  loci <- rownames(af[af$diff >= difference,])
  
  # make new genotype matrix with all parental individuals, only at loci that met threshold for subsampled parentals
  passed.loci.matrix <- matrix[loci,]
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- passed.loci.matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- passed.loci.matrix[, pm[pm$pop == p2, "id"]]
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Calculate allele frequency differences
  af_diff <- abs(af_p1 - af_p2)
  
  # Create a data frame with allele frequencies and differences
  af <- data.frame(p1 = af_p1, p2 = af_p2, diff = af_diff)
  
  # print statement
  s <- nrow(af)
  print(paste0(s, " sites passed allele frequency difference threshold"))
  return(af)
}




hybridIndex.matrix <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  
  if (!(identical(sort(unique(colnames(matrix))), sort(unique(pm$id))))) {
    stop("There is at least one individual in the matrix that is not in the popmap, or vice versa")
  }
  
  # get number of differences above threshold
  d <- nrow(matrix)
  print(paste0("calculating hybrid indices and heterozygosities based on ", d, " sites"))
  
  
  # Filter and subset the genotypes for the two populations
  p1.gts <- matrix[, pm[pm$pop == p1, "id"]]
  p2.gts <- matrix[, pm[pm$pop == p2, "id"]]
  
  # convert to numeric
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # Replace -9 with NA
  p1.gts[p1.gts == -9] <- NA
  p2.gts[p2.gts == -9] <- NA
  
  # Calculate allele frequencies for p1 and p2
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 2, na.rm = TRUE))) / (2 * rowSums(!is.na(p2.gts)))
  
  # Determine p1 and p2 allele based on allele frequencies
  p1.allele <- ifelse(af_p1 > af_p2, 2, 0)
  p2.allele <- ifelse(af_p2 > af_p1, 2, 0)
  
  # Create a matrix to store hybrid index scores
  n <- matrix(nrow = nrow(matrix), ncol = ncol(matrix))
  
  # Compare genotypes and assign scores
  n[matrix == p1.allele] <- 0
  n[matrix == 1] <- 1
  n[matrix == p2.allele] <- 2
  n[is.na(matrix)] <- NA
  n[matrix == -9] <- NA
  
  colnames(n) <- colnames(matrix)
  rownames(n) <- rownames(matrix)
  
  # Count alleles and non-missing genotypes for each individual
  counts <- colSums(n, na.rm = TRUE)
  #sites <- colSums(!is.na(n))
  sites <- colSums(!apply(n, MARGIN = 2, is.na))
  
  # Calculate hybrid index
  hi <- counts / (sites * 2)
  
  # Create a dataframe for the results
  tri <- data.frame(
    id = names(hi),
    pop = pm[match(names(hi), pm$id), "pop"],
    hybrid.index = hi,
    heterozygosity = colSums(matrix == 1, na.rm = TRUE) / colSums(!is.na(matrix)),
    perc.missing = colSums(is.na(matrix)) / nrow(matrix)
  )
  
  return(tri)
}




makePopmap <- function(file = NULL, vcf = NULL) {
  vcfR <- read.vcfR(vcf)
  num <- read_table(file)
  pm <- data.frame(matrix(nrow=(sum(num)+80), ncol=2))
  colnames(pm) <- c("id", "pop")
  pm$id <- c(colnames(vcfR@gt)[-1], paste0(rep("f1.",20),1:20), paste0(rep("f2.",20),1:20), paste0(rep("bc1.",20),1:20), paste0(rep("bc2.",20),1:20) )
  pm$pop <- c(rep("p0", num[,1]), rep("p20", num[,2]), rep("f1", 20), rep("f2", 20), rep("bc1", 20), rep("bc2", 20))
  return(pm)
}




####################
#   Make hybrids   #
####################
D750.pm <- makePopmap(file = "D750/D750.num.parentals.end.of.div.txt", vcf = "D750/D750.end.of.div.vcf")
D750.manual <- makeHybrids(vcf = "D750/D750.end.of.div.vcf", samples.p1 = sum(D750.pm$pop=="p0"), samples.p2 = sum(D750.pm$pop=="p20"))
write.table(D750.pm[D750.pm$pop%in%c("p0","p20"),], file = "D750/D750.end.of.div.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)

D1000.pm <- makePopmap(file = "D1000/D1000.num.parentals.end.of.div.txt", vcf = "D1000/D1000.end.of.div.vcf")
D1000.manual <- makeHybrids(vcf = "D1000/D1000.end.of.div.vcf", samples.p1 = sum(D1000.pm$pop=="p0"), samples.p2 = sum(D1000.pm$pop=="p20"))
write.table(D1000.pm[D1000.pm$pop%in%c("p0","p20"),], file = "D1000/D1000.end.of.div.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)

D2000.pm <- makePopmap(file = "D2000/D2000.num.parentals.end.of.div.txt", vcf = "D2000/D2000.end.of.div.vcf")
D2000.manual <- makeHybrids(vcf = "D2000/D2000.end.of.div.vcf", samples.p1 = sum(D2000.pm$pop=="p0"), samples.p2 = sum(D2000.pm$pop=="p20"))
write.table(D2000.pm[D2000.pm$pop%in%c("p0","p20"),], file = "D2000/D2000.end.of.div.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)





################################################
#   Calculate hi and het from available data   #
################################################

# >1000 individuals in parental pops, 3 different afdts
for (m in c("D750", "D1000", "D2000")) {
  af.diffs <- alleleFreqDiff.matrix(matrix = get(paste0(m,".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", differences = c(1,0.75,0.5))
  for (a in 1:3) {
    d <- nrow(af.diffs[a][[1]])
    temp <- hybridIndex.matrix(matrix = af.diffs[a][[1]], pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20")
    tri <- triangle.plot(temp, colors = c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837"), alpha = 0.8)
    tri <- tri + annotate(geom="text", x=-0.05, y=1.05, hjust = 0, label=paste0(d, " SNPs"), size=4) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_x_continuous(breaks = c(0,0.5,1)) +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
      
    if (a == 1) {
      stack <- tri
    } else {
      stack <- stack / tri
    }
    assign(paste0(m, ".afdts"), stack)
  }
  rm(stack)
}

D750.afdts
D1000.afdts
D2000.afdts


afdts <- ggarrange(D750.afdts[[1]], D1000.afdts[[1]], D2000.afdts[[1]], 
                   D750.afdts[[2]], D1000.afdts[[2]], D2000.afdts[[2]],
                   D750.afdts[[3]], D1000.afdts[[3]], D2000.afdts[[3]],
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                  font.label = list(size = 14),
                  #label.y = -0.2,
                  label.x = 0.2,
                  ncol = 3, nrow = 3)

ggsave("Figure5fromR.pdf", plot = afdts, path = "Figures/", width = 7, height = 6)

# accuracy and precision
for (m in c("D750", "D1000", "D2000")) {
  if (m == "D750") {sample.size.ap.afdt <- data.frame(matrix(nrow = 0, ncol = 8))}
  af.diffs <- alleleFreqDiff.matrix(matrix = get(paste0(m,".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", differences = c(1,0.75,0.5))
  for (a in 1:3) {
    if (a==1) {i<-1}
    if (a==2) {i<-0.75}
    if (a==3) {i<-0.5}
    d <- nrow(af.diffs[a][[1]])
    temp <- hybridIndex.matrix(matrix = af.diffs[a][[1]], pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20")
    for (class in c("f1", "f2", "bc1", "bc2")) {
      #expected values
      if (class=="f1") {e.hi <- 0.5}
      if (class=="f1") {e.het <- 1}
      if (class=="f2") {e.hi <- 0.5}
      if (class=="f2") {e.het <- 0.5}
      if (class=="bc1") {e.hi <- 0.25}
      if (class=="bc1") {e.het <- 0.5}
      if (class=="bc2") {e.hi <- 0.75}
      if (class=="bc2") {e.het <- 0.5}
      #observed values
      o.hi <- mean(temp[temp$pop==class,]$hybrid.index)
      o.het <- mean(temp[temp$pop==class,]$heterozygosity)
      #accuracy
      a.hi <- 1 - abs(e.hi-o.hi)/e.hi
      a.het <- 1 - abs(e.het-o.het)/e.het
      #precision
      p.hi <- mean(abs(o.hi - temp[temp$pop==class,]$hybrid.index))
      p.het <- mean(abs(o.het - temp[temp$pop==class,]$heterozygosity))
      #put in dataframe
      sample.size.ap.afdt  <- rbind(sample.size.ap.afdt , c(m, i, class, "hybrid index", e.hi, o.hi, a.hi, p.hi))
      sample.size.ap.afdt  <- rbind(sample.size.ap.afdt , c(m, i, class, "heterozygosity", e.het, o.het, a.het, p.het))
    }
    
  }
}
colnames(sample.size.ap.afdt) <- c("divergence", "afdt", "class", "metric", "expected", "observed", "accurary", "precision")
sample.size.ap.afdt$divergence <- c(rep("Low", 24), rep("Med", 24), rep("High", 24))
sample.size.ap.afdt$accurary <- as.numeric(sample.size.ap.afdt$accurary)
sample.size.ap.afdt$precision <- as.numeric(sample.size.ap.afdt$precision)
sample.size.ap.afdt[sample.size.ap.afdt$precision==0,"precision"]<-0.001

hi.a.plot.afdt <- ggplot(data = sample.size.ap.afdt[sample.size.ap.afdt$metric=="hybrid index",],  aes(fill=factor(afdt, levels = c("1", "0.75", "0.5")), y=accurary, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "afdt") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.a.plot.afdt <- ggplot(data = sample.size.ap.afdt[sample.size.ap.afdt$metric=="heterozygosity",],  aes(fill=factor(afdt, levels = c("1", "0.75", "0.5")), y=accurary, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6")) +
  ggtitle("Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "afdt") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

hi.p.plot.afdt <- ggplot(data = sample.size.ap.afdt[sample.size.ap.afdt$metric=="hybrid index",],  aes(fill=factor(afdt, levels = c("1", "0.75", "0.5")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.14), breaks = c(0, 0.07, 0.14)) +
  #theme(legend.position = "none")
  labs(fill = "afdt") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.p.plot.afdt <- ggplot(data = sample.size.ap.afdt[sample.size.ap.afdt$metric=="heterozygosity",],  aes(fill=factor(afdt, levels = c("1", "0.75", "0.5")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d")) +
  ggtitle("Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.14), breaks = c(0, 0.07, 0.14)) +
  #theme(legend.position = "none")
  labs(fill = "afdt") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

sample.size.ap.plot.afdt <- ggarrange(hi.a.plot.afdt, het.a.plot.afdt, hi.p.plot.afdt, het.p.plot.afdt,
                                 labels = c("A", "B", "C", "D"),
                                 font.label = list(size = 14),
                                 #label.y = -0.2,
                                 label.x = 0.1,
                                 ncol = 2, nrow = 2)

ggsave("FigureX2fromR.pdf", plot = sample.size.ap.plot.afdt, path = "Figures/", width = 10, height = 8)





# Fewer individuals per parental pop, afdt of 1
for (m in c("D750", "D1000", "D2000")) {
  set.seed(6412)
  for (i in c(20, 10, 5, 2)) {
    pmt <- get(paste0(m, ".pm"))
    ids <- c(sample(pmt[pmt$pop=="p0","id"], i, replace = F), sample(pmt[pmt$pop=="p20","id"], i, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
    pmt <- pmt[pmt$id %in% ids,]
    temp <- get(paste0(m, ".manual"))[,pmt$id]
    temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = 1)
    d <- nrow(temp[1][[1]])
    temp <- hybridIndex.matrix(matrix = temp[1][[1]], pm = pmt, p1 = "p0", p2 = "p20")
    tri <- triangle.plot(temp, colors = c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837"))
    tri <- tri + annotate(geom="text", x=-0.05, y=1.05, hjust = 0, label=paste0(d, " SNPs"), size=4) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_x_continuous(breaks = c(0,0.5,1)) +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
    
    if (i == 20) {
      stack <- tri
    } else {
      stack <- stack / tri
    }
    assign(paste0(m, ".inds"), stack)
  }
  rm(stack)
}

inds.1 <- ggarrange(D750.inds[[1]], D1000.inds[[1]], D2000.inds[[1]], 
                    D750.inds[[2]], D1000.inds[[2]], D2000.inds[[2]], 
                    D750.inds[[3]], D1000.inds[[3]], D2000.inds[[3]], 
                    D750.inds[[4]], D1000.inds[[4]], D2000.inds[[4]], 
                       labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                       font.label = list(size = 14),
                       #label.y = -0.2,
                       label.x = 0.1,
                       ncol = 3, nrow = 4)

ggsave("Figure4fromR.pdf", plot = inds.1, path = "Figures/", width = 7, height = 8)


# Accuracy and precision 
for (m in c("D750", "D1000", "D2000")) {
  set.seed(6412)
  if (m == "D750") {sample.size.ap <- data.frame(matrix(nrow = 0, ncol = 8))}
  for (i in c(20, 10, 5, 2)) {
    pmt <- get(paste0(m, ".pm"))
    ids <- c(sample(pmt[pmt$pop=="p0","id"], i, replace = F), sample(pmt[pmt$pop=="p20","id"], i, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
    pmt <- pmt[pmt$id %in% ids,]
    temp <- get(paste0(m, ".manual"))[,pmt$id]
    temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = 1)
    d <- nrow(temp[1][[1]])
    temp <- hybridIndex.matrix(matrix = temp[1][[1]], pm = pmt, p1 = "p0", p2 = "p20")
    for (class in c("f1", "f2", "bc1", "bc2")) {
      #expected values
      if (class=="f1") {e.hi <- 0.5}
      if (class=="f1") {e.het <- 1}
      if (class=="f2") {e.hi <- 0.5}
      if (class=="f2") {e.het <- 0.5}
      if (class=="bc1") {e.hi <- 0.25}
      if (class=="bc1") {e.het <- 0.5}
      if (class=="bc2") {e.hi <- 0.75}
      if (class=="bc2") {e.het <- 0.5}
      #observed values
      o.hi <- mean(temp[temp$pop==class,]$hybrid.index)
      o.het <- mean(temp[temp$pop==class,]$heterozygosity)
      #accuracy
      a.hi <- 1 - abs(e.hi-o.hi)/e.hi
      a.het <- 1 - abs(e.het-o.het)/e.het
      #precision
      p.hi <- mean(abs(o.hi - temp[temp$pop==class,]$hybrid.index))
      p.het <- mean(abs(o.het - temp[temp$pop==class,]$heterozygosity))
      #put in dataframe
      sample.size.ap <- rbind(sample.size.ap, c(m, i, class, "hybrid index", e.hi, o.hi, a.hi, p.hi))
      sample.size.ap <- rbind(sample.size.ap, c(m, i, class, "heterozygosity", e.het, o.het, a.het, p.het))
    }
  }
}
colnames(sample.size.ap) <- c("divergence", "sample.size", "class", "metric", "expected", "observed", "accurary", "precision")
sample.size.ap$divergence <- c(rep("Low", 32), rep("Med", 32), rep("High", 32))
sample.size.ap$accurary <- as.numeric(sample.size.ap$accurary)
sample.size.ap$precision <- as.numeric(sample.size.ap$precision)

hi.a.plot <- ggplot(data = sample.size.ap[sample.size.ap$metric=="hybrid index",],  aes(fill=factor(sample.size, levels = c("20", "10", "5", "2")), y=accurary, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6", "#6baed6")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "Sample size") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.a.plot <- ggplot(data = sample.size.ap[sample.size.ap$metric=="heterozygosity",],  aes(fill=factor(sample.size, levels = c("20", "10", "5", "2")), y=accurary, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6", "#6baed6")) +
  ggtitle("Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "Sample size") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

hi.p.plot <- ggplot(data = sample.size.ap[sample.size.ap$metric=="hybrid index",],  aes(fill=factor(sample.size, levels = c("20", "10", "5", "2")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d", "#74c476")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  #theme(legend.position = "none")
  labs(fill = "Sample size") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.p.plot <- ggplot(data = sample.size.ap[sample.size.ap$metric=="heterozygosity",],  aes(fill=factor(sample.size, levels = c("20", "10", "5", "2")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d", "#74c476")) +
  ggtitle("Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.05, 0.1)) +
  #theme(legend.position = "none")
  labs(fill = "Sample size") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

sample.size.ap.plot <- ggarrange(hi.a.plot, het.a.plot, hi.p.plot, het.p.plot,
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 14),
          #label.y = -0.2,
          label.x = 0.1,
          ncol = 2, nrow = 2)

ggsave("FigureX1fromR.pdf", plot = sample.size.ap.plot, path = "Figures/", width = 10, height = 8)


##############################################
#          Allele freq diff spectrum         #
##############################################
for (m in c("D750", "D1000", "D2000")) {
  temp <- specFreqDiff(matrix = get(paste0(m, ".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", difference = 0)
  if (m == "D750") {
    a.freq.diff.spec <- ggplot(data = temp,  aes(x=diff)) +
      geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 15) +
      #scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) +
      xlab(paste("Allele Frequency Difference")) +
      ylab(paste("Count")) +
      ylim(0,4000) +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  } else {
    a.freq.diff.spec <- a.freq.diff.spec + ggplot(data = temp,  aes(x=diff)) +
      geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 15) +
      #scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) +
      xlab(paste("Allele Frequency Difference")) +
      ylab(paste("Count")) +
      ylim(0,4000) +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
    
  }
}

a.freq.diff.spec
ggsave("Figure2fromR.pdf", plot = a.freq.diff.spec, path = "Figures/", width = 12, height = 6)


##############################################
#     Distribution of allele frequencies     #
##############################################

for (m in c("D750", "D1000", "D2000")) {
  for (d in c(0.9, 0.7, 0.5)) {
    temp <- aimFreqDist.matrix(matrix = get(paste0(m, ".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", difference = d)
    temp <- temp[,c("p1", "p2")]
    temp <- melt(temp, variable_name = "pop")
    temp[temp$value==0.5 & temp$pop=="p1","value"] <- 0.4999
    temp[temp$value==0.5 & temp$pop=="p2","value"] <- 0.5001
    
    dist <- ggplot(data = temp,  aes(x=value, fill=pop)) +
      geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 100) +
      scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) +
      xlab(paste("Frequency of Parental Allele")) +
      ylab(paste("Count")) +
      theme_classic() +
      scale_y_continuous(limits = c(0,20)) +
      theme(legend.position = "none")
    
    if (d == 0.9) {
      stack <- dist
    } else {
      stack <- stack / dist
    }
    print(stack)
  }
  assign(paste0(m, ".dist.par.allele"), stack)
  rm(stack)
}

D750.dist.par.allele
D1000.dist.par.allele
D2000.dist.par.allele


dist.par.allele <- ggarrange(D750.dist.par.allele, D1000.dist.par.allele, D2000.dist.par.allele, 
                       labels = c("A", "B", "C"),
                       font.label = list(size = 14),
                       #label.y = -0.2,
                       label.x = 0.1,
                       ncol = 3, nrow = 1)


###########################################################################
#     Distribution of allele frequencies after downsampling parentals     #
###########################################################################
for (m in c("D750", "D1000", "D2000")) {
  if (m == "D750") {false.pos.pro <- data.frame(matrix(nrow = 0, ncol = 5))}
  if (m == "D750") {breaks <- c(0,2000,4000)}
  if (m == "D1000") {breaks <- c(0,4000,8000)}
  if (m == "D2000") {breaks <- c(0,30000,60000)}
  for (n in c(20, 10, 5, 2)) {
    for (i in 1:200) {
      if (i == 1) {
        rep <- trueFreqDist.matrix(matrix = get(paste0(m, ".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", 
                                   difference = 1, num.pars = n, seed = "rep")
        next
      }
      temp <- trueFreqDist.matrix(matrix = get(paste0(m, ".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", 
                                  difference = 1, num.pars = n, seed = "rep")
      rep <- rbind(rep, temp)
    }
    true.freq <- ggplot(data = rep,  aes(x=diff)) +
      geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 30) +
      #scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) +
      xlab(paste("True Allele Frequency Difference")) +
      ylab(paste("Count")) +
      theme_classic() +
      xlim(c(-0.05,1.05)) +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
      scale_y_continuous(breaks = breaks) +
      geom_vline(aes(xintercept = quantile(diff, probs = c(0.05))), color = "#000000", size = 0.1, linetype = "dashed") +
      geom_vline(aes(xintercept = quantile(diff, probs = c(0.1))), color = "#91bfdb", size = 0.1, linetype = "solid") +
      geom_vline(aes(xintercept = quantile(diff, probs = c(0.25))), color = "#fc8d59", size = 0.1, linetype = "solid")
    
    false.pos.pro <- rbind(false.pos.pro, c(m,n,sum(rep$diff<1)/nrow(rep),sum(rep$diff<0.95)/nrow(rep), mean(rep$diff), quantile(rep$diff, probs = c(0.05)), quantile(rep$diff, probs = c(0.1)), quantile(rep$diff, probs = c(0.25))))

    
    if (n == 20) {
      stack <- true.freq
    } else {
      stack <- stack / true.freq
    }
    print(stack)
  }
  assign(paste0(m, ".true.freq"), stack)
  rm(stack)
}

colnames(false.pos.pro) <- c("sim", "sample", "false_pos_1", "false_pos_0.95", "accuracy", "0.05_quantile", "0.1_quantile", "0.25_quantile")
false.pos.pro[false.pos.pro[,"sim"]=="D750","sim"] <- "low"
false.pos.pro[false.pos.pro[,"sim"]=="D1000","sim"] <- "med"
false.pos.pro[false.pos.pro[,"sim"]=="D2000","sim"] <- "high"
round(false.pos.pro)
write.table(false.pos.pro, file = "Figures/SFig1.txt", append = F, quote = F, sep = "\t", col.names = T, row.names = F)


D750.true.freq
D1000.true.freq
D2000.true.freq


true.freq <- ggarrange(D750.true.freq[[1]], D1000.true.freq[[1]], D2000.true.freq[[1]], 
                       D750.true.freq[[2]], D1000.true.freq[[2]], D2000.true.freq[[2]], 
                       D750.true.freq[[3]], D1000.true.freq[[3]], D2000.true.freq[[3]], 
                       D750.true.freq[[4]], D1000.true.freq[[4]], D2000.true.freq[[4]], 
                       labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                       font.label = list(size = 14),
                       #label.y = -0.2,
                       label.x = 0.25,
                       ncol = 3, nrow = 4)


ggsave("Figure3_fromR.pdf", plot = true.freq, path = "Figures/", width = 9, height = 6)






# Fewer SNPs when divergence is low and d=0.5
for (m in c("D750")) {
  set.seed(12039458)
  nAIMs.ap <- data.frame(matrix(nrow = 0, ncol = 8))
  af.diffs <- alleleFreqDiff.matrix(matrix = get(paste0(m,".manual")), pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20", differences = c(0.5))
  z <- 2
  for (n in c(10, 20, 40, 60, 80, 100)) {
    af.diffs[z][[1]] <- af.diffs[1][[1]][ceiling(runif(n, 0, 315)),]
    d <- nrow(af.diffs[z][[1]])
    temp <- hybridIndex.matrix(matrix = af.diffs[z][[1]], pm = get(paste0(m, ".pm")), p1 = "p0", p2 = "p20")
    for (class in c("f1", "f2", "bc1", "bc2")) {
      #expected values
      if (class=="f1") {e.hi <- 0.5}
      if (class=="f1") {e.het <- 1}
      if (class=="f2") {e.hi <- 0.5}
      if (class=="f2") {e.het <- 0.5}
      if (class=="bc1") {e.hi <- 0.25}
      if (class=="bc1") {e.het <- 0.5}
      if (class=="bc2") {e.hi <- 0.75}
      if (class=="bc2") {e.het <- 0.5}
      #observed values
      o.hi <- mean(temp[temp$pop==class,]$hybrid.index)
      o.het <- mean(temp[temp$pop==class,]$heterozygosity)
      #accuracy
      a.hi <- 1 - abs(e.hi-o.hi)/e.hi
      a.het <- 1 - abs(e.het-o.het)/e.het
      #precision
      p.hi <- mean(abs(o.hi - temp[temp$pop==class,]$hybrid.index))
      p.het <- mean(abs(o.het - temp[temp$pop==class,]$heterozygosity))
      #put in dataframe
      nAIMs.ap <- rbind(nAIMs.ap, c(m, n, class, "hybrid index", e.hi, o.hi, a.hi, p.hi))
      nAIMs.ap <- rbind(nAIMs.ap, c(m, n, class, "heterozygosity", e.het, o.het, a.het, p.het))
    }
    tri <- triangle.plot(temp, colors = c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837"), alpha = 0.8)
    tri <- tri + annotate(geom="text", x=-0.05, y=1.05, hjust = 0, label=paste0(d, " SNPs"), size=4) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_x_continuous(breaks = c(0,0.5,1)) +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
    assign(paste0("tri.",n),tri)
    z <- z + 1
  }
}


D750.0.5.few.AIMs <- ggarrange(tri.10, tri.20, tri.40, tri.60, tri.80, tri.100,
                                 labels = c("A", "B", "C", "D", "E", "F"),
                                 font.label = list(size = 14),
                                 #label.y = -0.2,
                                 label.x = 0.1,
                                 ncol = 3, nrow = 2)

ggsave("SFig1fromR.pdf", plot = D750.0.5.few.AIMs, path = "Figures/", width = 7, height = 4)

colnames(nAIMs.ap) <- c("divergence", "num.AIMs", "class", "metric", "expected", "observed", "accurary", "precision")
nAIMs.ap$divergence <- rep("Low", 48)
nAIMs.ap$accurary <- as.numeric(nAIMs.ap$accurary)
nAIMs.ap$precision <- as.numeric(nAIMs.ap$precision)

nAIMs.a.plot <- ggplot(data = nAIMs.ap,  aes(fill=factor(num.AIMs, levels = c("10", "20", "40", "60", "80", "100")), y=accurary, x=factor(metric, levels = c("hybrid index", "heterozygosity")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#08306b", "#08519c", "#2171b5", "#4292c6", "#6baed6", "#9ecae1")) +
  ggtitle("") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "Number of AIMs") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))


nAIMs.p.plot <- ggplot(data = nAIMs.ap,  aes(fill=factor(num.AIMs, levels = c("10", "20", "40", "60", "80", "100")), y=precision, x=factor(metric, levels = c("hybrid index", "heterozygosity")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ab5d", "#74c476", "#a1d99b")) +
  ggtitle("") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.14), breaks = c(0, 0.07, 0.14)) +
  #theme(legend.position = "none")
  labs(fill = "Number of AIMs") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))


nAIMs.ap.plot <- ggarrange(nAIMs.a.plot, nAIMs.p.plot,
                                 labels = c("A", "B"),
                                 font.label = list(size = 14),
                                 #label.y = -0.2,
                                 label.x = 0.1,
                                 ncol = 1, nrow = 2)

ggsave("FigureS8fromR.pdf", plot = nAIMs.ap.plot, path = "Figures/", width = 8, height = 8)





################################################
#  Mean absolute error of triangulaR and bgchm #
################################################



# Average error for hi and het estimates
triangulaR_error <- data.frame(matrix(nrow = 0, ncol = 6))
triangulaR_hihet <- data.frame(matrix(nrow = 0, ncol = 7))

for (m in c("D750", "D1000", "D2000")) {
  for (i in c("even_5", "even_20")) {
    set.seed(312564)
    if (i == "even_20") {
      j <- c(20,20)
    }
    if (i == "even_5") {
      j <- c(5,5)
    }
    if (i == "uneven") {
      j <- c(20,5)
    }
    print(j)
    pmt <- get(paste0(m, ".pm"))
    ids <- c(sample(pmt[pmt$pop=="p0","id"], j[1] + 20, replace = F), sample(pmt[pmt$pop=="p20","id"], j[2] + 20, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
    pmt <- pmt[pmt$id %in% ids,]
    pmt[pmt$pop=="p0","pop"] <- c(rep("p0",j[1]), rep("P1",20))
    pmt[pmt$pop=="p20","pop"] <- c(rep("p20",j[2]), rep("P2",20))
    temp <- get(paste0(m, ".manual"))[,pmt$id]
    temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = 1)
    d <- nrow(temp[1][[1]])
    hihet <- hybridIndex.matrix(matrix = temp[1][[1]], pm = pmt, p1 = "p0", p2 = "p20")
    
    for (class in c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2")) {
    
    if (class == "F1") {
      syn <- c("F1", "f1")
      hi.mean <- 0.5
      het.mean <- 1
    }
    if (class == "F2") {
      syn <- c("F2", "f2")
      hi.mean <- 0.5
      het.mean <- 0.5
    }
    if (class == "P1xF1") {
      syn <- c("P1xF1", "bc1")
      hi.mean <- 0.25
      het.mean <- 0.5
    }
    if (class == "P2xF1") {
      syn <- c("P2xF1", "bc2")
      hi.mean <- 0.75
      het.mean <- 0.5
    }
    if (class == "P1") {
      syn <- c("P1")
      hi.mean <- 0
      het.mean <- 0
    }
    if (class == "P2") {
      syn <- c("P2")
      hi.mean <- 1
      het.mean <- 0
    }
      
      hi.error <- abs(hihet[hihet$pop %in% syn,]$hybrid.index - hi.mean)
      het.error <- abs(hihet[hihet$pop %in% syn,]$heterozygosity - het.mean)
      
      triangulaR_error <- rbind(triangulaR_error, cbind(hi.error, het.error, class, "triangulaR", m, i))
      
      triangulaR_hihet <- rbind(triangulaR_hihet, cbind(hihet, m, i))
    }
    triangulaR_error$class <- factor(triangulaR_error$class, levels = c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2"))
    print(triangle.plot(hihet))
    
  }
}

test <- data.frame(hi.even5 = triangulaR_hihet[triangulaR_hihet$i=="even_5" & !triangulaR_hihet$pop %in% c("p0", "20") & triangulaR_hihet$m %in% c("D750"),]$hybrid.index,
           hi.uneven = triangulaR_hihet[triangulaR_hihet$i=="uneven" & !triangulaR_hihet$pop %in% c("p0", "20") & triangulaR_hihet$m %in% c("D750"),]$hybrid.index)
ggplot(test) +
  geom_point(aes(x=hi.even5, y=hi.uneven)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black")




ggplot(triangulaR_error, aes(x = class, y = as.numeric(hi.error), fill = i)) +
  geom_boxplot(aes(color = i), position = position_dodge(width = 0.75), alpha = 0) +
  geom_jitter(
    aes(color = i), 
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
    alpha = 1, size = 3, show.legend = FALSE
  ) +
  facet_wrap(~factor(m, levels = c("D750", "D1000", "D2000"))) +
  #scale_fill_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ab5d", "#74c476", "#a1d99b")) +
  ggtitle("") +
  xlab(paste("")) +
  ylab(paste("Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  #scale_y_continuous(limits = c(0, 0.14), breaks = c(0, 0.07, 0.14)) +
  #theme(legend.position = "none")
  #labs(fill = "Number of AIMs") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1))


ggplot(triangulaR_error, aes(x = class, y = as.numeric(het.error), fill = i)) +
  geom_boxplot(aes(color = i), position = position_dodge(width = 0.75), alpha = 0) +
  geom_jitter(
    aes(color = i), 
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
    alpha = 1, size = 3, show.legend = FALSE
  ) +
  facet_wrap(~factor(m, levels = c("D750", "D1000", "D2000"))) +
  #scale_fill_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ab5d", "#74c476", "#a1d99b")) +
  ggtitle("") +
  xlab(paste("")) +
  ylab(paste("Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  #scale_y_continuous(limits = c(0, 0.14), breaks = c(0, 0.07, 0.14)) +
  #theme(legend.position = "none")
  #labs(fill = "Number of AIMs") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1))


bgchm_error <- data.frame(matrix(nrow = 0, ncol = 6))
bgchm_hihet <- data.frame(matrix(nrow = 0, ncol = 7))



for (m in c("D750", "D1000", "D2000")) {
  for (i in c("even_5", "even_20")) {
    set.seed(312564)
    if (i == "even_20") {
      j <- c(20,20)
    }
    if (i == "even_5") {
      j <- c(5,5)
    }
    if (i == "uneven") {
      j <- c(20,5)
    }
    print(j)
    pmt <- get(paste0(m, ".pm"))
    ids <- c(sample(pmt[pmt$pop=="p0","id"], j[1] + 20, replace = F), sample(pmt[pmt$pop=="p20","id"], j[2] + 20, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
    pmt <- pmt[pmt$id %in% ids,]
    pmt[pmt$pop=="p0","pop"] <- c(rep("p0",j[1]), rep("P1",20))
    pmt[pmt$pop=="p20","pop"] <- c(rep("p20",j[2]), rep("P2",20))
    temp <- get(paste0(m, ".manual"))[,pmt$id]
    temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = 1)

    temp <- t(temp[1][[1]])
    identical(row.names(temp), pmt$id)
    temp <- apply(temp, 2, as.numeric)
    row.names(temp) <- pmt$id
    
    example.P1 <- temp[pmt[pmt$pop=="p0",]$id,]
    example.P2 <- temp[pmt[pmt$pop=="p20",]$id,]
    
    example_q_out<-est_Q(Gx=temp,G0=example.P1,G1=example.P2,model="genotype",ploidy="diploid")
    tri_plot(hi=example_q_out$hi[,1],Q10=example_q_out$Q10[,1],pdf=FALSE,pch=19)
    segments(example_q_out$hi[,1],example_q_out$Q10[,2],example_q_out$hi[,1],example_q_out$Q10[,3])
    segments(example_q_out$hi[,2],example_q_out$Q10[,1],example_q_out$hi[,3],example_q_out$Q10[,1])
    
    
    hihet <- data.frame(cbind(pmt$id, pmt$pop, example_q_out$hi[,1],example_q_out$Q10[,1]))
    colnames(hihet) <- c("id", "pop", "hybrid.index", "heterozygosity")
    hihet$hybrid.index <- as.numeric(hihet$hybrid.index)
    hihet$heterozygosity <- as.numeric(hihet$heterozygosity)
    
    for (class in c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2")) {
      if (class == "F1") {
        syn <- c("F1", "f1")
        hi.mean <- 0.5
        het.mean <- 1
      }
      if (class == "F2") {
        syn <- c("F2", "f2")
        hi.mean <- 0.5
        het.mean <- 0.5
      }
      if (class == "P1xF1") {
        syn <- c("P1xF1", "bc1")
        hi.mean <- 0.25
        het.mean <- 0.5
      }
      if (class == "P2xF1") {
        syn <- c("P2xF1", "bc2")
        hi.mean <- 0.75
        het.mean <- 0.5
      }
      if (class == "P1") {
        syn <- c("P1")
        hi.mean <- 0
        het.mean <- 0
      }
      if (class == "P2") {
        syn <- c("P2")
        hi.mean <- 1
        het.mean <- 0
      }
  
      hi.error <- abs(hihet[hihet$pop %in% syn,]$hybrid.index - hi.mean)
      het.error <- abs(hihet[hihet$pop %in% syn,]$heterozygosity - het.mean)
      
      bgchm_error <- rbind(bgchm_error, cbind(hi.error, het.error, class, "bgchm", m, i))
      bgchm_hihet <- rbind(bgchm_hihet, cbind(hihet, m, i))
      
    }
    bgchm_error$class <- factor(bgchm_error$class, levels = c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2"))
    


  }

}

all_error <- rbind(triangulaR_error, bgchm_error)
all_error$box <- paste(all_error$V4, all_error$i, sep = "_")
all_error$box <- factor(all_error$box, levels = c("bgchm_even_5", "triangulaR_even_5", "bgchm_even_20", "triangulaR_even_20"))




hi.error <- ggplot(all_error, aes(x = class, y = as.numeric(hi.error), fill = box)) +
  geom_boxplot(aes(color = box), width = 0.5,
               position = position_dodge(width = 0.75), alpha = 1, outlier.shape = 16) +
  #geom_jitter(
  #  aes(color = box), 
  #  position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
  #  alpha = 1, size = 2, show.legend = FALSE
  #) +
  facet_wrap(~factor(m, levels = c("D750", "D1000", "D2000"))) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.16), breaks = c(0, 0.05, 0.1, 0.15)) +
  scale_color_manual(values=c("#b2182b", "#2166ac", "#f4a582", "#92c5de")) +
  scale_fill_manual(values=c("#b2182b", "#2166ac", "#f4a582", "#92c5de")) +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1))

het.error <- ggplot(all_error, aes(x = class, y = as.numeric(het.error), fill = box)) +
  geom_boxplot(aes(fill = box, color = box), width = 0.5,
               position = position_dodge(width = 0.7), alpha = 1, outlier.shape = 16) +
  #geom_jitter(
  #  aes(color = box), 
  #  position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3), 
  #  alpha = 1, size = 2, show.legend = FALSE
  #) +
  facet_wrap(~factor(m, levels = c("D750", "D1000", "D2000"))) +
  ggtitle("Interclass Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.28), breaks = c(0, 0.10, 0.20)) +
  scale_color_manual(values=c("#b2182b", "#2166ac", "#f4a582", "#92c5de")) +
  scale_fill_manual(values=c("#b2182b", "#2166ac", "#f4a582", "#92c5de")) +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1))

tri_bgchm_error <-ggarrange(hi.error, het.error,
                            ncol = 1, nrow = 2)

ggsave("Figure4newfromR.pdf", plot = tri_bgchm_error, path = "Figures/", width = 15, height = 8)



colnames(bgchm_hihet) <- paste("bgchm", colnames(bgchm_hihet), sep = ".")
colnames(triangulaR_hihet) <- paste("triangulaR", colnames(triangulaR_hihet), sep = ".")
all_hihet <- cbind(bgchm_hihet, triangulaR_hihet)


ggplot(all_hihet, aes(x=as.numeric(bgchm.hybrid.index), y=as.numeric(triangulaR.hybrid.index))) +
  geom_point(aes(color = bgchm.i)) +
  facet_wrap(~factor(bgchm.m, levels = c("D750", "D1000", "D2000"))) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
  ggtitle("") +
  xlab(paste("bgchm")) +
  ylab(paste("triangulaR")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggplot(all_hihet, aes(x=as.numeric(bgchm.heterozygosity), y=as.numeric(triangulaR.heterozygosity))) +
  geom_point(aes(color = bgchm.i)) +
  facet_wrap(~factor(bgchm.m, levels = c("D750", "D1000", "D2000"))) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
  ggtitle("") +
  xlab(paste("bgchm")) +
  ylab(paste("triangulaR")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))






# 5 individuals in parental pops, 3 different afdts

for (m in c("D750", "D1000", "D2000")) {
  set.seed(7825061)
  pmt <- get(paste0(m, ".pm"))
  ids <- c(sample(pmt[pmt$pop=="p0","id"], 5 + 20, replace = F), sample(pmt[pmt$pop=="p20","id"], 5 + 20, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
  pmt <- pmt[pmt$id %in% ids,]
  pmt[pmt$pop=="p0","pop"] <- c(rep("p0",5), rep("P1",20))
  pmt[pmt$pop=="p20","pop"] <- c(rep("p20",5), rep("P2",20))
  temp <- get(paste0(m, ".manual"))[,pmt$id]
  af.diffs <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = c(1,0.75,0.5))

  for (a in 1:3) {
    d <- nrow(af.diffs[a][[1]])
    temp <- hybridIndex.matrix(matrix = af.diffs[a][[1]], pm = pmt, p1 = "p0", p2 = "p20")
    tri <- triangle.plot(temp, colors = c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "black", "#762a83", "#1b7837", "black"), alpha = 0.8)
    tri <- tri + annotate(geom="text", x=-0.05, y=0.95, hjust = 0, label=paste0(d, " SNPs"), size=4) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_x_continuous(breaks = c(0,0.5,1)) +
      theme(legend.position = "none", 
            axis.title.y = element_text(size = 14), 
            axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
    
    assign(paste0(m, ".20.afdts.", a), tri)
  }
}

afdts.20 <- ggarrange(D750.20.afdts.1, D1000.20.afdts.1, D2000.20.afdts.1, 
                   D750.20.afdts.2, D1000.20.afdts.2, D2000.20.afdts.2,
                   D750.20.afdts.3, D1000.20.afdts.3, D2000.20.afdts.3,
                   #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                   #font.label = list(size = 14),
                   #label.y = -0.2,
                   #label.x = 0.2,
                   ncol = 3, nrow = 3)

ggsave("FigureXXXXXXfromR.pdf", plot = afdts.20, path = "Figures/", width = 7, height = 6)











# bgchm 5 parentals 0.5 d
bgchm_error_0.5 <- data.frame(matrix(nrow = 0, ncol = 6))
bgchm_hihet_0.5 <- data.frame(matrix(nrow = 0, ncol = 7))

for (m in c("D750", "D1000", "D2000")) {
  for (i in c(5)) {
    set.seed(7825061)
    pmt <- get(paste0(m, ".pm"))
    ids <- c(sample(pmt[pmt$pop=="p0","id"], i + 20, replace = F), sample(pmt[pmt$pop=="p20","id"], i + 20, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
    pmt <- pmt[pmt$id %in% ids,]
    pmt[pmt$pop=="p0","pop"] <- c(rep("p0",i), rep("P1",20))
    pmt[pmt$pop=="p20","pop"] <- c(rep("p20",i), rep("P2",20))
    temp <- get(paste0(m, ".manual"))[,pmt$id]
    temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = 0.5)
    
    temp <- t(temp[1][[1]])
    identical(row.names(temp), pmt$id)
    temp <- apply(temp, 2, as.numeric)
    row.names(temp) <- pmt$id
    
    example.P1 <- temp[pmt[pmt$pop=="p0",]$id,]
    example.P2 <- temp[pmt[pmt$pop=="p20",]$id,]
    
    print(Sys.time())
    example_q_out<-est_Q(Gx=temp,G0=example.P1,G1=example.P2,model="genotype",ploidy="diploid")
    print(Sys.time())
    tri_plot(hi=example_q_out$hi[,1],Q10=example_q_out$Q10[,1],pdf=FALSE,pch=19)
    segments(example_q_out$hi[,1],example_q_out$Q10[,2],example_q_out$hi[,1],example_q_out$Q10[,3])
    segments(example_q_out$hi[,2],example_q_out$Q10[,1],example_q_out$hi[,3],example_q_out$Q10[,1])
    title(m)
    
    hihet <- data.frame(cbind(pmt$id, pmt$pop, example_q_out$hi[,1],example_q_out$Q10[,1]))
    colnames(hihet) <- c("id", "pop", "hybrid.index", "heterozygosity")
    hihet$hybrid.index <- as.numeric(hihet$hybrid.index)
    hihet$heterozygosity <- as.numeric(hihet$heterozygosity)
    
    for (class in c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2")) {
      if (class == "F1") {
        syn <- c("F1", "f1")
        hi.mean <- 0.5
        het.mean <- 1
      }
      if (class == "F2") {
        syn <- c("F2", "f2")
        hi.mean <- 0.5
        het.mean <- 0.5
      }
      if (class == "P1xF1") {
        syn <- c("P1xF1", "bc1")
        hi.mean <- 0.25
        het.mean <- 0.5
      }
      if (class == "P2xF1") {
        syn <- c("P2xF1", "bc2")
        hi.mean <- 0.75
        het.mean <- 0.5
      }
      if (class == "P1") {
        syn <- c("P1")
        hi.mean <- 0
        het.mean <- 0
      }
      if (class == "P2") {
        syn <- c("P2")
        hi.mean <- 1
        het.mean <- 0
      }
      
      hi.error <- abs(hihet[hihet$pop %in% syn,]$hybrid.index - hi.mean)
      het.error <- abs(hihet[hihet$pop %in% syn,]$heterozygosity - het.mean)
      
      bgchm_error_0.5 <- rbind(bgchm_error_0.5, cbind(hi.error, het.error, class, "bgchm", m, i))
      bgchm_hihet_0.5 <- rbind(bgchm_hihet_0.5, cbind(hihet, m, i))
      
    }
    bgchm_error_0.5$class <- factor(bgchm_error_0.5$class, levels = c("F1", "F2", "P1xF1", "P2xF1", "P1", "P2"))
    
    
    
  }
  
}


###########################################
###   Time for triangulaR and bgchm   #####
###########################################
tri_runtime <- function(m=NULL, i=NULL, d=NULL, extra.parentals=NULL, expand.factor=NULL) {
# m is simulation, i is number of parentals to use as reference, d is allele freq difference
# extra.parentals is number to take and not assign to reference, just for testing 160 vs 320 total inds
# expand.factor is whether or not to duplicate sites, just to increase the number of sites
  pmt <- get(paste0(m, ".pm"))
  ids <- c(pmt[pmt$pop=="p0","id"][1:i], pmt[pmt$pop=="p0","id"][(i+1):(i+extra.parentals)], pmt[pmt$pop=="p20","id"][1:i], pmt[pmt$pop=="p20","id"][(i+1):(i+extra.parentals)],  pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])

  #ids <- c(sample(pmt[pmt$pop=="p0","id"], i, replace = F), sample(pmt[pmt$pop=="p20","id"], i, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])
  #ids <- c(ids, sample(pmt[pmt$pop=="p0","id"], extra.parentals, replace = F), sample(pmt[pmt$pop=="p20","id"], extra.parentals, replace = F))
  pmt <- pmt[pmt$id %in% ids,]
  pmt[pmt$pop=="p0","pop"] <- c(rep("p0",i), rep("P1",extra.parentals))
  pmt[pmt$pop=="p20","pop"] <- c(rep("p20",i), rep("P2",extra.parentals))
  
  
  temp <- get(paste0(m, ".manual"))[,pmt$id]
  temp <- alleleFreqDiff.matrix(matrix = temp, pm = pmt, p1 = "p0", p2 = "p20", differences = d)
  
  temp <- t(temp[1][[1]])
  identical(row.names(temp), pmt$id)
  temp <- apply(temp, 2, as.numeric)
  row.names(temp) <- pmt$id
  
  if (expand.factor) {
    temp <- cbind(temp, temp, temp, temp, temp, temp)
  }
  
  # number of inds
  print("inds:")
  print(nrow(temp))
  # number of sites
  print("sites:")
  print(ncol(temp))
  
  example.P1 <- temp[pmt[pmt$pop=="p0",]$id,]
  example.P2 <- temp[pmt[pmt$pop=="p20",]$id,]
  
  print("triangulaR")
  print(Sys.time())
  print(system.time(hihet <- hybridIndex.matrix(t(temp), pmt, "p0", "p20")))
  print(Sys.time())
  print(triangle.plot(hihet, jitter = 0.002))
  
  
  print("bgchm")
  print(Sys.time())
  print(system.time(example_q_out<-est_Q(Gx=temp,G0=example.P1,G1=example.P2,model="genotype",ploidy="diploid")))
  print(Sys.time())
  
  hihet <- data.frame(cbind(pmt$id, pmt$pop, example_q_out$hi[,1],example_q_out$Q10[,1]))
  colnames(hihet) <- c("id", "pop", "hybrid.index", "heterozygosity")
  hihet$hybrid.index <- as.numeric(hihet$hybrid.index)
  hihet$heterozygosity <- as.numeric(hihet$heterozygosity)
  print(triangle.plot(hihet, jitter = 0.002))
}


# 160 inds 346 AIMs
tri_runtime(m = "D750", i = 20, d = 0.5, extra.parentals = 20, expand.factor = F)

# 160 inds 2076 AIMs
tri_runtime(m = "D750", i = 20, d = 0.5, extra.parentals = 20, expand.factor = T)

# 360 inds 346 AIMs
tri_runtime(m = "D750", i = 20, d = 0.5, extra.parentals = 120, expand.factor = F)

# 360 inds ~ 2076 AIMs
tri_runtime(m = "D750", i = 20, d = 0.5, extra.parentals = 120, expand.factor = T)








########################################
###            Low depth             ###
########################################

depth_simulator <- function(m=NULL, i=NULL, d=NULL, depth=NULL) {
  # m is simulation, i is number of parentals to use as reference, d is allele freq difference
  pmt <- get(paste0(m, ".pm"))
  ids <- c(sample(pmt[pmt$pop=="p0","id"], i + 20, replace = F), sample(pmt[pmt$pop=="p20","id"], i + 20, replace = F), pmt[pmt$pop %in% (c("f1", "f2", "bc1", "bc2")),"id"])

  pmt <- pmt[pmt$id %in% ids,]
  pmt[pmt$pop=="p0","pop"] <- c(rep("p0",i), rep("P1",20))
  pmt[pmt$pop=="p20","pop"] <- c(rep("p20",i), rep("P2",20))
  
  
  temp <- get(paste0(m, ".manual"))[,pmt$id]
  r.names <- row.names(temp)
  
  temp <- apply(temp, 2, as.numeric)
  temp2 <- temp
  
  
  for (i in 1:nrow(temp)) {
    for (j in 1:ncol(temp)) {
      if (temp[i,j]==1) {
        alleles <- rbinom(depth, 1, 0.5)
        het <- all(c(0 %in% alleles, 1 %in% alleles))
        if (het) {
          next
        }
        if (all(alleles==0)) {
          temp2[i,j] <- 0
        } else if(all(alleles==1)) {
          temp2[i,j] <- 2
        } else {
          stop("")
        }
      }
    }
  }
  
  rownames(temp2) <- r.names
  temp2 <- alleleFreqDiff.matrix(matrix = temp2, pm = pmt, p1 = "p0", p2 = "p20", differences = d)
  
  temp2 <- temp2[1][[1]]

  # number of inds
  print("inds:")
  print(ncol(temp2))
  # number of sites
  print("sites:")
  print(nrow(temp2))
  
  
  hihet <- hybridIndex.matrix(temp2, pmt, "p0", "p20")
  print(triangle.plot(hihet, jitter = 0.002))
  return(hihet)
  
}

set.seed(123841)
depth_simulator(m="D2000", i=20, d=1, depth=6)





# Accurary and precision 
for (m in c("D750", "D1000", "D2000")) {
  set.seed(6412)
  if (m == "D750") {depth.ap <- data.frame(matrix(nrow = 0, ncol = 8))}
  for (depth in c(2, 4, 6, 8, 10)) {
    temp <- depth_simulator(m=m, i=20, d=1, depth=depth)
    
    for (class in c("f1", "f2", "bc1", "bc2")) {
      #expected values
      if (class=="f1") {e.hi <- 0.5}
      if (class=="f1") {e.het <- 1}
      if (class=="f2") {e.hi <- 0.5}
      if (class=="f2") {e.het <- 0.5}
      if (class=="bc1") {e.hi <- 0.25}
      if (class=="bc1") {e.het <- 0.5}
      if (class=="bc2") {e.hi <- 0.75}
      if (class=="bc2") {e.het <- 0.5}
      #observed values
      o.hi <- mean(temp[temp$pop==class,]$hybrid.index)
      o.het <- mean(temp[temp$pop==class,]$heterozygosity)
      #accuracy
      a.hi <- 1 - abs(e.hi-o.hi)/e.hi
      a.het <- 1 - abs(e.het-o.het)/e.het
      #precision
      p.hi <- mean(abs(o.hi - temp[temp$pop==class,]$hybrid.index))
      p.het <- mean(abs(o.het - temp[temp$pop==class,]$heterozygosity))
      #put in dataframe
      depth.ap <- rbind(depth.ap, c(m, depth, class, "hybrid index", e.hi, o.hi, a.hi, p.hi))
      depth.ap <- rbind(depth.ap, c(m, depth, class, "heterozygosity", e.het, o.het, a.het, p.het))
    }
  }
}
colnames(depth.ap) <- c("divergence", "depth", "class", "metric", "expected", "observed", "accuracy", "precision")
depth.ap$divergence <- c(rep("Low", 40), rep("Med", 40), rep("High", 40))
depth.ap$accuracy <- as.numeric(depth.ap$accuracy)
depth.ap$precision <- as.numeric(depth.ap$precision)

hi.depth.a <- ggplot(data = depth.ap[depth.ap$metric=="hybrid index",],  aes(fill=factor(depth, levels = c("10", "8", "6", "4", "2")), y=accuracy, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "Depth") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.depth.a <- ggplot(data = depth.ap[depth.ap$metric=="heterozygosity",],  aes(fill=factor(depth, levels = c("10", "8", "6", "4", "2")), y=accuracy, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1")) +
  ggtitle("Interclass Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Accuracy")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  #theme(legend.position = "none")
  labs(fill = "Depth") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))


hi.depth.p <- ggplot(data = depth.ap[depth.ap$metric=="hybrid index",],  aes(fill=factor(depth, levels = c("10", "8", "6", "4", "2")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d", "#74c476", "#a1d99b")) +
  ggtitle("Hybrid Index") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.12), breaks = c(0, 0.06, 0.12)) +
  #theme(legend.position = "none")
  labs(fill = "Depth") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

het.depth.p <- ggplot(data = depth.ap[depth.ap$metric=="heterozygosity",],  aes(fill=factor(depth, levels = c("10", "8", "6", "4", "2")), y=precision, x=factor(divergence, levels = c("Low", "Med", "High")))) +
  geom_bar(alpha=1, position = "dodge", stat = "identity") +
  facet_wrap(~factor(class, levels = c("f1", "f2", "bc1", "bc2"))) +
  scale_fill_manual(values=c("#005a32", "#238b45", "#41ab5d", "#74c476", "#a1d99b")) +
  ggtitle("Interclass Heterozygosity") +
  xlab(paste("")) +
  ylab(paste("Precision")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 0.12), breaks = c(0, 0.06, 0.12)) +
  #theme(legend.position = "none")
  labs(fill = "Depth") +
  theme(panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    size = 1))

depth.ap.plot <- ggarrange(hi.depth.a, het.depth.a, hi.depth.p, het.depth.p,
                                 labels = c("A", "B", "C", "D"),
                                 font.label = list(size = 14),
                                 #label.y = -0.2,
                                 label.x = 0.1,
                                 ncol = 2, nrow = 2)

ggsave("FigureSXXfromR.pdf", plot = depth.ap.plot, path = "Figures/", width = 10, height = 8)

min(depth.ap[depth.ap$metric=="hybrid index" & depth.ap$depth==2, "accuracy"])
mean(depth.ap[depth.ap$metric=="hybrid index" & depth.ap$depth==2, "accuracy"])

mean(depth.ap[depth.ap$metric=="heterozygosity" & depth.ap$depth==6, "accuracy"])


# average allele frequency difference between p0 and p20 in each simulation
D750.afd <- specFreqDiff(matrix = D750.manual, pm = D750.pm, p1 = "p0", p2 = "p20", difference = 0)
D1000.afd <- specFreqDiff(matrix = D1000.manual, pm = D1000.pm, p1 = "p0", p2 = "p20", difference = 0)
D2000.afd <- specFreqDiff(matrix = D2000.manual, pm = D2000.pm, p1 = "p0", p2 = "p20", difference = 0)

round(mean(D750.afd$diff), digits = 2)
round(mean(D1000.afd$diff), digits = 2)
round(mean(D2000.afd$diff), digits = 2)





# mean and sd of pi from sims
pis <- c(0.0563105370493329, 0.0580125236491738, 0.0561969974346078, 0.0553544930338812, 0.0553514228744697, 0.0566716867091329)
round(mean(pis), digits = 3)
round(sd(pis), digits = 3)









