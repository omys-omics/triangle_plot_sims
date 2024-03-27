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


# Accurary and precision 
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





defFst <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  af <- specFreqDiff(matrix, pm = pm, p1 = p1, p2 = p2, difference = 0)
  
  # calculate average frequency of p1 allele
  o_p1_bar <- rowMeans(af[,c("p1","p2")])
  # calculate average frequency of p2 allele
  o_p2_bar <- 1-o_p1_bar
  # calculate variance
  o_var <- 0
  o_var <- o_var + (af$p1-o_p1_bar)^2
  o_var <- o_var + (af$p2-o_p1_bar)^2
  o_var <- o_var/2
  
  # calculate Fst
  obs.fst <- o_var/(o_p1_bar*o_p2_bar)
  hist(obs.fst)
  mean.obs.fst <- mean(obs.fst, na.rm = T)
  return(mean.obs.fst)
}

defFst(D750.manual, pm = D750.pm, p1 = "p0", p2 = "p20")
defFst(D1000.manual, pm = D1000.pm, p1 = "p0", p2 = "p20")
defFst(D2000.manual, pm = D2000.pm, p1 = "p0", p2 = "p20")




# pi for p0 and p20 for each simulation
pixy_pis <- c(0.056310537,0.058012524,0.056196997,0.055354493,0.055351423,0.056671687)
mean(pixy_pis)
sd(pixy_pis)

# average allele frequency difference between p0 and p20 in each simulation
D750.afd <- specFreqDiff(matrix = D750.manual, pm = D750.pm, p1 = "p0", p2 = "p20", difference = 0)
D1000.afd <- specFreqDiff(matrix = D1000.manual, pm = D1000.pm, p1 = "p0", p2 = "p20", difference = 0)
D2000.afd <- specFreqDiff(matrix = D2000.manual, pm = D2000.pm, p1 = "p0", p2 = "p20", difference = 0)

round(mean(D750.afd$diff), digits = 2)
round(mean(D1000.afd$diff), digits = 2)
round(mean(D2000.afd$diff), digits = 2)



########################
#     Calculate Fst    #
########################


vcfR <- read.vcfR("D750/D750.end.of.div.vcf")
vcfR <- vcfR[samples = c("i0", "i1", "i2", "i3", "i2241", "i2242", "i2243", "i2244")]
gi <- vcfR2genind(vcfR)
gi@pop <- as.factor(c("p0","p0","p0","p0","p20","p20","p20","p20"))
df <- genind2df(gi)
df2 <- data.frame(sapply(df, as.numeric))
df2[df2==0] <- 22
df2[df2==10] <- 23
df2[df2==1] <- 23
df2[df2==11] <- 33
df2[,1] <- c(1,1,1,1,2,2,2,2)
pairwise.WCfst(gi)
pairwise.WCfst(df2)


wrightFst <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  af <- specFreqDiff(matrix, pm = pm, p1 = p1, p2 = p2, difference = 0)
  
  # calculate expected heterozygosity under HWE for each pop separately
  af$p1.He <- 2*(af$p1 * (1-af$p1))
  af$p2.He <- 2*(af$p2 * (1-af$p2))
  
  
  # calculate expected heterozygosity under HWE if everything is one population
  af$ave.freq <- ((af$p1 * sum(pm$pop==p1)) + (af$p2 * sum(pm$pop==p2))) / (sum(pm$pop==p1) + sum(pm$pop==p2))
  af$Ht <- 2*(af$ave.freq * (1-af$ave.freq))
  
  # calulate fst by averaging across SNPs in two different ways
  # ratio of averages is considered better
  ratio.of.aves <- 1-mean(colMeans(af)[4:5])/colMeans(af[7])
  ave.of.ratios <- mean(1-((af$p1.He + af$p2.He) / 2) / af$Ht)
  fst <- c(ratio.of.aves, ave.of.ratios)
  names(fst) <- c("ratio.of.aves", "ave.of.ratios")
  return(fst)
}

wrightFst(D750.manual, pm = D750.pm, p1 = "p0", p2 = "p20")
wrightFst(D1000.manual, pm = D1000.pm, p1 = "p0", p2 = "p20")
wrightFst(D2000.manual, pm = D2000.pm, p1 = "p0", p2 = "p20")

defFst <- function(matrix = NULL, pm = NULL, p1 = NULL, p2 = NULL) {
  af <- specFreqDiff(matrix, pm = pm, p1 = p1, p2 = p2, difference = 0)
  
  # calculate average frequency of p1 allele
  o_p1_bar <- rowMeans(af[,c("p1","p2")])
  # calculate average frequency of p2 allele
  o_p2_bar <- 1-o_p1_bar
  # calculate variance
  o_var <- 0
  o_var <- o_var + (af$p1-o_p1_bar)^2
  o_var <- o_var + (af$p2-o_p1_bar)^2
  o_var <- o_var/2
  
  # calculate Fst
  obs.fst <- o_var/(o_p1_bar*o_p2_bar)
  hist(obs.fst)
  mean.obs.fst <- mean(obs.fst, na.rm = T)
  return(mean.obs.fst)
}

defFst(D750.manual, pm = D750.pm, p1 = "p0", p2 = "p20")
defFst(D1000.manual, pm = D1000.pm, p1 = "p0", p2 = "p20")
defFst(D2000.manual, pm = D2000.pm, p1 = "p0", p2 = "p20")









vcfR <- read.vcfR("D750/D750.end.of.div.invariants.vcf")
test <- extract.gt(test)

vcfR <- read.vcfR("D750/D750.end.of.div.vcf")
vcfR <- vcfR[samples = c("i0", "i1", "i2", "i3", "i2245", "i2246", "i2247", "i2248")]
#vcfR <- min.mac(vcfR, min.mac = 1)
gt <- extract.gt(vcfR)
gt <- t(gt)
lnames <- colnames(gt)
gt <- data.frame(gt)
colnames(gt) <- lnames

gst <- genetic_diff(vcfR, as.factor(c("p0","p0","p0","p0","p20","p20","p20","p20")), method = "nei")
mean(gst$Gst, na.rm = T)


gi3 <- df2genind(gt, sep = "\\|", ncode = 1, pop = as.factor(c("p0","p0","p0","p0","p20","p20","p20","p20")),
          type = "codom")

gl <- gi2gl(gi3)
gl2vcf(gl)

df3 <- genind2df(gi3)
identical(gi3, df3)


gen <- vcfR2genlight(vcfR)
gen@pop <- as.factor(sort(rep(c(0,20), 1000)))

gi <- vcfR2genind(vcfR)
gi@pop <- as.factor(c("p0","p0","p0","p0","p20","p20","p20","p20"))
gi@pop <- as.factor(c(rep(0,sum(D750.pm$pop=="p0")), rep (20,sum(D750.pm$pop=="p20"))))
df <- genind2df(gi)
df2 <- data.frame(sapply(df, as.numeric))
df <- data.frame(sapply(df, as.numeric))
df2[df2==0] <- 22
df2[df2==10] <- 23
df2[df2==1] <- 23
df2[df2==11] <- 33
df2[,1] <- c(1,1,1,1,2,2,2,2)

pairwise.neifst(gi)
pairwise.WCfst(gi)
pairwise.WCfst(df)
pairwise.WCfst(df2)
pairwise.WCfst(gi3)
pairwise.WCfst(df[c(1:100,1901:2000),])

wc.output <- wc(df2)
wc.output$per.loc
colMeans(wc.output$per.loc, na.rm = T)
wc.output$FST


mean(hierfstat::Hs(gi))
mean(adegenet::Hs(gi))

gi.no.pop <- gi
gi.no.pop@pop <- as.factor(rep(1, 8))
Hs(gi.no.pop)[1]

1 - (mean(Hs(gi)) / Hs(gi.no.pop)[1])

gen@pop <- as.factor(sort(rep(c(0,20), 4)))

heat <- stamppFst(gen)

#extract the pairwise matrix
m<-heat$Fsts
#fill in upper triangle of the matrix
m[upper.tri(m)] <- t(m)[upper.tri(m)]

#melt to tidy format for ggplotting
heat <- reshape::melt(m)


ggplot(data = heat, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()+
  geom_text(data=heat,aes(label=round(value, 2)))+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1))
