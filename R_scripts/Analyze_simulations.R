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


########################################
###         DEFINE FUNCTIONS       #####
########################################


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
    p <- round((sum(sfs < min.mac)/length(sfs)) * 100, 2)
    message(p, "% of SNPs fell below a minor allele count of ", 
            min.mac, " and were removed from the VCF")
    vcfR <- vcfR[sfs >= min.mac, ]
    return(vcfR)
  }
}



makePopmap <- function(file = NULL, vcf = NULL) {
  vcfR <- read.vcfR(vcf)
  num <- read_table(file)
  pm <- data.frame(matrix(nrow=(sum(num)), ncol=2))
  colnames(pm) <- c("id", "pop")
  pm$id <- colnames(vcfR@gt)[-1]
  pm$pop <- c(rep("p0", num[,1]), rep("p20", num[,2]))
  return(pm)
}

known.freqs <- function(file = NULL, fixed.diffs = NULL) {
  vcfR <- read.vcfR(file)
  
  # Subset vcfR file
  indices <- match(fixed.diffs, rownames(extract.gt(vcfR)))
  vcfR <- vcfR[indices,]
  
  # Make popmap and assign populations
  sim.popmap <- data.frame(matrix(nrow = length(colnames(extract.gt(vcfR)))))
  colnames(sim.popmap) <- "id"
  sim.popmap$id <- colnames(extract.gt(vcfR))
  num.pops <- nrow(sim.popmap) / 20
  sim.popmap$pop <- sort(c(rep(0:(num.pops-1),20)))
  
  # Calculate metrics at known fixed differences
  hi.het <- hybridIndex(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20")
  return(hi.het)
}


makePlots <- function(s = NULL, gen = NULL, difference = NULL) {
  v <- paste0(s, "/gen.", gen, ".vcf.gz")
  vcfR <- read.vcfR(v)
  vcfR <- min.mac(vcfR, min.mac = 1)
  vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = difference)
  num <- nrow(vcfR.diff@gt)
  hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
  tri <- triangulaR::triangle.plot(data = hi.het, colors = colors, cex = 3, alpha = 0.5, jitter = 0.01) + 
          labs(title = paste0(s, ", gen: ", gen)) + 
          annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0(num, " SNPs"), size=4) +
          theme(legend.position = "none")
  box <- ggplot(data = hi.het,  aes(x=pop, y=hybrid.index, fill=as.factor(pop), group=pop)) +
          geom_boxplot(alpha=1) +
          scale_fill_manual(values=colors) +
          xlab(paste("Population")) +
          ylab(paste("Hybrid Index")) +
          theme_classic() +
          scale_y_continuous(limits = c(0,1)) +
          theme(legend.position = "none") +
          labs(title = paste0(s, ", gen: ", gen))
  return(list(tri, box))
}





#######################################
#            Name of AIMs             #
#######################################

D750.beg.of.con.pm <- makePopmap(file = "D750/D750.num.parentals.beg.of.con.txt", vcf = "D750/D750.beg.of.con.vcf")
D750.beg.of.con.vcfR <- read.vcfR("D750/D750.beg.of.con.vcf")
D750.beg.of.con.vcfR <- min.mac(D750.beg.of.con.vcfR, min.mac = 1)

D1000.beg.of.con.pm <- makePopmap(file = "D1000/D1000.num.parentals.beg.of.con.txt", vcf = "D1000/D1000.beg.of.con.vcf")
D1000.beg.of.con.vcfR <- read.vcfR("D1000/D1000.beg.of.con.vcf")
D1000.beg.of.con.vcfR <- min.mac(D1000.beg.of.con.vcfR, min.mac = 1)

D2000.beg.of.con.pm <- makePopmap(file = "D2000/D2000.num.parentals.beg.of.con.txt", vcf = "D2000/D2000.beg.of.con.vcf")
D2000.beg.of.con.vcfR <- read.vcfR("D2000/D2000.beg.of.con.vcf")
D2000.beg.of.con.vcfR <- min.mac(D2000.beg.of.con.vcfR, min.mac = 1)

D750.fixed.diffs <- AIMnames(vcfR = D750.beg.of.con.vcfR, pm = D750.beg.of.con.pm, p1 = "p0", p2 = "p20", difference = 1)
write.table(D750.fixed.diffs, file = "D750/D750.fixed.diffs.txt", quote = F, sep = "\t", row.names = F, col.names = F)
D1000.fixed.diffs <- AIMnames(vcfR = D1000.beg.of.con.vcfR, pm = D1000.beg.of.con.pm, p1 = "p0", p2 = "p20", difference = 1)
write.table(D1000.fixed.diffs, file = "D1000/D1000.fixed.diffs.txt", quote = F, sep = "\t", row.names = F, col.names = F)
D2000.fixed.diffs <- AIMnames(vcfR = D2000.beg.of.con.vcfR, pm = D2000.beg.of.con.pm, p1 = "p0", p2 = "p20", difference = 1)
write.table(D2000.fixed.diffs, file = "D2000/D2000.fixed.diffs.txt", quote = F, sep = "\t", row.names = F, col.names = F)


# Write out popmaps
write.table(D750.beg.of.con.pm, file = "D750/D750.beg.of.con.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(D1000.beg.of.con.pm, file = "D1000/D1000.beg.of.con.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(D2000.beg.of.con.pm, file = "D2000/D2000.beg.of.con.pm.txt", quote = F, sep = "\t", row.names = F, col.names = F)




########################################
###      Known fixed differences     ###
########################################


color_ramp <- colorRampPalette(c("#313695", "yellow3", "#a50026"))
colors <- color_ramp(21)

for (f in paste0("D2000/gen.", seq(6000,6000,by=200), ".vcf.gz")) {
  hi.het <- known.freqs(file = f, fixed.diffs = D2000.fixed.diffs)
  print(triangulaR::triangle.plot(data = hi.het, colors = colors, cex = 3, alpha = 0.5, jitter = 0.01) + 
    annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("Gen: ", f), size=4))
  
  print(ggplot(data = hi.het,  aes(x=pop, y=hybrid.index, fill=as.factor(pop), group=pop)) +
    geom_boxplot(alpha=1) +
    scale_fill_manual(values=colors) +
    xlab(paste("Population")) +
    ylab(paste("Hybrid Index")) +
    theme_classic() +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = "none")
  )
}


# clines using known differences
for (s in c("D750", "D1000", "D2000")) {
  for (i in c(0,200,1000,2000,3000,4000,5000,6000)) {
    f <- paste0(s, "/gen.", i, ".vcf.gz")
    if (i == 0) {
      hi.het <- known.freqs(file = f, fixed.diffs = get(paste0(s, ".fixed.diffs")))
      hi.het$gen <- i
    } else {
      temp <- known.freqs(file = f, fixed.diffs = get(paste0(s, ".fixed.diffs")))
      temp$gen <- i
      hi.het <- rbind(hi.het, temp)
    }
  }
  curve <- ggplot(data = hi.het, aes(x=pop, y=hybrid.index, color=gen, group=gen)) +
    stat_summary(geom="line",fun="mean",size=1) +
    #scale_fill_manual(values=colors) +
    xlab(paste("Population")) +
    ylab(paste("Hybrid Index")) +
    theme_classic() +
    scale_y_continuous(limits = c(0,1)) +  
    #labs(title = s) +
    theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

  
  assign(paste0(s, ".known.curve"), curve)
}


known.curves <- ggarrange(D750.known.curve, D1000.known.curve, D2000.known.curve, 
                   labels = c("A", "B", "C"),
                   font.label = list(size = 14),
                   #label.y = -0.2,
                   label.x = 0.1,
                   ncol = 3, nrow = 1)










############################################
###         Find fixed differences       ###
############################################
# Make popmap and assign populations
sim.popmap <- data.frame(matrix(nrow = 420))
colnames(sim.popmap) <- "id"
sim.popmap$id <- paste0("i", 0:419)
num.pops <- nrow(sim.popmap) / 20
sim.popmap$pop <- sort(c(rep(0:(num.pops-1),20)))
write.table(sim.popmap, file = "sim.popmap.txt", quote = F, sep = "\t", row.names = F, col.names = F)



# Clines using fixed differences
for (s in c("D750", "D1000", "D2000")) {
  for (i in c(0,200,1000,2000,3000,4000,5000,6000)) {
    v <- paste0(s, "/gen.", i, ".vcf.gz")
    if (i == 0) {
      vcfR <- read.vcfR(v)
      vcfR <- min.mac(vcfR, min.mac = 1)
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)
      hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
      hi.het$gen <- i
    } else {
      vcfR <- read.vcfR(v)
      vcfR <- min.mac(vcfR, min.mac = 1)
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)
      temp <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
      temp$gen <- i
      hi.het <- rbind(hi.het, temp)
    }
    print(i)
  }
  curve <- ggplot(data = hi.het, aes(x=pop, y=hybrid.index, color=gen, group=gen)) +
    stat_summary(geom="line",fun="mean",size=1) +
    #scale_fill_manual(values=colors) +
    xlab(paste("Population")) +
    ylab(paste("Hybrid Index")) +
    theme_classic() +
    scale_y_continuous(limits = c(0,1)) +
    #labs(title = s) +
    theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  
  
  assign(paste0(s, ".find.curve.1"), curve)
}


find.curve.1 <- ggarrange(D750.find.curve.1, D1000.find.curve.1, D2000.find.curve.1, 
                          labels = c("A", "B", "C"),
                          font.label = list(size = 14),
                          #label.y = -0.2,
                          label.x = 0.1,
                          ncol = 3, nrow = 1)



# Clines using threshold of 0.75
for (s in c("D750", "D1000", "D2000")) {
  for (i in c(0,200,1000,2000,3000,4000,5000,6000)) {
    v <- paste0(s, "/gen.", i, ".vcf.gz")
    if (i == 0) {
      vcfR <- read.vcfR(v)
      vcfR <- min.mac(vcfR, min.mac = 1)
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)
      hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
      hi.het$gen <- i
    } else {
      vcfR <- read.vcfR(v)
      vcfR <- min.mac(vcfR, min.mac = 1)
      vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)
      temp <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
      temp$gen <- i
      hi.het <- rbind(hi.het, temp)
    }
    print(i)
  }
  curve <- ggplot(data = hi.het, aes(x=pop, y=hybrid.index, color=gen, group=gen)) +
    stat_summary(geom="line",fun="mean",size=1) +
    #scale_fill_manual(values=colors) +
    xlab(paste("Population")) +
    ylab(paste("Hybrid Index")) +
    theme_classic() +
    scale_y_continuous(limits = c(0,1)) +
    #labs(title = s) +
    theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  
  
  assign(paste0(s, ".find.curve.0.75"), curve)
}


find.curve.0.75 <- ggarrange(D750.find.curve.0.75, D1000.find.curve.0.75, D2000.find.curve.0.75, 
                          labels = c("A", "B", "C"),
                          font.label = list(size = 14),
                          #label.y = -0.2,
                          label.x = 0.2,
                          ncol = 3, nrow = 1)



all.curves <- ggarrange(D750.known.curve, D1000.known.curve, D2000.known.curve, 
                        D750.find.curve.1, D1000.find.curve.1, D2000.find.curve.1, 
                        D750.find.curve.0.75, D1000.find.curve.0.75, D2000.find.curve.0.75, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          font.label = list(size = 14),
          #label.y = -0.2,
          label.x = 0.25,
          ncol = 3, nrow = 3)

ggsave("Figure6fromR.pdf", plot = all.curves, path = "Figures/", width = 8, height = 6)



#   GET LEGEND
fig6scale <- ggplot(data = hi.het, aes(x=pop, y=hybrid.index, color=gen, group=gen)) +
  stat_summary(geom="line",fun="mean",size=1) +
  xlab(paste("Population")) +
  ylab(paste("Hybrid Index")) +
  scale_fill_manual() +
  scale_colour_continuous(breaks=c(0,200,1000,2000,3000,4000,5000,6000)) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "right") +
  theme(legend.key.height= unit(2, 'cm'), legend.key.width= unit(0.5, 'cm'))

ggsave("Figure6scale_fromR.pdf", plot = fig6scale, path = "Figures/", width = 8, height = 8)








# Compare D750 gen 1000 to D2000 gen 5000
D750.gen1000.1 <- makePlots(s = "D750", gen = "1000", difference = 1)
D750.gen1000.0.75 <- makePlots(s = "D750", gen = "1000", difference = 0.75)
D750.gen1000.0.5 <- makePlots(s = "D750", gen = "1000", difference = 0.5)
D2000.gen5000.1 <- makePlots(s = "D2000", gen = "5000", difference = 1)
D2000.gen5000.0.75 <- makePlots(s = "D2000", gen = "5000", difference = 0.75)
D2000.gen5000.0.5 <- makePlots(s = "D2000", gen = "5000", difference = 0.5)


figsss <- ggarrange(D750.gen1000.1[[1]], D750.gen1000.1[[2]], D2000.gen5000.1[[1]], D2000.gen5000.1[[2]],
          D750.gen1000.0.75[[1]], D750.gen1000.0.75[[2]], D2000.gen5000.0.75[[1]], D2000.gen5000.0.75[[2]],
          D750.gen1000.0.5[[1]], D750.gen1000.0.5[[2]], D2000.gen5000.0.5[[1]], D2000.gen5000.0.5[[2]],
          #labels = c("A", "B", "C"),
          font.label = list(size = 14),
          #label.y = -0.2,
          label.x = 0.1,
          ncol = 4, nrow = 3)

ggsave("FigureSSS_fromR.pdf", plot = figsss, path = "Figures/", width = 12, height = 8)




#####################     MAKE FIGURE WITH ALL TRIANGLE PLOTS HERE

for (i in c(1000)) {
  v <- paste0("D750/gen.", i, ".vcf.gz")
  vcfR <- read.vcfR(v)
  #fds <- AIMnames(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)
  
  vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)
  
  hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
  
  print(triangulaR::triangle.plot(data = hi.het, colors = colors, cex = 3, alpha = 0.5, jitter = 0.01) + 
          annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("Gen: ", v), size=4))

  )
}













##############################################
#  Spectrum of allele frequency differences  #
##############################################
for (i in c(0,200,1000,2000,3000,4000,5000,6000,6600)) {
  v <- paste0("D2000/gen.", i, ".vcf.gz")
  vcfR <- read.vcfR(v)
  af <- specFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0)
  hist(af$diff)
}









for (i in c(1000)) {
  v <- paste0("D750/gen.", i, ".vcf.gz")
  vcfR <- read.vcfR(v)
  #fds <- AIMnames(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)
  
  vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)
  
  hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
  
  print(triangulaR::triangle.plot(data = hi.het, colors = colors, cex = 3, alpha = 0.5, jitter = 0.01) + 
          annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("Gen: ", v), size=4))
  print(ggplot(data = hi.het,  aes(x=pop, y=hybrid.index, fill=as.factor(pop), group=pop)) +
          geom_boxplot(alpha=1) +
          scale_fill_manual(values=colors) +
          xlab(paste("Population")) +
          ylab(paste("Hybrid Index")) +
          theme_classic() +
          scale_y_continuous(limits = c(0,1)) +
          theme(legend.position = "none")
  )
}




for (i in c(6000)) {
  v <- paste0("D2000/gen.", i, ".vcf.gz")
  vcfR <- read.vcfR(v)
  #fds <- AIMnames(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)
  
  vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)
  
  hi.het <- hybridIndex(vcfR = vcfR.diff, pm = sim.popmap, p1 = "0", p2 = "20")
  
  print(triangulaR::triangle.plot(data = hi.het, colors = colors, cex = 3, alpha = 0.5, jitter = 0.01) + 
          annotate(geom="text", x=-0.05, y=1.05, hjust =0, label=paste0("Gen: ", v), size=4))
  print(ggplot(data = hi.het,  aes(x=pop, y=hybrid.index, fill=as.factor(pop), group=pop)) +
          geom_boxplot(alpha=1) +
          scale_fill_manual(values=colors) +
          xlab(paste("Population")) +
          ylab(paste("Hybrid Index")) +
          theme_classic() +
          scale_y_continuous(limits = c(0,1)) +
          theme(legend.position = "none")
  )
}


















AIMnames(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 1)

afd <- aimFreqDist(vcfR = vcfR, pm = sim.popmap, p1 = "0", p2 = "20", difference = 0.75)



afd <- afd[,c("p1", "p2")]
afd <- melt(afd, variable_name = "pop")
afd[afd$value==0.5 & afd$pop=="p1","value"] <- 0.4999
afd[afd$value==0.5 & afd$pop=="p2","value"] <- 0.5001
ggplot(data = afd,  aes(x=value, fill=pop)) +
  geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 50) +
  scale_fill_manual(values=c("#af8dc3", "#7fbf7b")) +
  xlab(paste("Frequency of Parental Allele")) +
  ylab(paste("Count")) +
  theme_classic() +
  scale_y_continuous(limits = c(0,20)) +
  theme(legend.position = "none")
