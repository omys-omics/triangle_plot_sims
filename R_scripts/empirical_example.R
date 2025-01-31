setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/Triangle_Plots_best_practices/empirical_example")


# Load packages
library(triangulaR)
library(vcfR)
library(SNPfiltR)
library(bgchm)
library(ggpubr)



# Read in empirical data
fox_sparrow_90 <- read.vcfR("snps_3_90.vcf.gz")
fox_sparrow_100 <- read.vcfR("snps_3_100.vcf.gz")
fox_sparrow_pm <- read.table("passerella_samples_3.txt", header = T)


############################################################
##   Uneven parental sampling (18 & 5) vs even (5 & 5)    ##
###########################################################
d_18_5 <- specFreqDiff(fox_sparrow_100, fox_sparrow_pm, p1 = "ili", p2 = "una")
mean(d_18_5$diff)
sum(d_18_5$diff==1)

d_18_5_spec_diffs <- ggplot(data = d_18_5,  aes(x=diff)) +
  geom_histogram(color="#FFFFFF", alpha=1, position = 'identity', bins = 15) +
  xlab(paste("Observed Allele Frequency Difference")) +
  ylab(paste("Count")) +
  ylim(0,8000) +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

ggsave("d_18_5_spec_diffs.pdf", plot = d_18_5_spec_diffs, path = "./", width = 4, height = 6)


aims_tri_18_5 <- alleleFreqDiff(vcfR = fox_sparrow_100, pm = fox_sparrow_pm, 
                                   p1 = "ili", p2 = "una", difference = 1)
hi_het_tri_18_5 <-hybridIndex(vcfR = aims_tri_18_5, pm = fox_sparrow_pm, p1 = "ili", p2 = "una")
triangle.plot(hi_het_tri_18_5)


# randomly choose 5 ili parentals
fox_sparrow_pm_5 <- fox_sparrow_pm
ili_pars <- fox_sparrow_pm_5[fox_sparrow_pm$pop=="ili",]$id
set.seed(894612)
fox_sparrow_pm_5[fox_sparrow_pm$id %in% sample(ili_pars, 5, replace = F),"pop"] <- "ili_ref"



aims_tri_5_5 <- alleleFreqDiff(vcfR = fox_sparrow_100, pm = fox_sparrow_pm_5, 
                                   p1 = "ili_ref", p2 = "una", difference = 1)


hi_het_tri_5_5 <-hybridIndex(vcfR = aims_tri_5_5, pm = fox_sparrow_pm_5, p1 = "ili_ref", p2 = "una")
triangle.plot(hi_het_tri_5_5)



# convert to genotype matrix
aims_bgchm_18_5 <- t(extract.gt(aims_tri_18_5))
identical(row.names(aims_bgchm_18_5), fox_sparrow_pm$id)
aims_bgchm_18_5[aims_bgchm_18_5=="0/0"] <- 0
aims_bgchm_18_5[aims_bgchm_18_5=="0/1"] <- 1
aims_bgchm_18_5[aims_bgchm_18_5=="1/0"] <- 1
aims_bgchm_18_5[aims_bgchm_18_5=="1/1"] <- 2

#convert to numeric and add row names
aims_bgchm_18_5 <- apply(aims_bgchm_18_5, 2, as.numeric)
row.names(aims_bgchm_18_5) <- fox_sparrow_pm$id

# make matrices for parental individuals
P1_bgchm_18_5 <- aims_bgchm_18_5[fox_sparrow_pm[fox_sparrow_pm$pop=="ili",]$id,]
P2_bgchm_18_5 <- aims_bgchm_18_5[fox_sparrow_pm[fox_sparrow_pm$pop=="una",]$id,]

# run bgchm
date()
q_out_18_5<-est_Q(Gx=aims_bgchm_18_5,G0=P1_bgchm_18_5,G1=P2_bgchm_18_5,model="genotype",ploidy="diploid")
date()
tri_plot(hi=q_out_18_5$hi[,1],Q10=q_out_18_5$Q10[,1],pdf=FALSE,pch=19)
segments(q_out_18_5$hi[,1],q_out_18_5$Q10[,2],q_out_18_5$hi[,1],q_out_18_5$Q10[,3])
segments(q_out_18_5$hi[,2],q_out_18_5$Q10[,1],q_out_18_5$hi[,3],q_out_18_5$Q10[,1])

# make hi het 
hi_het_bgchm_18_5 <- data.frame(bgchm_hi_18_5 = q_out_18_5$hi[,1],
                                bgchm_het_18_5 = q_out_18_5$Q10[,1])


# convert to genotype matrix
aims_bgchm_5_5 <- t(extract.gt(aims_tri_5_5))
identical(row.names(aims_bgchm_5_5), fox_sparrow_pm_5$id)
aims_bgchm_5_5[aims_bgchm_5_5=="0/0"] <- 0
aims_bgchm_5_5[aims_bgchm_5_5=="0/1"] <- 1
aims_bgchm_5_5[aims_bgchm_5_5=="1/0"] <- 1
aims_bgchm_5_5[aims_bgchm_5_5=="1/1"] <- 2

#convert to numeric and add row names
aims_bgchm_5_5 <- apply(aims_bgchm_5_5, 2, as.numeric)
row.names(aims_bgchm_5_5) <- fox_sparrow_pm_5$id

# make matrices for parental individuals
P1_bgchm_5_5 <- aims_bgchm_5_5[fox_sparrow_pm_5[fox_sparrow_pm_5$pop=="ili_ref",]$id,]
P2_bgchm_5_5 <- aims_bgchm_5_5[fox_sparrow_pm_5[fox_sparrow_pm_5$pop=="una",]$id,]

# run bgchm
date()
q_out_5_5<-est_Q(Gx=aims_bgchm_5_5,G0=P1_bgchm_5_5,G1=P2_bgchm_5_5,model="genotype",ploidy="diploid")
date()
tri_plot(hi=q_out_5_5$hi[,1],Q10=q_out_5_5$Q10[,1],pdf=FALSE,pch=19)
segments(q_out_5_5$hi[,1],q_out_5_5$Q10[,2],q_out_5_5$hi[,1],q_out_5_5$Q10[,3])
segments(q_out_5_5$hi[,2],q_out_5_5$Q10[,1],q_out_5_5$hi[,3],q_out_5_5$Q10[,1])

# make hi het 
hi_het_bgchm_5_5 <- data.frame(bgchm_hi_5_5 = q_out_5_5$hi[,1],
                                bgchm_het_5_5 = q_out_5_5$Q10[,1])

hi_het_all <- data.frame(id = hi_het_tri_18_5$id,
                         pop = hi_het_tri_18_5$pop,
                         tri_18_5_hi = hi_het_tri_18_5$hybrid.index,
                         tri_18_5_het = hi_het_tri_18_5$heterozygosity,
                         tri_5_5_hi = hi_het_tri_5_5$hybrid.index,
                         tri_5_5_het = hi_het_tri_5_5$heterozygosity,
                         bgchm_18_5_hi = hi_het_bgchm_18_5$bgchm_hi_18_5,
                         bgchm_18_5_het = hi_het_bgchm_18_5$bgchm_het_18_5,
                         bgchm_5_5_hi = hi_het_bgchm_5_5$bgchm_hi_5_5,
                         bgchm_5_5_het = hi_het_bgchm_5_5$bgchm_het_5_5)

compare.hi <- ggplot(hi_het_all, ) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
  geom_point(aes(x=bgchm_5_5_hi, y=tri_5_5_hi, colour = "5_5"), cex = 1.5) +
  geom_point(aes(x=bgchm_18_5_hi, y=tri_18_5_hi, colour = "18_5"), cex = 1.5) +
  ggtitle("Hybrid Index Estimates") +
  xlab(paste("bgchm")) +
  ylab(paste("triangulaR")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  ylim(c(0, 1.05)) + 
  xlim(c(0, 1.05)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

compare.het <- ggplot(hi_het_all, ) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
  geom_point(aes(x=bgchm_5_5_het, y=tri_5_5_het, colour = "5_5"), cex = 1.5) +
  geom_point(aes(x=bgchm_18_5_het, y=tri_18_5_het, colour = "18_5"), cex = 1.5) +
  ggtitle("Interclass Heterozygosity Estimates") +
  xlab(paste("bgchm")) +
  ylab(paste("triangulaR")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  ylim(c(0, 1.05)) + 
  xlim(c(0, 1.05)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


set.seed(1231423)
tri_18_5 <- triangle.plot(hi_het_tri_18_5, jitter = 0, colors = c("black", "#785EF0", "#FFB000"), cex = 2) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10)) 
  

set.seed(1231423)
tri_5_5 <- triangle.plot(hi_het_tri_5_5, jitter = 0, colors = c("black", "#648FFF", "#785EF0", "#FFB000"), cex = 2) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10)) 




####################################################
#    Genetic and morphological hybrid indices      #
####################################################
passerella_plumage <- read.table("passerella_plumage_scores.csv", sep = ",", header = T)

gen_morph <- hi_het_tri_18_5[hi_het_tri_18_5$pop=="con",]

# make id just the catalog number
gen_morph$id <- sapply(strsplit(gen_morph$id, split = "_"), '[', 3)
gen_morph$id <- sapply(strsplit(gen_morph$id, split = "b"), '[', 1)
gen_morph$id <- sapply(strsplit(gen_morph$id, split = "a"), '[', 1)

# make sure all gen ids are in the passerella plumage dataset
all(gen_morph$id %in% passerella_plumage$catalog_no)

# add plumage score
gen_morph$scaled_plumage_score <- NA
for (i in gen_morph$id) {
  if (!i %in% passerella_plumage$catalog_no) {
    next
  }
  gen_morph[gen_morph$id==i, "scaled_plumage_score"] <- passerella_plumage[passerella_plumage$catalog_no==i, "scaled_plumage_score"]
}


gen_morph_plot <- ggplot(gen_morph, ) +
  #geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black") +
  geom_smooth(aes(x=scaled_plumage_score, y=hybrid.index), method = "lm", color = "black") +
  geom_point(aes(x=scaled_plumage_score, y=hybrid.index), cex = 2) +
  annotate(geom="text", x=1, y=0.08, hjust = 1, label=expression(R^2 == 0.97), size=4) +
  annotate(geom="text", x=1, y=0, hjust = 1, label=expression(p < 2.2e-16), size=4) +
  ggtitle("") +
  xlab(paste("Plumage Score")) +
  ylab(paste("Hybrid Index")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

summary(lm(hybrid.index ~ scaled_plumage_score, gen_morph))




sparrows <- ggarrange(gen_morph_plot, tri_18_5, compare.hi, 
                      gen_morph_plot, tri_5_5, compare.het, 
                      nrow = 2, ncol = 3,
                      widths = c(1,1.2,1))

ggsave("sparrows_fromR.pdf", plot = sparrows, path = "./", width = 9, height = 5.8)




#################################
##   d and missing sampling    ##
#################################

for (i in c("fox_sparrow_100", "fox_sparrow_90")) {
  for (j in c(1,0.75,0.5)) {
    fox_sparrow_aims <- alleleFreqDiff(vcfR = get(i), pm = fox_sparrow_pm_5, 
                                       p1 = "ili_ref", p2 = "una", difference = j)
    d <- nrow(extract.gt(fox_sparrow_aims))
    fox_sparrow_hi_het <-hybridIndex(vcfR = fox_sparrow_aims, pm = fox_sparrow_pm_5, p1 = "ili_ref", p2 = "una")
    p <- triangle.plot(fox_sparrow_hi_het, colors = c("black", "#648FFF", "#785EF0", "#FFB000"))  +
      annotate(geom="text", x=-0.05, y=0.95, hjust = 0, label=paste0(d, " SNPs"), size=4) +
      scale_y_continuous(breaks = c(0,0.5,1)) +
      scale_x_continuous(breaks = c(0,0.5,1)) +
      theme(legend.position = "none", axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
    assign(paste0(i,"_",j,"_triplot"), p)
  }
}

library(ggpubr)
tris <- ggarrange(fox_sparrow_100_1_triplot, fox_sparrow_90_1_triplot, 
          fox_sparrow_100_0.75_triplot, fox_sparrow_90_0.75_triplot, 
          fox_sparrow_100_0.5_triplot, fox_sparrow_90_0.5_triplot,
          #labels = c("A", "B", "C", "D", "E", "F"),
          #font.label = list(size = 14),
          #label.y = -0.2,
          #label.x = 0.1,
          ncol = 2, nrow = 3)


ggsave("fox_sparrow_tris_fromR.pdf", plot = tris, path = "./", width = 4.6666667, height = 6)





