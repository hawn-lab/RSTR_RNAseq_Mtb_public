library(tidyverse)
library(cowplot)
set.seed(4389)

#### Uganda ####
load("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")

#Calculate PCA
PCA <- as.data.frame(dat.combined.voom$E) %>% 
  t() %>% 
  prcomp()

PC1.label <- paste("PC1 (", summary(PCA)$importance[2,1]*100, "%)", sep="")
PC2.label <-paste("PC2 (", summary(PCA)$importance[2,2]*100, "%)", sep="")

# Extract PC values
PCA.dat <- as.data.frame(PCA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.combined.voom$targets, by="libID")

PCA1 <- ggplot(PCA.dat, aes(PC1, PC2)) +
    geom_point(aes(color=condition),
               size=2) +
    #Beautify
    theme_classic() +
    labs(x=PC1.label, y=PC2.label, 
         title="Uganda PCA by Mtb stim") +
    coord_fixed(ratio=1) +
    guides(color=guide_legend(title.position="top", 
                              title.hjust = 0.5)) +
  scale_color_manual(name = "Stim", values=c("#7CAE00","#C77CFF"), 
                     labels = c("Media", "Mtb"))
PCA1

PCA2 <- ggplot(PCA.dat, aes(PC1, PC2)) +
  geom_point(aes(color=Sample_Group),
             size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label, y=PC2.label, 
       title="Uganda PCA by clinical phenotype") +
  coord_fixed(ratio=1) +
  guides(color=guide_legend(title.position="top", 
                            title.hjust = 0.5)) + 
  scale_color_discrete(name = "Phenotype")
PCA2

#### SAfrica ####
load("../SAfrica/data_clean/RSTR_SA_dat_clean.RData")

#Calculate PCA
PCA.SA <- as.data.frame(dat.abund.norm.voom$E) %>% 
  t() %>% 
  prcomp()

PC1.label2 <- paste("PC1 (", summary(PCA.SA)$importance[2,1]*100, "%)", sep="")
PC2.label2 <-paste("PC2 (", summary(PCA.SA)$importance[2,2]*100, "%)", sep="")

# Extract PC values
PCA.dat2 <- as.data.frame(PCA.SA$x) %>% 
  rownames_to_column("libID") %>%
  # Select PCs
  dplyr::select(libID, PC1:PC3) %>% 
  # Merge with metadata
  left_join(dat.abund.norm.voom$targets, by="libID")

PCA3 <- ggplot(PCA.dat2, aes(PC1, PC2)) +
  geom_point(aes(color=condition),
             size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label2, y=PC2.label2, 
       title="South Africa PCA by Mtb stim") +
  coord_fixed(ratio=1) +
  guides(color=guide_legend(title.position="top", 
                            title.hjust = 0.5)) +
  scale_color_manual(name = "Stim", values=c("#7CAE00","#C77CFF"),
                     labels = c("Media", "Mtb"))
PCA3

PCA4 <- ggplot(PCA.dat2, aes(PC1, PC2)) +
  geom_point(aes(color=Sample_Group),
             size=2) +
  #Beautify
  theme_classic() +
  labs(x=PC1.label2, y=PC2.label2, 
       title="South Africa PCA by clinical phenotype") +
  coord_fixed(ratio=1) +
  guides(color=guide_legend(title.position="top", 
                            title.hjust = 0.5)) +
  scale_color_discrete(name = "Phenotype")
  
PCA4

#### Save ####
ggsave("FigS2.PCA.pdf", 
       plot_grid(PCA1,PCA2,PCA3,PCA4, align = "hv", axis="lrtb",
                 ncol=2,labels = c("a)","b)","c)","d)"), label_fontface = "plain"),
       height=6, width=8)
tiff("FigS2.PCA.tiff", height=6, width=8, units="in", res=300)
plot_grid(PCA1,PCA2,PCA3,PCA4, align = "hv", axis="lrtb",
          ncol=2,labels = c("a)","b)","c)","d)"), label_fontface = "plain")
dev.off()