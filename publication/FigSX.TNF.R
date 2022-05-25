library(tidyverse)
library(kimma)
library(ggpubr)

#### Expression data ####
genes.OI <- c("TNF")

attach("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")
voomU <- dat.combined.voom
voomU$targets <- voomU$targets %>% 
  mutate(Sample_Group = factor(Sample_Group, c("LTBI","RSTR")))

attach("../SAfrica/data_clean/RSTR_SA_dat_clean.RData")
voomSA <- dat.abund.norm.voom
voomSA$targets <- voomSA$targets %>% 
  mutate(Sample_Group = factor(Sample_Group, c("LTBI","RSTR")))

kin <- read_csv("data_raw/kinship_Hawn_all.csv") %>% 
  column_to_rownames() %>% as.matrix()

#### Model ####
resultU <- kmFit(voomU, kin, patientID="FULLIDNO",
                 model="~condition*Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + experiment + (1|FULLIDNO)",
                 run.lmekin = TRUE, 
                 run.contrast = TRUE, contrast.var = "condition:Sample_Group",
                 subset.genes = genes.OI)

resultSA <- kmFit(voomSA, patientID="FULLIDNO",
                 model="~condition*Sample_Group + age + (1|FULLIDNO)",
                 run.lme = TRUE, 
                 run.contrast = TRUE, contrast.var = "condition:Sample_Group",
                 subset.genes = genes.OI)

#### Plot data ####
datU <- voomU$E %>% 
  rownames_to_column("gene")%>% 
  filter(gene %in% c(genes.OI)) %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(select(voomU$targets, libID, Sample_Group)) %>% 
  mutate(dataset="Uganda")

datSA <- as.data.frame(voomSA$E) %>% 
  rownames_to_column("gene")%>% 
  filter(gene %in% c(genes.OI)) %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(select(voomSA$targets, libID, Sample_Group)) %>% 
  mutate(dataset="South Africa")

dat.all <- datU %>% 
  full_join(datSA) %>% 
  separate(libID, into=c("ID","condition"), sep="_", remove=FALSE) %>% 
  mutate(condition = recode(condition, "MEDIA"="media","TB"="+Mtb"),
         x = paste(condition, Sample_Group, sep="\n"),
         x = factor(x, levels=c("media\nLTBI","media\nRSTR",
                                "+Mtb\nLTBI","+Mtb\nRSTR"))) %>% 
  mutate(dataset = factor(dataset, levels=c("Uganda","South Africa")))

pval <- resultU$lmekin.contrast %>% 
  full_join(resultSA$lme.contrast) %>% 
  mutate(contrast = recode(contrast, 
                           "MEDIA LTBI - MEDIA RSTR"="MEDIA_RSTR - MEDIA_LTBI",
                           "TB LTBI - TB RSTR"="TB_RSTR - TB_LTBI")) %>% 
  filter(contrast %in% c("MEDIA_RSTR - MEDIA_LTBI","TB_RSTR - TB_LTBI")) %>% 
  separate(contrast, into=c("group1", "group2"), sep=" - ") %>% 
  mutate(across(c(group1, group2),
                ~gsub("MEDIA","media",gsub("^TB","+Mtb",gsub("_","\n",.))))) %>% 
  mutate(pval = paste("p =", signif(pval, digits=2))) %>% 
  arrange(group1,model) %>% 
  mutate(y.position = c(11.5,11.5,6,6),
         dataset=c("South Africa","Uganda","South Africa","Uganda"))%>% 
  mutate(dataset = factor(dataset, levels=c("Uganda","South Africa")))

#### Plot ####
plot1 <- dat.all %>% 
  filter(gene %in% genes.OI) %>% 
  
  ggplot(aes(x=x, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes( color=Sample_Group),
             position=position_jitterdodge(), size=1) +
  theme_classic() +
  facet_wrap(dataset~gene) +
  labs(y="Normalized log2 expression", x="", fill="") +
  theme(legend.position = "none") +
  stat_pvalue_manual(data=pval, 
                   label="pval", xmin="group1", xmax="group2") +
  lims(y=c(2.4,12))

plot1

#### Save ####
ggsave(plot1, 
       filename="FigS.TNF.boxplot.pdf", width=5, height=4)

# tiff("FigS5.gene.boxplot.tiff", width=7, height=6, units="in", res=300)
# plot_grid(plot1,row2, labels = c("a)",""), 
#           label_fontface = "plain", nrow=2)
# dev.off()