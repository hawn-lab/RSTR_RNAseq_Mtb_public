library(tidyverse)
library(ggrepel)
library(viridis)
library(cowplot)

#### Fold change ####
attach("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")
voomU <- dat.combined.voom

#Calculate mean Mtb-media FC
FC.U <- as.data.frame(voomU$E) %>% 
  #add metadata
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "libID") %>% 
  left_join(voomU$targets) %>% 
  #calculate TB/MEDIA fold change
  select(rowname, FULLIDNO, Sample_Group, condition, value) %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(FC = TB-MEDIA) %>% 
  #summarise mean FC within RSTR group
  group_by(rowname, Sample_Group) %>% 
  summarise(meanFC = mean(FC, na.rm=TRUE)) %>% 
  rename(group=Sample_Group) %>% 
  #Calculate RSTR-LTBI of mean Mtb-media FC
  pivot_wider(names_from = group, values_from = meanFC) %>% 
  rowwise() %>% 
  mutate(FC.rstr = RSTR-LTBI)

#Add FDR
FC.fdr.U <- read_csv("../Uganda/results/gene_level/RSTR.Mtb.model.results.anno.csv") %>% 
  filter(grepl(":",variable)) %>% 
  select(gene, FDR) %>% 
  inner_join(FC.U, by=c("gene"="rowname")) %>% 
  mutate(group="Uganda")

attach("../SAfrica/data_clean/RSTR_SA_dat_clean.RData")
voomSA <- dat.abund.norm.voom

#Calculate mean Mtb-media FC
FC.SA <- as.data.frame(voomSA$E) %>% 
  #add metadata
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "libID") %>% 
  left_join(voomSA$targets) %>% 
  #calculate TB/MEDIA fold change
  select(rowname, FULLIDNO, Sample_Group, condition, value) %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(FC = TB-MEDIA) %>% 
  #summarise mean FC within RSTR group
  group_by(rowname, Sample_Group) %>% 
  summarise(meanFC = mean(FC, na.rm=TRUE)) %>% 
  rename(group=Sample_Group) %>% 
  #Calculate RSTR-LTBI of mean Mtb-media FC
  pivot_wider(names_from = group, values_from = meanFC) %>% 
  rowwise() %>% 
  mutate(FC.rstr = RSTR-LTBI)

#Add FDR
FC.fdr.SA <- read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.model.results.anno.csv") %>% 
  filter(grepl(":",variable)) %>% 
  select(gene, FDR) %>% 
  inner_join(FC.U, by=c("gene"="rowname"))%>% 
  mutate(group="South Africa")
  
#### Combine U and SA ####
dat.all <- full_join(FC.fdr.U, FC.fdr.SA) %>% 
  #create color variable for significance
  mutate(col.group = ifelse(FDR > 0.2, "FDR > 0.2",
                            ifelse(FDR <= 0.2 & FDR > 0.05, "0.05 < FDR < 0.2",
                                   ifelse(FDR <= 0.05, "FDR < 0.05",
                                          NA)))) %>% 
  #create gene label group
  mutate(lab.group = ifelse(group == "Uganda" &
             gene %in% c("CIITA","CXCL9","IFNG","IL2","DDX10","UMPS"), 1, 
    ifelse(group == "South Africa" & 
             gene %in% c("ASNA1","GANAB","MRPL54","MTUS2","STAU2"),1,
           NA))) %>% 
  #dataset order
  mutate(group = factor(group, level=c("Uganda","South Africa")))

#### Plot ####
plot <- dat.all %>% 
  
  ggplot(aes(x=FC.rstr, y=-log10(FDR))) +
  #Interaction signif < 0.2
  geom_point(data=filter(dat.all, col.group == "0.05 < FDR < 0.2"),
             aes(color=-log10(FDR)), size=1) +  
  #Interaction not signif 
  geom_point(data=filter(dat.all, col.group == "FDR > 0.2"),
             aes(fill=col.group), shape=21, size=1, color="#d9d9d9") +
  #FDR < 0.05
  geom_point(data=filter(dat.all, col.group == "FDR < 0.05"),
             aes(fill=col.group), shape=21, size=1, color="red") +
  scale_fill_manual(values = c("red","#d9d9d9"),
                    guide = guide_legend(override.aes = list(color="white",
                                                             size=3) ) ) +

  #Label genes
  geom_text_repel(dat=filter(dat.all, !is.na(lab.group)),
                  aes(label=gene), min.segment.length = unit(0, 'lines')) +
  facet_wrap(~group,ncol=1) +

  geom_hline(yintercept = -log10(0.2), lty="dashed") +
  geom_vline(xintercept = 0, lty="dashed") +
  
  theme_classic() +
  labs(x="Mean log2 Mtb:RSTR fold change\n(Mtb - media in RSTR) / (Mtb - media in LTBI)", 
       y="-log10( Mtb:RSTR FDR )",
       color="-log10( FDR )", fill="") 
#plot

#### Save ####

ggsave("Fig2A.volcano.plot.pdf", plot_grid(plot, labels = "a)", label_fontface = "plain"),
       width=4.2, height=6)

tiff("Fig2A.volcano.plot.tiff", units="in", width=4.2, height=6, res=300)
plot_grid(plot, labels = "a)", label_fontface = "plain")
dev.off()
