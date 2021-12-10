library(tidyverse)
library(cowplot)
library(kableExtra)

#### Expression data ####
genes.OI <- c("IFNG","IL2","CXCL9")
genes.OI2 <- c("CIITA","BTN2A2")

attach("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")
voomU <- dat.combined.voom

datU <- voomU$E %>% 
  rownames_to_column("gene")%>% 
  filter(gene %in% c(genes.OI,genes.OI2)) %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(select(voomU$targets, libID, Sample_Group)) %>% 
  mutate(dataset="Uganda")

dat.all <- datU %>% 
  separate(libID, into=c("ID","condition"), sep="_", remove=FALSE) %>% 
  mutate(condition = recode(condition, "MEDIA"="media","TB"="+Mtb"),
         x = paste(condition, Sample_Group, sep="\n"),
         x = factor(x, levels=c("media\nLTBI","media\nRSTR",
                                "+Mtb\nLTBI","+Mtb\nRSTR"))) %>% 
  mutate(dataset = factor(dataset, levels=c("Uganda","South Africa")))

#### Plot ####
plot1 <- dat.all %>% 
  filter(gene %in% genes.OI) %>% 
  
  ggplot(aes(x=x, y=value,)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes( color=Sample_Group),
             position=position_jitterdodge(), size=1) +
  theme_classic() +
  facet_wrap(~gene, scales="free") +
  labs(y="Normalized log2 expression", x="", fill="") +
  theme(legend.position = "none")
plot1

plot2 <- dat.all %>% 
  filter(gene %in% genes.OI2) %>% 
  
  ggplot(aes(x=x, y=value,)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes( color=Sample_Group),
             position=position_jitterdodge(), size=1) +
  theme_classic() +
  facet_wrap(~gene, scales="free") +
  labs(y="Normalized log2 expression", x="", fill="") +
  theme(legend.position = "none")
plot2

#### Save ####
row2 <- plot_grid(plot2, NULL, labels = c("b)",""), 
                  label_fontface = "plain",
                  rel_widths = c(2,1))
ggsave(plot_grid(plot1,row2, labels = c("a)",""), 
                 label_fontface = "plain", nrow=2), 
       filename="FigS5.gene.boxplot.pdf", width=7, height=6)

tiff("FigS5.gene.boxplot.tiff", width=7, height=6, units="in", res=300)
plot_grid(plot1,row2, labels = c("a)",""), 
          label_fontface = "plain", nrow=2)
dev.off()
