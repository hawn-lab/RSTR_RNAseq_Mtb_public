library(tidyverse)
fdr.cutoff <- 0.1

#### Data ####
hU <- read_csv("../Uganda/results/GSEA/h_GSEA.result.csv")
hSA <- read_csv("../SAfrica/results/GSEA/h_GSEA.result.csv")
  
#### Significant terms ####
# Signif for at least 1 TB and 1 RSTR contrast
## RSTR
hU.signifR <- hU %>% 
  filter(group %in% c("RSTRinMEDIA","RSTRinTB") &
           fgsea.FDR <= fdr.cutoff) %>% 
  distinct(pathway) %>% unlist(use.names = FALSE)
hSA.signifR <- hSA %>% 
  filter(group %in% c("RSTRinMEDIA","RSTRinTB") &
           fgsea.FDR <= fdr.cutoff) %>% 
  distinct(pathway) %>% unlist(use.names = FALSE)

## TB
hU.signifTB <- hU %>% 
  filter(group %in% c("TBinLTBI","TBinRSTR") &
           fgsea.FDR <= fdr.cutoff) %>% 
  distinct(pathway) %>% unlist(use.names = FALSE)
hSA.signifTB <- hSA %>% 
  filter(group %in% c("TBinLTBI","TBinRSTR") &
           fgsea.FDR <= fdr.cutoff) %>% 
  distinct(pathway) %>% unlist(use.names = FALSE)

# Signifi in both datasets
signif.all <- intersect(intersect(hU.signifR, hU.signifTB),
                        intersect(hSA.signifR, hSA.signifTB))

# Combine and filter data
h.signif <- mutate(hU, dataset = "Uganda") %>% 
  bind_rows(mutate(hSA, dataset = "South Africa")) %>% 
  filter(pathway %in% signif.all)

## Save list for use in fold change plot
save(h.signif, file="GSEA.signif.RData")

#### Format data ####
h.signif.format <- h.signif %>% 
  #dataset order
  mutate(dataset = factor(dataset, level=c("Uganda","South Africa"))) %>% 
  #color by significance
  mutate(Significance = ifelse(fgsea.FDR <= fdr.cutoff, 
                               paste("FDR <", fdr.cutoff),
                               "NS")) %>% 
  #Beautify labels
  mutate(pathway = gsub("HALLMARK_","",pathway)) %>% 
  ##lowercase and then correction
  mutate(pathway = gsub("_"," ", tolower(pathway)),
         pathway = gsub("tnfa","TNFA",pathway),
         pathway = gsub("nfkb","NF-kB",pathway),
         pathway = gsub("interferon","IFN",pathway),
         pathway = gsub("il2","IL2",pathway),
         pathway = gsub("stat5","STAT5",pathway),
         pathway = gsub("kras","KRAS",pathway),
         pathway = gsub(" dn"," down",pathway)) %>% 
  ## order
  mutate(pathway = factor(pathway, 
                          levels=c(#Split to up R
                                   "coagulation", 
                                   
                                   #Split R
                                   "epithelial mesenchymal transition",
                                   "hypoxia",
                                   "KRAS signaling down",
                                   "IL2 STAT5 signaling",
                                   #All R down
                                   "IFN alpha response",
                                   "allograft rejection",
                                   "IFN gamma response",
                                   
                                   #Down with TB
                                   "adipogenesis",
                                   "oxidative phosphorylation",
                                   #Down to up R
                                   "inflammatory response",
                                   "TNFA signaling via NF-kB"))) %>% 
  mutate(group2 = ifelse(pathway %in% c("adipogenesis",
                                        "oxidative phosphorylation"), "ii",
                         ifelse(pathway %in% c("IFN alpha response",
                                               "allograft rejection",
                                               "IFN gamma response"), "iii",
                                ifelse(pathway %in% c("inflammatory response",
                                                      "TNFA signaling via NF-kB"), "i", "iv"))))

#### Plot param ####
#facet labels
facet_lab <- c('RSTRinMEDIA'="Down in RSTR <-  -> Up in RSTR",
               'RSTRinTB'="Down in RSTR <-  -> Up in RSTR",
               'TBinLTBI'="Down in +Mtb <-  -> Up in +Mtb",
               'TBinRSTR'="Down in +Mtb <-  -> Up in +Mtb",
               'i'='i','ii'='ii','iii'='iii','iv'='iv')

#Enrichment score limits
plot.lim <- max(abs(h.signif$fgsea.NES))+0.15

#### Plot ####
plot1 <- h.signif.format %>% 
  filter(group %in% c("RSTRinMEDIA","RSTRinTB")) %>% 
  ggplot() +
  geom_segment(aes(x=pathway, xend=pathway, 
                   y=0, yend=fgsea.NES),
               size=0.5) +
  geom_point(aes(x=pathway, y=fgsea.NES,
                 fill = Significance, size = dataset, shape=dataset),
             stroke=0.7) +
  geom_hline(yintercept = 0) +
    
  scale_fill_manual(values=c("FDR < 0.1"="#fdae61",
                               "NS"="grey")) +
  scale_shape_manual(values=c(22,24)) +
  scale_size_manual(values=c(2.8,2)) +
  lims(y=c(-plot.lim,plot.lim)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized enrichment score (NES)",
         fill = "FGSEA significance", shape="",
       title="a) RSTR vs LTBI in media   b) RSTR vs LTBI in +Mtb") + 
  facet_grid(group2 ~ group, scales="free_y", 
             labeller = as_labeller(facet_lab),
                space = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        strip.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #force legend color/shape
  guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "black")),
           size = "none") 
  
plot1

plot2<- h.signif.format %>% 
  filter(group %in% c("TBinRSTR","TBinLTBI")) %>% 
  ggplot() +
  geom_segment(aes(x=pathway, xend=pathway, 
                   y=0, yend=fgsea.NES),
               size=0.5) +
  geom_point(aes(x=pathway, y=fgsea.NES,
                 fill = Significance, size = dataset, shape=dataset),
             stroke=0.7) +
  geom_hline(yintercept = 0) +
  
  scale_fill_manual(values=c("FDR < 0.1"="#fdae61",
                             "NS"="grey")) +
  scale_shape_manual(values=c(22,24)) +
  scale_size_manual(values=c(2.8,2)) +
  lims(y=c(-plot.lim,plot.lim)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized enrichment score (NES)",
       fill = "FGSEA significance", shape="",
       title="a) +Mtb vs media in LTBI     b) +Mtb vs media in RSTR") + 
  facet_grid(group2 ~ group, scales="free_y", 
             labeller = as_labeller(facet_lab),
             space = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        strip.text.y = element_blank(),
        strip.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #force legend color/shape
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")),
         size = "none") 

#plot2

#### Save ####
ggsave(filename="Fig3.GSEA.RSTR.pdf",
         plot1,
         width = 6.8, height=4.5)
ggsave(filename="FigS3.GSEA.MTB.pdf",
       plot2,
       width = 6.7, height=4.5)

ggsave(filename="Fig3.GSEA.RSTR.tiff",
       plot1,
       width = 6.8, height=4.5)
ggsave(filename="FigS3.GSEA.MTB.tiff",
       plot2,
       width = 6.7, height=4.5)

