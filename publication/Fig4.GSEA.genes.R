#### Setup ####
library(tidyverse)
library(ggrepel)
library(cowplot)
set.seed(4389)
fdr.cutoff <- 0.2

#### Data ####
#List genes in each Hallmark term of interest
myGO <- fgsea::gmtPathways("../Uganda/data_clean/Broad.gene.sets/h.all.v7.2.symbols.gmt")
GO.df <- plyr::ldply(myGO, rbind) %>% 
  rename(term = `.id`) %>% 
  pivot_longer(-term, values_to="hgnc_symbol") %>% 
  select(-name)

#### Fold change ####
attach("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")
voomU <- dat.combined.voom

attach("../SAfrica/data_clean/RSTR_SA_dat_clean.RData")
voomSA <- dat.abund.norm.voom

for(dat in c("voomU", "voomSA")){
  # TB vs MEDIA
  FC.tb <- as.data.frame(get(dat)$E) %>% 
    #add metadata
    rownames_to_column() %>% 
    pivot_longer(-rowname, names_to = "libID") %>% 
    left_join(get(dat)$targets) %>% 
    #calculate TB/MEDIA fold change
    select(rowname, FULLIDNO, Sample_Group, condition, value) %>% 
    pivot_wider(names_from = condition) %>% 
    mutate(FC = TB-MEDIA) %>% 
    #summarise mean FC within RSTR group
    group_by(rowname, Sample_Group) %>% 
    summarise(meanFC = mean(FC, na.rm=TRUE)) %>% 
    rename(group=Sample_Group)
  
  # RSTR vs LTBI
  FC.rstr <- as.data.frame(get(dat)$E) %>% 
    #add metadata
    rownames_to_column() %>% 
    pivot_longer(-rowname, names_to = "libID") %>% 
    left_join(get(dat)$targets) %>% 
    #Caclulate mean expression for RSTR vs LTCI
    select(rowname, FULLIDNO, Sample_Group, condition, value) %>% 
    group_by(rowname, Sample_Group, condition) %>% 
    summarise(value = mean(value, na.rm=TRUE), .groups="drop") %>% 
    #summarise mean FC within TB condition
    pivot_wider(names_from = Sample_Group) %>% 
    mutate(meanFC = RSTR-LTBI) %>% 
    select(rowname, condition, meanFC) %>% 
    rename(group=condition)
  
  #combine for facet plot
  FC <- bind_rows(FC.tb, FC.rstr) 
  
  assign(paste(dat,"FC",sep="."), FC, envir=.GlobalEnv)
}
# Significant terms to include (determined in GSEA.R)
load("GSEA.signif.RData")

#### Format data ####
#Select genes of interest and format
h.toPlot <- mutate(voomU.FC, dataset="Uganda") %>% 
  bind_rows(mutate(voomSA.FC, dataset="South Africa")) %>% 
  #Genes in Hallmark terms
  inner_join(GO.df, by=c("rowname"="hgnc_symbol")) %>% 
  #hallmark terms of interest 
  filter(term %in% h.signif$pathway) %>% 
  #direction color groups
  mutate(col.group = ifelse(meanFC <0, "down", "up")) %>% 
  #Beautify labels
  mutate(pathway = gsub("HALLMARK_","",term)) %>% 
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
                            "coagulation", #
                            
                            #Split R
                            "epithelial mesenchymal transition",#
                            "hypoxia",#
                            "KRAS signaling down",#
                            "IL2 STAT5 signaling",#
                            #Down to up R
                            "inflammatory response",#
                            "TNFA signaling via NF-kB",#
                            
                            #All R down
                            "IFN alpha response",#
                            "allograft rejection",#
                            "IFN gamma response",
                            
                            #Down with TB
                            "adipogenesis",#
                            "oxidative phosphorylation"))) 

#### List DEGs to label ####
degU <- read_csv("../Uganda/results/model_selection/RSTR.interaction_age.sex.batch.model.results.csv.gz") %>% 
  filter(model=="lmekin" & FDR <= fdr.cutoff & 
           variable == "conditionTB:Sample_GroupRSTR") %>% 
  distinct(gene) %>% 
  mutate(DEG = "y", dataset = "Uganda")

degSA <- read_csv("../SAfrica/results/model_selection/SA_RSTR.interaction_age.model.results.csv.gz") %>% 
  filter(grepl("Sample_Group", variable) & FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% 
  mutate(DEG = "y", dataset = "South Africa")

#add to data
h.toPlot.lab <- h.toPlot %>% 
  left_join(degU, by = c("rowname"="gene", "dataset")) %>% 
  left_join(degSA, by = c("rowname"="gene", "dataset", "DEG"))%>% 
  ungroup()
## Note: none of the SA DEGs are in Hallmark terms

#### Leading edge labels ####
edgeU <- read_csv("../Uganda/results/GSEA/h_GSEA.result.csv") %>% 
  dplyr::select(group, pathway, fgsea.leadingEdge) %>% 
  separate(fgsea.leadingEdge, into=as.character(c(1:200)), sep=";") %>% 
  pivot_longer(as.character(c(1:200)), values_to = "rowname", names_to = "rank") %>% 
  drop_na(rowname) %>% 
  mutate(group = gsub("TBin|RSTRin","",group)) %>% 
  mutate(dataset = "Uganda") %>% 
  mutate(LE = "Leading edge") %>% 
  select(-rank) %>% 
  rename(term=pathway)

edgeSA <- read_csv("../SAfrica//results/GSEA/h_GSEA.result.csv") %>% 
  dplyr::select(group, pathway, fgsea.leadingEdge) %>% 
  separate(fgsea.leadingEdge, into=as.character(c(1:200)), sep=";") %>% 
  pivot_longer(as.character(c(1:200)), values_to = "rowname", names_to = "rank") %>% 
  drop_na(rowname) %>% 
  mutate(group = gsub("TBin|RSTRin","",group))  %>% 
  mutate(dataset = "South Africa") %>% 
  mutate(LE = "Leading edge") %>% 
  select(-rank) %>% 
  rename(term=pathway)

#Groups for color by leading edge

h.toPlot.lab2 <- h.toPlot.lab %>% 
  left_join(full_join(edgeU, edgeSA)) %>% 
  mutate(col.group = ifelse(meanFC < 0 & !is.na(LE), "Negative leading edge",
                            ifelse(meanFC > 0 & !is.na(LE), "Positive leading edge",
                                   ifelse(meanFC < 0 & is.na(LE), "Negative",
                                          ifelse(meanFC > 0 & is.na(LE), "Positive",
                                                 NA))))) %>% 
  mutate(col.group = factor(col.group, levels=c("Positive leading edge","Positive",
                                                "Negative","Negative leading edge")))

#### Plot parameters ####
#Set jitter
pos <- position_jitter(width = 0.3, seed = 589, height = 0)
facet_lab <- c('MEDIA'="Down in RSTR <-  -> Up in RSTR",
               'TB'="Down in RSTR <-  -> Up in RSTR",
               'LTBI'="Down in +Mtb <-  -> Up in +Mtb",
               'RSTR'="Down in +Mtb <-  -> Up in +Mtb")

#### TNF Plot ####
# #Overlap LE
# temp <- edgeU %>% 
#   filter(term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" &
#            group %in% c("MEDIA","TB")) %>% 
#   rename("U"="LE") %>% 
#   select(group, term, rowname, U)
# 
# connect <- edgeSA %>% 
#   filter(term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" &
#            group %in% c("MEDIA","TB")) %>% 
#   rename("SA"="LE") %>% 
#   select(group, term, rowname, SA) %>% 
#   full_join(temp) %>% 
#   filter(!is.na(U) & !is.na(SA)) %>% 
#   select(group, rowname) %>% 
#   mutate(connect=rowname)


tnf.toPlot <- h.toPlot.lab2 %>% 
  filter(term == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" &
           group %in% c("MEDIA","TB")) 

plot2 <- tnf.toPlot %>% 
  ggplot(aes(x = dataset, y = meanFC)) +
  geom_violin() +
  #Add non-labeled points with jitter
  geom_jitter(data=filter(tnf.toPlot, is.na(DEG) ),
              aes(color=col.group), position = pos, alpha=0.5) +
  #Add labeled points
  geom_point(data=filter(tnf.toPlot, !is.na(DEG) ),
             aes(color=col.group), alpha=0.5) +
  #Add DEG labels
  geom_text_repel(data=filter(tnf.toPlot, !is.na(DEG)),
                  aes(label=rowname, color=col.group), direction="both",
                  nudge_x=1.2, min.segment.length = unit(0, 'lines'),
                  show.legend = FALSE, size=3, max.overlaps = 100,
                  segment.size=0.3) +
  facet_grid(~ group,
             labeller = as_labeller(facet_lab)) +
  coord_flip(ylim = c(-1.21, 1.21)) +
  scale_x_discrete(expand = expansion(mult = c(0, 2))) +
  #Beautify
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="", y="Gene mean log2 fold change",
       title="a) RSTR vs LTBI in media                  b) RSTR vs LTBI in +Mtb") +
  scale_color_manual(values=c("Positive leading edge"="#67001f","Positive"="#d6604d",
                              "Negative"="#4393c3","Negative leading edge"="#053061"),
                     name="") 
#plot2

#### Genes plot ####
datU <- voomU$E %>% 
  rownames_to_column("gene")%>% 
  filter(gene %in% c("ABCA1","DUSP2","NR4A2","TNFAIP3")) %>% 
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(select(voomU$targets, libID, Sample_Group)) %>% 
  mutate(dataset="Uganda") %>% 
  separate(libID, into=c("ID","condition"), sep="_", remove=FALSE) %>% 
  mutate(condition = recode_factor(factor(condition), "MEDIA"="media","TB"="+Mtb"),
         x = paste(condition, Sample_Group, sep="\n"),
         x = factor(x, levels=c("media\nLTBI","media\nRSTR","+Mtb\nLTBI","+Mtb\nRSTR"))) %>% 
  mutate(dataset = factor(dataset, levels=c("Uganda","South Africa")))

plot1 <- datU %>% 
  mutate(lab = paste(condition, Sample_Group, sep='\n'),
         lab = factor(lab, levels=c("media\nLTBI","media\nRSTR",
                                    "+Mtb\nLTBI","+Mtb\nRSTR"))) %>% 
  ggplot(aes(x=lab, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color=Sample_Group), position=position_jitterdodge(),
             size=1) +
  theme_classic() +
  facet_wrap(~gene, scales="free", nrow=1) +
  labs(y="Normalized log2 expression", x="", color="") +
  theme(legend.position = "none")
#plot1

#### Save ####
ggsave("Fig4.GSEA.TNFA.genes.pdf", plot_grid(plot2, plot1, 
                                             labels = c("","c)"), nrow=2,
                                             label_fontface = "plain"), 
       width = 9, height = 6)
tiff("Fig4.GSEA.TNFA.genes.tiff",width = 9, height = 6, units = "in", res=300)
plot_grid(plot2, plot1, 
          labels = c("","c)"), nrow=2,
          label_fontface = "plain")
dev.off()
