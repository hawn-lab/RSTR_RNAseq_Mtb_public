#### Setup ####
library(tidyverse)
library(ggtext)
library(ggrepel)
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
                                                "Negative","Negative leading edge"))) %>% 
  #facet labels
  mutate(facet.group = recode_factor(factor(group),
                                   'MEDIA'="RSTR vs LTBI in media<br><span style='font-size:9pt'>  Down in RSTR <- -> Up in RSTR</span>",
                                   'TB'="RSTR vs LTBI in +Mtb<br><span style='font-size:9pt'>  Down in RSTR <- -> Up in RSTR</span>",
                                   'LTBI'="+Mtb vs media in LTBI<br><span style='font-size:9pt'>  Down in +Mtb <- -> Up in +Mtb</span>",
                                   'RSTR'="+Mtb vs media in RSTR<br><span style='font-size:9pt'>  Down in +Mtb <- -> Up in +Mtb</span>"))

#### Plot parameters ####
#Set jitter
pos <- position_jitter(width = 0.3, seed = 589, height = 0)

#### Plot ####
plot1 <- h.toPlot.lab2 %>% 
  filter(dataset=="Uganda") %>% 
  ggplot(aes(x = pathway, y = meanFC)) +
  geom_violin() +
  #Add non-labeled points with jitter
  geom_jitter(data=filter(h.toPlot.lab2, is.na(DEG)),
              aes(color=col.group), position = pos, alpha=0.5) +
  #Add labeled points
  geom_point(data=filter(h.toPlot.lab2, !is.na(DEG)),
             aes(color=col.group), alpha=0.5) +
  #Add DEG labels
  geom_text_repel(data=filter(h.toPlot.lab2, !is.na(DEG)),
                  aes(label=rowname, color=col.group), direction="both",
                  nudge_x=-0.3, min.segment.length = unit(0, 'lines'),
                  show.legend = FALSE, size=3, max.overlaps = 100, segment.size=0.3) +
  facet_grid(group2 ~ facet.group, 
             scales="free",space="free_y") +
  coord_flip() +
  #Beautify
  theme_bw() +
  labs(x="", y="Gene mean log2 fold change") +
  scale_color_manual(values=c("Positive leading edge"="#b2182b","Positive"="#d6604d",
                              "Negative"="#4393c3","Negative leading edge"="#053061"),
                      name="") +
  scale_y_continuous(limits = function(x){c(min(min(x),-max(x)), 
                                            max(-min(x),max(x)))}) +
  theme(panel.grid.major.y = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_markdown(size = 11, lineheight = 1.1))
#plot1

#### Save ####
ggsave("FigSX.GSEA.genes.Uganda.pdf", plot1, 
       width = 20, height = 12)

