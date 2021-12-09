#### Setup ####
library(tidyverse)
library(limma)
set.seed(4389)

#### Data ####
load("data_clean/RSTR_SA_dat_clean.RData")

# List all genes to plot
GOI <- read_csv("results/gene_level/SA_RSTR.Mtb.model.results.anno.csv") %>% 
  filter(FDR <= 0.2 & variable %in% c("Sample_Group","condition:Sample_Group")) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

#### Extract expression data ####
dat.GOI <- as.data.frame(dat.abund.norm.voom$E) %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% GOI) %>% 
  pivot_longer(-gene, names_to="libID", values_to="expression") %>% 
  left_join(dat.abund.norm.voom$targets)

#### Format FDR data ####
library(ggpubr)
GOI.fdr <- read_csv("results/gene_level/SA_RSTR.Mtb.media.model.results.csv.gz") %>% 
  bind_rows(read_csv("results/gene_level/SA_RSTR.Mtb.tb.model.results.csv.gz")) %>% 
  bind_rows(read_csv("results/gene_level/SA_RSTR.Mtb.ltbi.model.results.csv.gz")) %>% 
  bind_rows(read_csv("results/gene_level/SA_RSTR.Mtb.rstr.model.results.csv.gz")) %>% 
  filter(gene %in% GOI &
           variable %in% c("Sample_GroupRSTR","conditionTB")) %>% 
  distinct(group, gene, variable, FDR) %>% 
  #Make symbols for plots
  mutate(symbol = ifelse(FDR <= 0.05,"***",
                         ifelse(FDR <= 0.1, "**",
                                ifelse(FDR <= 0.2, "*", NA)))) %>% 
  filter(!is.na(symbol)) %>% 
  #add start end groups
  mutate(group1 = ifelse(group == "SA_RSTR.Mtb.media", 1, 
                         ifelse(group == "SA_RSTR.Mtb.tb", 3,
                                ifelse(group == "SA_RSTR.Mtb.ltbi", 1,
                                       ifelse(group == "SA_RSTR.Mtb.rstr", 2,
                                              NA)))),
         group2 = ifelse(group == "SA_RSTR.Mtb.media", 2, 
                         ifelse(group == "SA_RSTR.Mtb.tb", 4,
                                ifelse(group == "SA_RSTR.Mtb.ltbi", 3,
                                       ifelse(group == "SA_RSTR.Mtb.rstr", 4,
                                              NA)))))

#Add y location for FDR based on max expression in plot
GOI.pe <- dat.GOI %>% 
  #Max expression per gene and experiment
  group_by(gene) %>% 
  summarise(max.e = max(expression, na.rm=TRUE)) %>% 
  ungroup() %>% 
  #Add to FDR data
  right_join(GOI.fdr) %>% 
  arrange(gene, group1, group2)

#first entry per gene
#Set y position to 1
first <- GOI.pe %>% 
  group_by(gene) %>% 
  slice(1) %>% 
  mutate(y.position1 = 1)

#Add first position data back and fill in remaining
GOI.pey <- GOI.pe %>% 
  full_join(first) %>% 
  
  #Fill in positions 2 - N
  group_by(gene) %>% 
  mutate(y.position2 = lag(y.position1)+1,
         y.position3 = lag(y.position2)+1,
         y.position4 = lag(y.position3)+1) %>% 
  #Collapse positions into 1 column
  mutate(y.position = ifelse(!is.na(y.position1),y.position1,
                             ifelse(!is.na(y.position2),y.position2,
                                    ifelse(!is.na(y.position3),y.position3,
                                           ifelse(!is.na(y.position4),y.position4,
                                                  NA))))) %>% 
  #Scale to max expression
  mutate(y.position = (2^max.e)*((y.position+2)/11+0.85)) %>% 
  ungroup()

#### Plot GSEA genes ####
plot.GOI <- dat.GOI %>% 
  droplevels() %>% 
  
  ggplot(aes(x=paste(condition,Sample_Group, sep="\n"), 
             y=2^expression, color=Sample_Group)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=..y.., ymin=..y..),
               color="black", width=0.5) +
  facet_wrap(~gene, scales="free") +
  #Add pval
  stat_pvalue_manual(data=filter(GOI.pey, gene %in% GOI), 
                     label="symbol", xmin="group1", xmax="group2") +
  # #Beautify
  theme_bw() +
  labs(x="", y="Normalized expression (CPM)") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background =element_rect(fill="white"))

#plot.GOI

#### Save ####
ggsave("figs/gene_level/DEG.all.pdf", plot.GOI,
       height=6, width=8)

