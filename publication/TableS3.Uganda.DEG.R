library(tidyverse)

#### Interaction model ####
#### FDR ####
model1 <- read_csv("../Uganda/results/gene_level/RSTR.Mtb.model.results.anno.csv")

DEG <- model1 %>% 
  filter(grepl(":", variable) & FDR <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

model1.DEG <- model1 %>% 
  select(gene,variable, pval, FDR) %>% 
  filter(gene %in% DEG & variable == "conditionTB:Sample_GroupRSTR") %>% 
  #rename variables
  mutate(variable = recode(variable,
                      "conditionTB:Sample_GroupRSTR"="stimMTB:phenotypeRSTR"))

#### Fold change #### 
load("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")

dat <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("gene") %>% 
  #DEGs
  filter(gene %in% DEG) %>% 
  #add metadata
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets)

#Main term: RSTR
# FC.tb <- dat %>% 
#   #calculate TB/MEDIA fold change
#   select(gene, FULLIDNO, condition, value) %>% 
#   pivot_wider(names_from = condition) %>% 
#   mutate(FC = TB-MEDIA) %>% 
#   #summarise mean FC
#   group_by(gene) %>% 
#   summarise(log2FC = mean(FC, na.rm=TRUE)) %>% 
#   mutate(variable = "stimMTB")

#Main term: MTB
# FC.rstr <- dat %>% 
#   #calculate RSTR/LTBI fold change
#   select(gene, FULLIDNO, Sample_Group, value) %>% 
#   group_by(gene, Sample_Group) %>% 
#   summarise(meanE = mean(value,na.rm=TRUE)) %>% 
#   pivot_wider(names_from = Sample_Group, values_from = meanE) %>% 
#   mutate(FC = RSTR-LTBI) %>% 
#   #summarise mean FC
#   group_by(gene) %>% 
#   summarise(log2FC = mean(FC, na.rm=TRUE)) %>% 
#   mutate(variable = "phenotypeRSTR")

#Interaction term: MTB:RSTR
FC.interact <- dat %>% 
  #calculate TB/MEDIA fold change
  select(gene, FULLIDNO, Sample_Group, condition, value) %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(FC = TB-MEDIA) %>% 
  #summarise mean FC within RSTR groups
  group_by(gene, Sample_Group) %>% 
  summarise(meanFC = mean(FC, na.rm=TRUE)) %>% 
  rename(group=Sample_Group) %>% 
  #Calculate RSTR-LTBI of mean Mtb-media FC
  pivot_wider(names_from = group, values_from = meanFC) %>% 
  rowwise() %>% 
  mutate(log2FC = RSTR-LTBI) %>% 
  select(gene, log2FC) %>% 
  mutate(variable="stimMTB:phenotypeRSTR")

#### Contrast model ####
#### FDR ####
models <- list.files(path="../Uganda/results/gene_level/",
                     pattern=".csv.gz", full.names = TRUE)

model2 <- data.frame()
for(file in models){
  model2 <- bind_rows(model2, read_csv(file))
}

model2.format <- model2 %>% 
  filter(grepl("TB|RSTR",variable)) %>% 
   mutate(subset = recode(group, 
                            "RSTR.Mtb.ltbi"="LTBI",
                            "RSTR.Mtb.media"="MEDIA",
                            "RSTR.Mtb.rstr"="RSTR",
                            "RSTR.Mtb.tb"="MTB")) %>% 
  mutate(variable = recode(variable, "conditionTB"="stimMTB",
                           "Sample_GroupRSTR"="phenotypeRSTR",
                    "conditionTB:Sample_GroupRSTR"="stimMTB:phenotypeRSTR")) %>% 
   select(gene, subset, variable, pval, FDR)

#### Fold change ####
# TB vs MEDIA
FC.contrast.tb <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("gene") %>% 
  #DEGs
  filter(gene %in% DEG) %>% 
  #add metadata
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  #calculate TB/MEDIA fold change
  select(gene, FULLIDNO, Sample_Group, condition, value) %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(FC = TB-MEDIA) %>% 
  #summarise mean FC within RSTR group
  group_by(gene, Sample_Group) %>% 
  summarise(log2FC = mean(FC, na.rm=TRUE)) %>% 
  rename(subset=Sample_Group) %>% 
  mutate(variable="stimMTB")

# RSTR vs LTBI
FC.contrast.rstr <- as.data.frame(dat.combined.voom$E) %>% 
  rownames_to_column("gene") %>% 
  #DEGs
  filter(gene %in% DEG) %>% 
  #add metadata
  pivot_longer(-gene, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  #Caclulate mean expression for RSTR vs LTCI
  select(gene, FULLIDNO, Sample_Group, condition, value) %>% 
  group_by(gene, Sample_Group, condition) %>% 
  summarise(value = mean(value, na.rm=TRUE), .groups="drop") %>% 
  #summarise mean FC within TB condition
  pivot_wider(names_from = Sample_Group) %>% 
  mutate(log2FC = RSTR-LTBI) %>% 
  select(gene, condition, log2FC) %>% 
  rename(subset=condition) %>% 
  mutate(variable="phenotypeRSTR")

#### Combine and save ####
#model1.all <- bind_rows(FC.rstr,FC.tb,FC.interact) %>% 
model1.all <- FC.interact %>% 
  full_join(model1.DEG) %>% 
  select(gene,variable,log2FC,pval,FDR)%>% 
  arrange(FDR)

model2.all <- bind_rows(FC.contrast.rstr,FC.contrast.tb) %>% 
  mutate(subset = recode(subset, "TB"="MTB")) %>% 
  full_join(model2.format) %>% 
  select(gene,subset,variable,log2FC,pval,FDR) %>% 
  arrange(gene,variable,subset)

library(xlsx)
write.xlsx(as.data.frame(model1.all), "TableS3.Uganda.DEG.xlsx", 
           sheetName = "interaction.model", 
           row.names = FALSE)
write.xlsx(as.data.frame(model2.all), "TableS3.Uganda.DEG.xlsx", 
           sheetName = "contrast.model", 
           row.names = FALSE, append = TRUE)
