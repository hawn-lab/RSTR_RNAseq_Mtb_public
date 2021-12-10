library(tidyverse)
library(matrixTests)
library(broom)
setwd("~/Documents/_Hawn/manuscripts/MS_RSTR_TB.RNAseq/RSTR_RNAseq_Mtb_public/")

#### Uganda ####
#Get all samples from orig and validation sequencing runs
orig <- read_csv("Uganda/data_raw/Hawn_UgandaMacrophage_counts.csv")
val <- read_csv("Uganda/data_raw/Hawn_UgandaValidation_counts.csv")

sampU <- data.frame(sampID = c(colnames(orig)[-1], colnames(val)[-1])) %>% 
  #Mark datasets
  mutate(dataset = ifelse(sampID %in% colnames(orig),"orig","val")) %>% 
  #get donor IDs
  separate(sampID, into=c("RS_SUB_ACCESSION_NO","condition"), 
           sep="_|-", remove = FALSE) %>% 
  distinct(RS_SUB_ACCESSION_NO,dataset)

#Add FULLIDNO and addtl demographic data
metaU <- read_csv("Uganda/data_raw/2020.11.20RSTR_Hawn_metadata.csv") %>% 
  distinct(RS_SUB_ACCESSION_NO, FULLIDNO, Sample_Group, 
           mono.RNAseq, mono.RNAseq.validation,
           M0_KCVAGE, KCHCA_AGE_YR_CURRENT, M0_KCVSEX, 
           avgBMI, HIVSTAT_CURRENT, KCB_BCGSCAR, RISK_SCORE, 
           family.3, family.1) %>% 
  right_join(sampU)

#### Total samples ####
#Original
metaU %>% 
  filter(dataset == "orig") %>% 
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)
#Low RNA or contamination
metaU %>% 
  filter(mono.RNAseq == "FAIL") %>%  
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)
#Pass-filter
metaU %>% 
  filter(mono.RNAseq == "PASS") %>%  
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)


#Validation
metaU %>% 
  filter(dataset == "val") %>% 
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)
#Low RNA or contamination
metaU %>% 
  filter(mono.RNAseq.validation == "FAIL") %>%  
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)
#Pass-filter
metaU %>% 
  filter(mono.RNAseq.validation == "PASS") %>%  
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)

#Combined
metaU %>% 
  filter(mono.RNAseq == "PASS" | mono.RNAseq.validation == "PASS") %>%  
  distinct(FULLIDNO, Sample_Group) %>% 
  count(Sample_Group)

#### Continuous vars ####
#final samples
attach("Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")
sampU.final <- dat.combined.voom$targets

metaU.final <- metaU %>% 
  filter(mono.RNAseq == "PASS" | mono.RNAseq.validation == "PASS" &
           RS_SUB_ACCESSION_NO %in% sampU.final$RSID) %>% 
  select(-mono.RNAseq, -mono.RNAseq.validation, -dataset, -RS_SUB_ACCESSION_NO) %>% 
  distinct()
  
## median + IQR + mean + sd
metaU.final %>% 
  group_by(Sample_Group) %>% 
  summarise(across(c(M0_KCVAGE, KCHCA_AGE_YR_CURRENT, avgBMI, RISK_SCORE), 
                   .fns=list(median=median, IQR=IQR, mean=mean, sd=sd), 
                   na.rm=TRUE)) %>% 
  as.matrix()

## NAs
metaU.final %>% 
  group_by(Sample_Group) %>% 
  summarise(across(c(M0_KCVAGE, KCHCA_AGE_YR_CURRENT, avgBMI, RISK_SCORE), 
                   ~ifelse(any(is.na(.)), sum(is.na(.)), 0))) %>% 
  rename_if(is.numeric, ~paste(., "NA", sep="_"))

## Mann-Whitney
wilcox.x <- metaU.final %>% 
  filter(Sample_Group == "LTBI") %>% 
  select(M0_KCVAGE, KCHCA_AGE_YR_CURRENT, avgBMI, RISK_SCORE) 
wilcox.y <- metaU.final %>% 
  filter(Sample_Group == "RSTR") %>% 
  select(M0_KCVAGE, KCHCA_AGE_YR_CURRENT, avgBMI, RISK_SCORE) 

col_wilcoxon_twosample(wilcox.x, wilcox.y, alternative = "two.sided",
                       mu = 0, exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#### Categorical vars ####
## % + counts
count(metaU.final, Sample_Group, M0_KCVSEX) %>% 
  #Percent male
  group_by(Sample_Group) %>% 
  mutate(perc = n[M0_KCVSEX=="M"]/sum(n[!is.na(M0_KCVSEX)]))

count(metaU.final, Sample_Group, HIVSTAT_CURRENT)

count(metaU.final, Sample_Group, KCB_BCGSCAR) %>% 
  #Percent Y
  filter(!is.na(KCB_BCGSCAR)) %>% 
  group_by(Sample_Group) %>% 
  mutate(perc = n[KCB_BCGSCAR=="Y"]/sum(n))

## Chi-squared
chisq.test(table(metaU.final$Sample_Group, metaU.final$M0_KCVSEX)) %>% 
  tidy()

chisq.test(table(metaU.final$Sample_Group, metaU.final$M0_KCVSEX)) %>% 
  tidy()

#### Family groups ####
fam3 <- metaU.final %>% 
  #fillin missing
  mutate(family.3 = ifelse(is.na(family.3), strsplit(FULLIDNO, split = "-")[[1]][1],
                           family.3)) %>% 
  group_by(Sample_Group, family.3) %>% 
  mutate(sum.3 = n()-1) %>% #total family groups size in dataset minus self
  ungroup()

fam3 %>% 
  group_by(Sample_Group) %>% 
  summarise(mean = mean(sum.3),
            sd=sd(sum.3),
            median=median(sum.3),
            IQR=IQR(sum.3))

fam1 <- metaU.final %>% 
  #fillin missing
  mutate(family.1 = ifelse(is.na(family.1), 
                           paste(strsplit(FULLIDNO, split = "-")[[1]][1],
                                 "A1", sep=""),
                           family.1),
         family.3 =as.character(family.3)) %>% 
  group_by(Sample_Group, family.1) %>% 
  mutate(sum.1 = n()-1) %>%  #total family groups size in dataset minus self
  ungroup()

fam1 %>% 
  group_by(Sample_Group) %>% 
  summarise(mean = mean(sum.1),
            sd=sd(sum.1),
            median=median(sum.1),
            IQR=IQR(sum.1))

## Mann-Whitney
wilcox.x <- full_join(fam1,fam3) %>% 
  filter(Sample_Group == "LTBI") %>% 
  select(sum.1, sum.3) 
wilcox.y <- full_join(fam1,fam3) %>% 
  filter(Sample_Group == "RSTR") %>% 
  select(sum.1, sum.3)

col_wilcoxon_twosample(wilcox.x, wilcox.y, alternative = "two.sided",
                       mu = 0, exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#### South Africa ####
#Get all samples from orig and validation sequencing runs
SA <- read_tsv("SAfrica/data_raw/Hawn_SAMacrophage_counts.txt")

sampSA <- data.frame(sampID = colnames(SA)[-1]) %>% 
  #get donor IDs
  separate(sampID, into=c("ptID","condition"), 
           sep="_", remove = FALSE) %>% 
  distinct(ptID)

#Add FULLIDNO and addtl demographic data
metaSA <- read_csv("SAfrica/data_raw/2020.11.05SA.RSTR_Hawn_metadata.csv") %>% 
  #Fill in missing RSTR
  mutate(Sample_Group = ifelse(is.na(Sample_Group) & status == "control",
                               "LTBI",
                               ifelse(is.na(Sample_Group) & status == "case",
                                      "RSTR", Sample_Group))) %>% 
  distinct(ptID, Sample_Group, 
           mono.RNAseq, 
           age, gender, 
           bmi, YearsWorkedUnderground, ethnic, bcgscar) %>% 
  right_join(sampSA)

#### Total samples ####
#Original
metaSA %>% 
  distinct(ptID, Sample_Group) %>% 
  count(Sample_Group)
#Low RNA or contamination
metaSA %>% 
  filter(is.na(mono.RNAseq)) %>%  
  distinct(ptID, Sample_Group) %>% 
  count(Sample_Group)
#Pass-filter
metaSA %>% 
  filter(mono.RNAseq == "PASS") %>%  
  distinct(ptID, Sample_Group) %>% 
  count(Sample_Group)
#nonWhite
metaSA %>% 
  filter(mono.RNAseq == "PASS" & !grepl("White", ethnic)) %>%  
  distinct(ptID, Sample_Group) %>% 
  count(Sample_Group)

#### Continuous vars ####
metaSA.final <- metaSA %>% 
  filter(mono.RNAseq == "PASS") %>% 
  distinct()

## median + IQR + mean + sd
metaSA.final %>% 
  group_by(Sample_Group) %>% 
  summarise(across(c(age, YearsWorkedUnderground, bmi), 
                   .fns=list(median=median, IQR=IQR), 
                   na.rm=TRUE)) %>% 
  as.matrix()

## NAs
metaSA.final %>% 
  group_by(Sample_Group) %>% 
  summarise(across(c(age, YearsWorkedUnderground, bmi), 
                   ~ifelse(any(is.na(.)), sum(is.na(.)), 0))) %>% 
  rename_if(is.numeric, ~paste(., "NA", sep="_"))

## Mann-Whitney
wilcox.x <- metaSA.final %>% 
  filter(Sample_Group == "LTBI") %>% 
  select(age, YearsWorkedUnderground, bmi) 
wilcox.y <- metaSA.final %>% 
  filter(Sample_Group == "RSTR") %>% 
  select(age, YearsWorkedUnderground, bmi) 

col_wilcoxon_twosample(wilcox.x, wilcox.y, alternative = "two.sided",
                       mu = 0, exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#### Categorical vars ####
## % + counts
count(metaSA.final, Sample_Group, gender) 

count(metaSA.final, Sample_Group, bcgscar) %>% 
  #Percent Y
  filter(!is.na(bcgscar)) %>% 
  group_by(Sample_Group) %>% 
  mutate(perc = n[bcgscar=="Yes"]/sum(n))

count(metaSA.final, Sample_Group, ethnic) %>% 
  group_by(Sample_Group) %>% 
  mutate(perc = n/sum(n))


## Chi-squared
chisq.test(table(metaSA.final$Sample_Group, metaSA.final$bcgscar)) %>% 
  tidy()

chisq.test(table(metaSA.final$Sample_Group, metaSA.final$ethnic)) %>% 
  tidy()
