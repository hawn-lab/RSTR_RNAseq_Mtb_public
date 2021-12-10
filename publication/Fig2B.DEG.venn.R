library(tidyverse)
library(venn)
fdr.cutoff <- 0.2

#### Load data ####
U <- read_csv("../Uganda/results/gene_level/RSTR.Mtb.model.results.anno.csv")
SA <- read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.model.results.anno.csv")

contrastU <- read_csv("../Uganda/results/gene_level/RSTR.Mtb.media.model.results.csv.gz") %>% 
  bind_rows(read_csv("../Uganda/results/gene_level/RSTR.Mtb.tb.model.results.csv.gz")) %>% 
  bind_rows(read_csv("../Uganda/results/gene_level/RSTR.Mtb.ltbi.model.results.csv.gz")) %>% 
  bind_rows(read_csv("../Uganda/results/gene_level/RSTR.Mtb.rstr.model.results.csv.gz")) 

contrastSA <- read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.media.model.results.csv.gz") %>% 
  bind_rows(read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.tb.model.results.csv.gz")) %>% 
  bind_rows(read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.ltbi.model.results.csv.gz")) %>% 
  bind_rows(read_csv("../SAfrica/results/gene_level/SA_RSTR.Mtb.rstr.model.results.csv.gz")) 

#### Uganda: Interaction model ####
venn.ls <- list()
venn.ls[["RSTR - LTBI"]] <- U %>% 
  filter(model=="lmekin" & variable == "Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["+Mtb - media"]] <- U %>% 
  filter(model=="lmekin" & variable == "conditionTB" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["Interaction"]] <- U %>% 
  filter(model=="lmekin" & variable == "conditionTB:Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

#### Uganda: Contrast model ####
venn.ls2 <- list()
venn.ls2[["RSTR - LTBI\nin media"]] <- contrastU %>% 
  filter(group=="RSTR.Mtb.media" & variable == "Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls2[["RSTR - LTBI\nin +Mtb"]] <- contrastU %>% 
  filter(group=="RSTR.Mtb.tb" & variable == "Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls2[["+Mtb - media\nin LTBI"]] <- contrastU %>% 
  filter(group=="RSTR.Mtb.ltbi" & variable == "conditionTB" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls2[["+Mtb - media\nin RSTR"]] <- contrastU %>% 
  filter(group=="RSTR.Mtb.rstr" & variable == "conditionTB" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

#### SAfrica: Interaction model ####
venn.ls3 <- list()
venn.ls3[["RSTR - LTBI"]] <- SA %>% 
  filter(model=="lme" & variable == "Sample_Group" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls3[["+Mtb - media"]] <- SA %>% 
  filter(model=="lme" & variable == "condition" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls3[["Interaction"]] <- SA %>% 
  filter(model=="lme" & variable == "condition:Sample_Group" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

#### Uganda: Contrast model ####
venn.ls4 <- list()
venn.ls4[["RSTR - LTBI\nin media"]] <- contrastSA %>% 
  filter(group=="SA_RSTR.Mtb.media" & variable == "Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls4[["RSTR - LTBI\nin +Mtb"]] <- contrastSA %>% 
  filter(group=="SA_RSTR.Mtb.tb" & variable == "Sample_GroupRSTR" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls4[["+Mtb - media\nin LTBI"]] <- contrastSA %>% 
  filter(group=="SA_RSTR.Mtb.ltbi" & variable == "conditionTB" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn.ls4[["+Mtb - media\nin RSTR"]] <- contrastSA %>% 
  filter(group=="SA_RSTR.Mtb.rstr" & variable == "conditionTB" &
           FDR <= fdr.cutoff) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

#### Plot ####
pdf("Fig2B.DEG.venn.unformatted.pdf", width=7, height=6)

par(mfrow=c(2,2))

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls, box=FALSE)
title(sub="Uganda\nMtb*RSTR interaction model", line = -1)
title(adj=0, main="b)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls2, box=FALSE)
title(sub="Uganda\nMtb:RSTR contrast model", line = -1)
title(adj=0, main="d)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls3, box=FALSE)
title(sub="South Africa\nMtb*RSTR interaction model", line = -1)
title( adj=0, main="c)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls4, box=FALSE)
title(sub="South Africa\nMtb:RSTR contrast model", line = -1)
title(adj=0, main="e)", cex=1.1, col="black", font=1, line=-1)

dev.off()

tiff("Fig2B.DEG.venn.unformatted.tiff", width=7, height=6, units = "in", res=300)

par(mfrow=c(2,2))

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls, box=FALSE)
title(sub="Uganda\nMtb*RSTR interaction model", line = -1)
title(adj=0, main="b)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls2, box=FALSE)
title(sub="Uganda\nMtb:RSTR contrast model", line = -1)
title(adj=0, main="d)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls3, box=FALSE)
title(sub="South Africa\nMtb*RSTR interaction model", line = -1)
title( adj=0, main="c)", cex=1.1, col="black", font=1, line=-1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.ls4, box=FALSE)
title(sub="South Africa\nMtb:RSTR contrast model", line = -1)
title(adj=0, main="e)", cex=1.1, col="black", font=1, line=-1)

dev.off()

