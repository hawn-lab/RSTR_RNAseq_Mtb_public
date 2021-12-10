library(tidyverse)
library(venn)

#### Data ####
HU <- read_csv("../Uganda/results/GSEA/h_GSEA.result.csv")

edgeU <- HU %>% 
  dplyr::select(group, pathway, fgsea.leadingEdge) %>% 
  separate(fgsea.leadingEdge, into=as.character(c(1:200)), sep=";") %>% 
  pivot_longer(as.character(c(1:200)), values_to = "gene", names_to = "rank") %>% 
  drop_na(gene) 

HSA <- read_csv("../SAfrica/results/GSEA/h_GSEA.result.csv")

# edgeSA <- HSA %>% 
#   dplyr::select(group, pathway, fgsea.leadingEdge) %>% 
#   separate(fgsea.leadingEdge, into=as.character(c(1:200)), sep=";") %>% 
#   pivot_longer(as.character(c(1:200)), values_to = "gene", names_to = "rank") %>% 
#   drop_na(gene) %>% 
#   mutate(DEG = NA)

#genes in sets
myGO <- fgsea::gmtPathways("../Uganda/data_clean/Broad.gene.sets/h.all.v7.2.symbols.gmt")

#### Uganda ####
venn.IFN1 <- list()
temp <- edgeU %>% 
  filter(group == "RSTRinMEDIA" & pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  distinct(gene)
venn.IFN1[["IFNg genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)
venn.IFN1[["IFNa genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)

venn.IFN2 <- list()
temp <- edgeU %>% 
  filter(group == "RSTRinTB" & pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% 
  distinct(gene)
venn.IFN2[["IFNg genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)
venn.IFN2[["IFNa genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)

venn.IFN3 <- list()
temp <- edgeU %>% 
  filter(group == "RSTRinMEDIA" & pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% 
  distinct(gene)
venn.IFN3[["IFNg genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)
venn.IFN3[["IFNa genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)

venn.IFN4 <- list()
temp <- edgeU %>% 
  filter(group == "RSTRinTB" & pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% 
  distinct(gene)
venn.IFN4[["IFNg genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)
venn.IFN4[["IFNa genes"]] <- temp %>% 
  filter(gene %in% myGO[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]]) %>% 
  unlist(use.names = FALSE)

#### Plot ####

pdf("FigS4.GSEA.IFN.LE.venn.pdf", width=4, height=5)
par(mfrow=c(2,2))

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN1, box=FALSE)
title(sub="61% (54/88)\nRSTR - LTBI leading-edge\nspecific to IFNg\nUganda media", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN2, box=FALSE)
title(sub="77% (30/39)\nRSTR - LTBI leading-edge\nspecific to IFNg\nUganda +Mtb", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN3, box=FALSE)
title(sub="28% (13/47)\nRSTR - LTBI leading-edge\nspecific to IFNa\nUganda media", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN4, box=FALSE)
title(sub="28% (9/32)\nRSTR - LTBI leading-edge\nspecific to IFNa\nUganda +Mtb", line = -1)

dev.off()

tiff("FigS4.GSEA.IFN.LE.venn.tiff", width=4, height=5, units = "in", res=300)
par(mfrow=c(2,2))

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN1, box=FALSE)
title(sub="61% (54/88)\nRSTR - LTBI leading-edge\nspecific to IFNg\nUganda media", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN2, box=FALSE)
title(sub="77% (30/39)\nRSTR - LTBI leading-edge\nspecific to IFNg\nUganda +Mtb", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN3, box=FALSE)
title(sub="28% (13/47)\nRSTR - LTBI leading-edge\nspecific to IFNa\nUganda media", line = -1)

venn(ilab=FALSE, ilcs=1, sncs=1, x=venn.IFN4, box=FALSE)
title(sub="28% (9/32)\nRSTR - LTBI leading-edge\nspecific to IFNa\nUganda +Mtb", line = -1)

dev.off()