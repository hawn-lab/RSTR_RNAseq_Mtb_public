library(tidyverse)
#Gene ontology
library(topGO) 
library(org.Hs.eg.db)
#STRING database
require(STRINGdb)
#Working with network objects
require(igraph) #vertex_attr()
#Graphing networks
require(ggraph)
require(scatterpie) #geom_scatterpie()
require(ggnetwork) #theme_blank()
require(scales)

# source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/STRING_network_fxn.R")
set.seed(8434)

#### Load data ####
#List DEGs
U <- read_csv("../Uganda/results/model_selection/RSTR.interaction_age.sex.batch.model.results.csv.gz") %>% 
  filter(model=="lmekin" & variable == "conditionTB:Sample_GroupRSTR") %>% 
  distinct(gene, FDR) 

#STRING db
string_db <- STRINGdb$new(version="11", species=9606,
                          score_threshold=700, input_directory="")

##### DEGs large STRING cluster #####
#DEGs in cluster
lrg <- c("CLNS1A","GEMIN5","GEMIN6","POLDIP3","SARNP","NUP210","HDAC6",
         "PHC2","SUMO1","EGFL8","GPSM3","GPR18",
         "P2RY13","CXCL9","EBF1","ZBTB17","IRF8",
         "FCGR1A","FCGR1B","CIITA","IRF1","IFNG",
         "ANPEP","SMARCB1","CSF1R","IL2",
         "RASGRP1","GRB2","GAB2","PIK3CD","VAV2",
         "AKT1","LIMS1","ARAF","PPP2R3C","PPP2CA",
         "TJP1","KIFC3","CTNNB1","AKT3","FOXO3",
         "CRTC2","TNFRSF1B","IKBKB","RAB3C","RAB5A",
         "RUFY1","DENND1C","RAB13","EXOC8","LEPROTL1",
         "C2CD5","ARAP1","AKAP13","AKAP1","KDELR1","TRAF3",
         "GNA12","NFKBIE","TNFAIP3","FCAR","TIMP1",
         "COPB1","TBK1","FCAMR","PEF1","TBKBP1",
         "CHST7","ADAMTS1","CYB5R3","VCAN","DNAJC3",
         "CST3","CFP","CANT1","UMPS","MTHFR")

##### Format color db with selected terms #####
#Custom groups
GOID.signif <- list()
#lrg
GOID.signif[["inflam"]] <- c("GO:0050863",
                             "GO:0060333",
                             "GO:0030217")
GOID.signif[["leuk"]] <- c("GO:0002444",
                           "GO:0002275",
                           "GO:0043299")
GOID.signif[["mig"]] <- c("GO:0045216",
                          "GO:0001776",
                          "GO:0046579")
GOID.signif[["other"]] <- c("GO:0051043",
                            "GO:0033138",
                            "GO:1901300")
GOID.signif[["PG"]] <- c("GO:0030166")
GOID.signif[["TNF"]] <- c("GO:0001817",
                          "GO:0051607",
                          "GO:0033209")

#Make topGO object
genes.vec <- U %>% 
  mutate(group = ifelse(gene %in% lrg, 1, 0),
         group = factor(group)) %>% 
  dplyr::select(group) %>% unlist(use.names=FALSE)
names(genes.vec) <- U$gene

GOdata <- new("topGOdata", ontology = "BP",
              allGenes = genes.vec,
              annot = annFUN.org, mapping="org.Hs.eg.db", ID = "symbol",
              nodeSize = 1)

#Get genes in terms of interest
db.GOBP <- plyr::ldply(genesInTerm(GOdata), data.frame, .id = "GOID") %>% 
  dplyr::rename("gene"="X..i..")

# size terms
db.size <- db.GOBP %>% 
  distinct(GOID,gene) %>% 
  count(GOID)

#Rank order all term groups
GOID.signif.df <- plyr::ldply(GOID.signif, data.frame, .id="group") %>% 
  rename("X..i.."="GOID") %>% 
  inner_join(db.size) %>% 
  arrange(-n) %>% 
  group_by(group) %>% 
  mutate(ID = row_number())

#Get most specific term per group for each gene
top.terms <- db.GOBP  %>% 
  #signif terms
  right_join(GOID.signif.df) %>%
  #DEGs
  filter(gene %in% lrg) %>%
  #Get most specific term per gene
  group_by(group, gene) %>%
  summarise(ID = max(ID))

# get GOID Descriptions
GO.term <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename("Description"="Term(GOTERM)") %>% 
  mutate(Description = gsub("_|-|[/]"," ",Description))

#Put it all together and force into format for fxn
db.GOBP.format <- db.GOBP %>% 
  #signif terms
  right_join(GOID.signif.df) %>%
  #DEGs
  filter(gene %in% lrg)%>% 
  #add desc
  inner_join(top.terms) %>% 
  inner_join(GO.term) %>% 
  #group other category
  mutate(Description = ifelse(group =="other","other", Description),
         GOID = ifelse(group =="other","other", GOID),
         ID = ifelse(group =="other",1, ID)) %>% 
  #format for fxn
  group_by(group, ID, GOID, Description) %>% 
  summarise(SYMBOLs = paste(unique(gene), collapse="/"),
            size.overlap.term=n()) %>%
  ungroup() %>% 
  mutate(p.adjust=0) %>%
  mutate(Description = factor(Description)) %>%
  arrange(Description) %>% 
  #order
  arrange(group, ID)

#Factor order within group
db.GOBP.format2 <- db.GOBP.format %>% 
  mutate(Description = factor(Description, levels=db.GOBP.format$Description),
         Description = fct_relevel(Description, "other", after=Inf)) %>% 
  arrange(Description) %>% 
  mutate(color = c('#fdd0a2','#fdae6b','#fd8d3c',#inflam orange
                   '#f1b6da', #leuk pink
                   '#bdd7e7','#6baed6','#3182bd',#mig blue
                   '#ffffb2', #PG yellow
                   '#bae4b3','#74c476','#31a354', #TNF green
                   '#969696' #other
         ))

#### Plot LARGE ####
#Format gene vector to matrix
genes.mat <- as.matrix(lrg)
colnames(genes.mat) <- "gene"

#Map genes to STRING
map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)

#Get enrichments
col.mat <- db.GOBP.format2 %>% 
  dplyr::select(Description, size.overlap.term, p.adjust,SYMBOLs)

#Format enrichment results for scatterpie plotting
col.mat.format <- col.mat %>% 
  #Split gene lists within terms
  dplyr::mutate(gene = strsplit(as.character(SYMBOLs), "/")) %>% 
  tidyr::unnest(gene) %>%
  #spread terms to columns
  distinct(gene, Description) %>% 
  dplyr::mutate(value=1) %>% 
  #Add string ID
  left_join(map, by = "gene") %>% 
  dplyr::select(-gene) %>% 
  dplyr::distinct() %>% 
  #Calculate terms per ID
  group_by(STRING_id) %>% 
  dplyr::mutate(total = sum(value)) %>%
  ungroup() %>% 
  #Calculate proportions within terms
  arrange(match(Description,unique(db.GOBP.format2$Description))) %>% 
  pivot_wider(names_from = Description, values_fill = 0) %>% 
  dplyr::mutate(across(-c(STRING_id,total), ~./total))

#Add to STRING data and create dummy group for genes without enrichment
map.unique <- map %>% 
  #collapse gene names
  group_by(STRING_id) %>% 
  dplyr::mutate(gene = paste(unique(gene), collapse=" / ")) %>% 
  #add enrichment color info
  full_join(col.mat.format, by = "STRING_id") %>%
  distinct() %>% 
  #no enrichment group
  dplyr::mutate(none = ifelse(is.na(total),1,0)) %>% 
  ungroup() %>% 
  #fill in 0
  dplyr::mutate(across(everything(), ~replace_na(., 0)))

# Create igraph object 
subgraph <- string_db$get_subnetwork(map.unique$STRING_id)

# Arrange metadata as in network
map.arrange <- map.unique %>% 
  dplyr::filter(STRING_id %in% vertex_attr(subgraph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(subgraph)$name)))

#Check order
identical(vertex_attr(subgraph)$name, map.arrange$STRING_id)
##gene names
V(subgraph)$symbol <- map.arrange$gene
##enrichment colors
for(term in colnames(map.arrange)[-c(1:3)]){
  vertex_attr(subgraph)[[term]] <- unlist(map.arrange[term])
}

#colors 
color.vec <- c(db.GOBP.format2$color,'#d9d9d9')
names(color.vec) <- c(as.character(db.GOBP.format2$Description),"none")


set.seed(8434)
##set layout
xy <- layout_with_lgl(subgraph) 

V(subgraph)$x <- xy[, 1]
V(subgraph)$y <- xy[, 2]

plot <- ggraph(subgraph, layout= "manual", 
               x = V(subgraph)$x, y = V(subgraph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey70") +
  scale_edge_width_continuous(range = c(0.2,2), name="STRING score")

#Add nodes
plot.col <- plot + 
    geom_scatterpie(data=as_data_frame(subgraph, "vertices"),
                    cols=colnames(map.arrange)[-c(1:3)], color=NA,
                    pie_scale = 0.8) +
    scale_fill_manual(values=color.vec, name="Gene ontology") +
    geom_nodetext(aes(x = V(subgraph)$x, y = V(subgraph)$y,
                      label=V(subgraph)$symbol), size=4) +
    theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

plot.col

##### Save #####
ggsave("../publication/STRING/Fig5.U.STRING.lrg.pdf", 
       plot.col, 
       height=15, width=8)

###########################################################################
##### DEGs small STRING clusters #####
#DEGs in cluster
med <- c("YTHDC2","MYBBP1A","DDX10","URB1","NOP56",
          "WDR3","RRP15","CIRH1A","DDX51","NOL8")
sm1<- c("CUL3","RNF6","FBXO4","NEDD4","ZNRF2")
sm2 <- c("HS3ST1","GCNT1","ST3GAL1","ST3GAL2","B4GALT1")
sm3 <- c("USB1","SYF2","SF3A3","HNRNPM","CD2BP2")
sm4 <- c("CEP164","TUBGCP6","MARK4")

all.small <- c(med,sm1,sm2,sm3,sm4)

##### Format color db with selected terms #####
#Custom groups
GOID.signif <- list()
#med and sm3
GOID.signif[["rna"]] <- c("GO:0022613",
                          #"GO:0042254", #Collapse with one above
                          "GO:0006396",
                          #"GO:0034470",
                          #"GO:0034660",
                          "GO:0016072",
                          "GO:0006364",
                          "GO:0008380",
                          "GO:0006397"
                          #"GO:0000375",
                          #"GO:0000377",
                          #"GO:0000398"
                          )
#sm1
GOID.signif[["ubiq"]] <- c("GO:0043632",
                           "GO:0019941",
                           "GO:0006511",
                           "GO:0070647",
                           "GO:0032446",
                           "GO:0016567",
                           "GO:0000209")
#sm2
GOID.signif[["glyco"]] <- c("GO:1903510",
                            #"GO:0006024",
                            "GO:0042339",
                            "GO:0018146",
                            "GO:0009101"
                            #"GO:0006486",
                            #"GO:0006493",
                            #"GO:0016266"
                            )

#Make topGO object
genes.vec <- U %>% 
  mutate(group = ifelse(gene %in% all.small, 1, 0),
         group = factor(group)) %>% 
  dplyr::select(group) %>% unlist(use.names=FALSE)
names(genes.vec) <- U$gene

GOdata <- new("topGOdata", ontology = "BP",
              allGenes = genes.vec,
              annot = annFUN.org, mapping="org.Hs.eg.db", ID = "symbol",
              nodeSize = 1)

#Get genes in terms of interest
db.GOBP <- plyr::ldply(genesInTerm(GOdata), data.frame, .id = "GOID") %>% 
  dplyr::rename("gene"="X..i..")

# size terms
db.size <- db.GOBP %>% 
  distinct(GOID,gene) %>% 
  count(GOID)

#Rank order all term groups
GOID.signif.df <- plyr::ldply(GOID.signif, data.frame, .id="group") %>% 
  rename("X..i.."="GOID") %>% 
  inner_join(db.size) %>% 
  arrange(-n) %>% 
  group_by(group) %>% 
  mutate(ID = row_number())

#Get most specific term per group for each gene
top.terms <- db.GOBP  %>% 
  #signif terms
  right_join(GOID.signif.df) %>%
  #DEGs
  filter(gene %in% all.small) %>%
  #Get most specific term per gene
  group_by(group, gene) %>%
  summarise(ID = max(ID))

# get GOID Descriptions
GO.term <- as.data.frame(Term(GOTERM)) %>% 
  rownames_to_column("GOID") %>% 
  dplyr::rename("Description"="Term(GOTERM)") %>% 
  mutate(Description = gsub("_|-|[/]"," ",Description))

#Put it all together and force into format for fxn
db.GOBP.format <- db.GOBP %>% 
  #signif terms
  right_join(GOID.signif.df) %>%
  #DEGs
  filter(gene %in% all.small)%>% 
  #add desc
  inner_join(top.terms) %>% 
  inner_join(GO.term) %>% 
  #format for fxn
  group_by(group, ID, GOID, Description) %>% 
  summarise(SYMBOLs = paste(unique(gene), collapse="/"),
            size.overlap.term=n()) %>%
  ungroup() %>% 
  mutate(p.adjust=0) %>%
  mutate(Description = factor(Description)) %>%
  arrange(Description) %>% 
  #order
  arrange(group, ID)

#Factor order within group
db.GOBP.format2 <- db.GOBP.format %>% 
  group_by(ID, GOID, Description, p.adjust) %>% 
  summarise(#group = group[1],
            SYMBOLs = paste0(SYMBOLs),
            size.overlap.term = sum(size.overlap.term)) %>% 
  ungroup() %>% 
  mutate(Description = factor(Description, levels=unique(db.GOBP.format$Description))) %>% 
  arrange(Description) %>% 
  mutate(color = c('#cbc9e2','#9e9ac8','#6a51a3', #rna purple
                   '#e31a1c', #uniq red
                    '#dfc27d','#bf812d' #glyco brown
  ))
    
#### Plot SMALL ####
#Format gene vector to matrix
genes.mat <- as.matrix(all.small)
colnames(genes.mat) <- "gene"

#Map genes to STRING
map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)

#Get enrichments
col.mat <- db.GOBP.format2 %>% 
  dplyr::select(Description, size.overlap.term, p.adjust,SYMBOLs)

#Format enrichment results for scatterpie plotting
col.mat.format <- col.mat %>% 
  #Split gene lists within terms
  dplyr::mutate(gene = strsplit(as.character(SYMBOLs), "/")) %>% 
  tidyr::unnest(gene) %>%
  #spread terms to columns
  distinct(gene, Description) %>% 
  dplyr::mutate(value=1) %>% 
  #Add string ID
  left_join(map, by = "gene") %>% 
  dplyr::select(-gene) %>% 
  dplyr::distinct() %>% 
  #Calculate terms per ID
  group_by(STRING_id) %>% 
  dplyr::mutate(total = sum(value)) %>%
  ungroup() %>% 
  #Calculate proportions within terms
  arrange(match(Description,unique(db.GOBP.format2$Description))) %>% 
  pivot_wider(names_from = Description, values_fill = 0) %>% 
  dplyr::mutate(across(-c(STRING_id,total), ~./total))

#Add to STRING data and create dummy group for genes without enrichment
map.unique <- map %>% 
  #collapse gene names
  group_by(STRING_id) %>% 
  dplyr::mutate(gene = paste(unique(gene), collapse=" / ")) %>% 
  #add enrichment color info
  full_join(col.mat.format, by = "STRING_id") %>%
  distinct() %>% 
  #no enrichment group
  dplyr::mutate(none = ifelse(is.na(total),1,0)) %>% 
  ungroup() %>% 
  #fill in 0
  dplyr::mutate(across(everything(), ~replace_na(., 0)))

# Create igraph object 
subgraph <- string_db$get_subnetwork(map.unique$STRING_id)

# Arrange metadata as in network
map.arrange <- map.unique %>% 
  dplyr::filter(STRING_id %in% vertex_attr(subgraph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(subgraph)$name)))

#Check order
identical(vertex_attr(subgraph)$name, map.arrange$STRING_id)
##gene names
V(subgraph)$symbol <- map.arrange$gene
##enrichment colors
for(term in colnames(map.arrange)[-c(1:3)]){
  vertex_attr(subgraph)[[term]] <- unlist(map.arrange[term])
}

#colors 
color.vec <- c(db.GOBP.format2$color,'#d9d9d9')
names(color.vec) <- c(as.character(db.GOBP.format2$Description),"none")


set.seed(8434)
##set layout
xy <- layout_with_fr(subgraph) 

V(subgraph)$x <- xy[, 1]
V(subgraph)$y <- xy[, 2]

plot2 <- ggraph(subgraph, layout= "manual", 
               x = V(subgraph)$x, y = V(subgraph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey70") +
  scale_edge_width_continuous(range = c(0.2,2), name="STRING score")

#Add nodes
plot.col2 <- plot2 + 
  geom_scatterpie(data=as_data_frame(subgraph, "vertices"),
                  cols=colnames(map.arrange)[-c(1:3)], color=NA,
                  pie_scale = 0.8) +
  scale_fill_manual(values=color.vec, name="Gene ontology") +
  geom_nodetext(aes(x = V(subgraph)$x, y = V(subgraph)$y,
                    label=V(subgraph)$symbol), size=3) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")

plot.col2

##### Save #####
ggsave("../publication/STRING/Fig5.U.STRING.sml.pdf", 
       plot.col2, 
       height=12, width=8)
