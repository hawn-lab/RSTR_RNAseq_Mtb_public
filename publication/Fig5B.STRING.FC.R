#Working with dataframes
require(tidyverse)
#STRING database
require(STRINGdb)
#Working with network objects
require(igraph) #vertex_attr()
#Graphing networks
require(ggraph)
#require(scatterpie) #geom_scatterpie()
require(ggnetwork) #theme_blank()
require(scales)
set.seed(8434)

#### Load data ####
#List DEGs in cluster
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

med <- c("YTHDC2","MYBBP1A","DDX10","URB1","NOP56",
         "WDR3","RRP15","CIRH1A","DDX51","NOL8")
sm1<- c("CUL3","RNF6","FBXO4","NEDD4","ZNRF2")
sm2 <- c("HS3ST1","GCNT1","ST3GAL1","ST3GAL2","B4GALT1")
sm3 <- c("USB1","SYF2","SF3A3","HNRNPM","CD2BP2")
sm4 <- c("CEP164","TUBGCP6","MARK4")

all.small <- c(med,sm1,sm2,sm3,sm4)

##### Fold change color #####
attach("../Uganda/data_clean/RSTR_RNAseq_data_combined_uniqueFULLID.RData")

log2FC.TB <- as.data.frame(dat.combined.voom$E) %>% 
  #DEGs
  rownames_to_column("geneName") %>% 
  filter(geneName %in% c(lrg, all.small)) %>% 
  #Add metadata
  pivot_longer(-geneName, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  #Calculate mean TB-media
  dplyr::select(geneName, value, RSID, Sample_Group, condition) %>% 
  pivot_wider(names_from = condition) %>% 
  mutate(delta=TB-MEDIA) %>% 
  #average delta within RSTR and LTBI
  group_by(geneName, Sample_Group) %>% 
  summarise(aveDelta = mean(delta, na.rm=TRUE)) %>% 
  #Calculate delta delta
  pivot_wider(names_from = Sample_Group, values_from = aveDelta) %>% 
  mutate(deltaDelta = RSTR-LTBI) %>% 
  #fix ORF gene name
  mutate(geneName = gsub("orf","ORF",geneName))

log2FC.RSTR <- as.data.frame(dat.combined.voom$E) %>% 
  #DEGs
  rownames_to_column("geneName") %>% 
  filter(geneName %in% c(lrg, all.small)) %>% 
  #Add metadata
  pivot_longer(-geneName, names_to = "libID") %>% 
  left_join(dat.combined.voom$targets) %>% 
  #Calculate mean RSTR-LTBI
  group_by(geneName, Sample_Group, condition) %>% 
  summarise(mean = mean(value, na.rm=TRUE)) %>% 
  pivot_wider(names_from = Sample_Group, values_from = mean) %>% 
  mutate(delta=RSTR-LTBI)%>% 
  #average delta within RSTR and LTBI
  group_by(geneName, condition) %>% 
  summarise(aveDelta = mean(delta, na.rm=TRUE)) %>% 
  #Calculate delta delta
  pivot_wider(names_from = condition, values_from = aveDelta) %>% 
  #fix ORF gene name
  mutate(geneName = gsub("orf","ORF",geneName))

#Combine
log2FC <- full_join(log2FC.TB,log2FC.RSTR)

#### STRING DB ####
#STRING database
string_db <- STRINGdb$new(version="11", species=9606,
                          score_threshold=700, input_directory="")

#### Large cluster ####
#Format gene vector to matrix
genes.mat <- as.matrix(lrg)
colnames(genes.mat) <- "gene"

#Map genes to STRING
map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)

#Collapse duplicate STRING ID
map.unique <- map %>% 
  group_by(STRING_id) %>% 
  dplyr::summarise(gene = paste(unique(gene), collapse = " / "), .groups="drop") %>% 
  #add fold change
  left_join(log2FC, by=c("gene"="geneName"))

#### Network ####
# Create igraph object 
subgraph <- string_db$get_subnetwork(map.unique$STRING_id)

# Arrange metadata as in network
map.arrange <- map.unique %>% 
  dplyr::filter(STRING_id %in% vertex_attr(subgraph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(subgraph)$name)))

# Set attributes
## Check order first
identical(vertex_attr(subgraph)$name, map.arrange$STRING_id)

##gene names
V(subgraph)$symbol <- map.arrange$gene
##FC colors
# V(subgraph)$FC <- map.arrange$deltaDelta
# V(subgraph)$rstrFC <- map.arrange$RSTR
# V(subgraph)$ltbiFC <- map.arrange$LTBI
# V(subgraph)$mediaFC <- map.arrange$MEDIA
V(subgraph)$mtbFC <- map.arrange$TB

#### Plot ####
#Get xy of nodes for manual layout
set.seed(8434)
xy <- layout_with_lgl(subgraph) 

V(subgraph)$x <- xy[, 1]
V(subgraph)$y <- xy[, 2]

plot <- ggraph(subgraph, layout= "manual", 
               x = V(subgraph)$x, y = V(subgraph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey70") +
  scale_edge_width(range = c(0.2,2), name="STRING score") 

plot.col <- plot + 
  geom_node_point(data=as_data_frame(subgraph, "vertices"),
                  aes(color=mtbFC), size=7) +
  scale_color_gradient2(low="darkblue", mid="white", high="darkred", 
                        name="log2 CPM\nRSTR - LTBI in +Mtb", limits=c(-3,3)) +
  geom_nodetext(aes(x = V(subgraph)$x, y = V(subgraph)$y,
                    label=V(subgraph)$symbol), size=3) +
  theme_blank() + coord_fixed() +
  theme(legend.position = "bottom", legend.direction = "vertical")
plot.col

ggsave("../publication/STRING/Fig5.U.STRING.lrg.FC.pdf", 
       plot.col, 
       height=15, width=8)

###################################################################
#### Small clusters ####
#Format gene vector to matrix
genes.mat <- as.matrix(all.small)
colnames(genes.mat) <- "gene"

#Map genes to STRING
map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)

#Collapse duplicate STRING ID
map.unique <- map %>% 
  group_by(STRING_id) %>% 
  dplyr::summarise(gene = paste(unique(gene), collapse = " / "), .groups="drop") %>% 
  #add fold change
  left_join(log2FC, by=c("gene"="geneName"))

#### Network ####
# Create igraph object 
subgraph <- string_db$get_subnetwork(map.unique$STRING_id)

# Arrange metadata as in network
map.arrange <- map.unique %>% 
  dplyr::filter(STRING_id %in% vertex_attr(subgraph)$name) %>% 
  arrange(match(STRING_id, c(vertex_attr(subgraph)$name)))

# Set attributes
## Check order first
identical(vertex_attr(subgraph)$name, map.arrange$STRING_id)

##gene names
V(subgraph)$symbol <- map.arrange$gene
##FC colors
# V(subgraph)$FC <- map.arrange$deltaDelta
# V(subgraph)$rstrFC <- map.arrange$RSTR
# V(subgraph)$ltbiFC <- map.arrange$LTBI
# V(subgraph)$mediaFC <- map.arrange$MEDIA
V(subgraph)$mtbFC <- map.arrange$TB

#### Plot ####
#Get xy of nodes for manual layout
set.seed(8434)
xy <- layout_with_fr(subgraph) 

V(subgraph)$x <- xy[, 1]
V(subgraph)$y <- xy[, 2]

plot2 <- ggraph(subgraph, layout= "manual", 
               x = V(subgraph)$x, y = V(subgraph)$y) +
  #Edges
  geom_edge_link(aes(width=combined_score), color="grey70") +
  scale_edge_width(range = c(0.2,2), name="STRING score") 

plot.col2 <- plot2 + 
  geom_node_point(data=as_data_frame(subgraph, "vertices"),
                  aes(color=mtbFC), size=7) +
  scale_color_gradient2(low="darkblue", mid="white", high="darkred", 
                        name="log2 CPM\nRSTR - LTBI in +Mtb", limits=c(-3,3)) +
  geom_nodetext(aes(x = V(subgraph)$x, y = V(subgraph)$y,
                    label=V(subgraph)$symbol), size=3) +
  theme_blank() + coord_fixed()+
  theme(legend.position = "bottom", legend.direction = "vertical")
plot.col2

ggsave("../publication/STRING/Fig5.U.STRING.sml.FC.pdf", 
       plot.col2, 
       height=12, width=8)
