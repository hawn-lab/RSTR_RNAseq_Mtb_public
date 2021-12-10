library(tidyverse)
library(cowplot)
library(venn)

##### UGANDA #####
#### Load data ####
#Model results for various co-variates
models <- list.files(path="../Uganda/results/model_selection/",
                     pattern=".csv.gz", full.names = TRUE)

result <- data.frame()
for(file in models){
  result <- bind_rows(result, read_csv(file))
}

#### Format data ####

result.format <- result %>% 
  #Keep interaction term results
  filter(grepl(":", variable)) %>% 
  #Spread FDR and sigma values by model
  distinct(group, model, gene, FDR, sigma) %>% 
  pivot_longer(FDR:sigma) %>% 
  mutate(name = paste(name, model, sep="_")) %>% 
  select(-model) %>% 
  pivot_wider() %>% 
  #Create groups for FDR significance
  mutate(fdr.group = ifelse(FDR_lme <= 0.2 & FDR_lmekin <= 0.2, "both",
                            ifelse(FDR_lme <= 0.2 & FDR_lmekin > 0.2, "lme",
                                   ifelse(FDR_lme > 0.2 & FDR_lmekin <= 0.2, "lmekin",
                                          "none")))) %>% 
  #Create groups for best fit model significance
  mutate(fit.group = ifelse(sigma_lme < sigma_lmekin, "+ age + sex + batch", 
                            "+ age + sex + batch + kinship")) 

#### Plot 1: kinship ####
# Model fits (sigma) colored by significance (FDR) without vs with kinship

plot1 <- result.format %>% 
  filter(group == "RSTR.interaction_age.sex.batch") %>% 
  arrange(fit.group) %>% 

  ggplot(aes(x=sigma_lmekin, y=sigma_lme, color=fit.group)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "+ age + sex + batch + kinship\nModel fit (sigma)",
       y = "Model fit (sigma)\n+ age + sex + batch",
       color = "Best fit", title="Kinship") +
  coord_fixed() +
  scale_color_manual(values=c("#d95f02", "#1b9e77"),
                     guide = guide_legend(reverse = TRUE))

#plot1

#### Plot 2: co-variates ####
#Format data for co-variate comparisons
result.format.covar <- result.format %>% 
  #Make best fit group by covariate
  select(group, gene, sigma_lmekin) %>% 
  pivot_wider(names_from = group, values_from = sigma_lmekin) %>% 
  # +/- age
  mutate(fit.group2 = ifelse(RSTR.interaction_sex.batch <= RSTR.interaction_age.sex.batch,
                             "+ sex + batch + kinship", 
                             "+ age + sex + batch + kinship")) %>% 
  # +/- sex
  mutate(fit.group3 = ifelse(RSTR.interaction_batch <= RSTR.interaction_age.sex.batch,
                             "+ batch + kinship", "+ age + sex + batch + kinship")) 

# +/- age
plot2 <- result.format.covar %>% 
  arrange(desc(fit.group2)) %>% 
  
  ggplot(aes(x=RSTR.interaction_age.sex.batch, 
             y=RSTR.interaction_sex.batch, 
             color=fit.group2)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "+ age + sex + batch + kinship\nModel fit (sigma)",
       y = "Model fit (sigma)\n+ sex + batch + kinship",
       color = "Best fit", title="Age") +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  coord_fixed()

#plot2

# +/- sex
plot3 <- result.format.covar %>% 
  arrange(desc(fit.group3)) %>% 
  
  ggplot(aes(x=RSTR.interaction_age.sex.batch, 
             y=RSTR.interaction_batch, 
             color=fit.group3)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  labs(x = "+ age + sex + batch + kinship\nModel fit (sigma)",
       y = "Model fit (sigma)\n+ batch + kinship",
       color = "Best fit", title="Sex") +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  coord_fixed()

#plot3

#### Save ####

ggsave("FigS1A.model.sigma.pdf", 
       plot_grid(plot1,plot2,plot3, 
                 labels = c("a)","",""), label_fontface = "plain",
                 nrow=1, align="hv", axis = "bt"),
       height=4, width=8)

tiff("FigS1A.model.sigma.tiff",  height=4, width=8, units = "in", res=300)
plot_grid(plot1,plot2,plot3, 
          labels = c("a)","",""), label_fontface = "plain",
          nrow=1, align="hv", axis = "bt")
dev.off()

#### Plot C: DEG venn ####
#Venns of differentially expressed genes
pdf("FigS1C.DEG.venn.pdf", width=4.5, height=4)
par(mfrow=c(1,1))
venn.ls <- list() 

venn.ls[["+ age + sex\n+ batch + kinship"]] <- result.format %>% 
    filter(group == "RSTR.interaction_age.sex.batch" & FDR_lmekin <= 0.2) %>% 
    distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["+ sex + batch + kinship"]] <- result.format %>% 
  filter(group == "RSTR.interaction_sex.batch" & FDR_lmekin <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["+ batch + kinship"]] <- result.format %>% 
  filter(group == "RSTR.interaction_batch" & FDR_lmekin <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
  
venn(ilab=FALSE, ilcs=1, sncs=1,
               x=venn.ls, box=FALSE)
  title(sub="Uganda total Mtb:RSTR DEGs", line = -1)
  title( adj=0, main="c)", font.main=1, cex=1.1, col="black", line=-1)
dev.off()

tiff("FigS1C.DEG.venn.tiff", width=4.5, height=4, units = "in", res=300)
par(mfrow=c(1,1))
venn.ls <- list() 

venn.ls[["+ age + sex\n+ batch + kinship"]] <- result.format %>% 
  filter(group == "RSTR.interaction_age.sex.batch" & FDR_lmekin <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["+ sex + batch + kinship"]] <- result.format %>% 
  filter(group == "RSTR.interaction_sex.batch" & FDR_lmekin <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)
venn.ls[["+ batch + kinship"]] <- result.format %>% 
  filter(group == "RSTR.interaction_batch" & FDR_lmekin <= 0.2) %>% 
  distinct(gene) %>% unlist(use.names = FALSE)

venn(ilab=FALSE, ilcs=1, sncs=1,
     x=venn.ls, box=FALSE)
title(sub="Uganda total Mtb:RSTR DEGs", line = -1)
title( adj=0, main="c)", font.main=1, cex=1.1, col="black", line=-1)
dev.off()

#### Best fit totals ####
result.format %>% 
  filter(group == "RSTR.interaction_age.sex.batch") %>% 
  count(fit.group)

result.format.covar %>% 
  count(fit.group2)

result.format.covar %>% 
  count(fit.group3)

##### SOUTH AFRICA #####
#### Load data ####
#Model results for various co-variates
models <- list.files(path="../SAfrica//results/model_selection/",
                     pattern=".csv.gz", full.names = TRUE)

result <- data.frame()
for(file in models){
  result <- bind_rows(result, read_csv(file))
}

#### Format data ####

result.format <- result %>% 
  #Keep interaction term results
  filter(grepl(":", variable)) %>% 
  mutate(group= recode_factor(factor(group), "SA_RSTR.interaction"="none",
                "SA_RSTR.interaction_age"="age")) %>% 
  #Spread FDR and sigma values by model
  distinct(group, gene, FDR, sigma) %>% 
  pivot_longer(FDR:sigma) %>% 
  mutate(name = paste(name, group, sep="_")) %>% 
  select(-group) %>% 
  pivot_wider() %>% 
  mutate(best = ifelse(sigma_none<sigma_age, "none", 
                       ifelse(sigma_age<sigma_none, "+ age", NA)))

#### Plot 4: co-variates ####
# +/- age
plot4 <- result.format %>% 
  arrange(desc(best)) %>% 
  
  ggplot(aes(x=sigma_age, 
             y=sigma_none, 
             color=best)) +
  geom_point(alpha=0.2) +
  geom_abline(slope=1, intercept=0) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "+ age\nModel fit (sigma)",
       y = "Model fit (sigma)\nnone",
       color = "Best fit", title="Age") +
  scale_color_manual(values=c("#1b9e77","#d95f02")) +
  coord_fixed()

plot4

#### Save ####

ggsave("FigS1B.model.sigma.pdf", 
       plot_grid(plot4, 
                 labels = c("b)","",""), label_fontface = "plain",
                 nrow=1, align="hv", axis = "bt"),
       height=4, width=3.2)
tiff("FigS1B.model.sigma.tiff",  height=4, width=3.2, units = "in", res=300)
plot_grid(plot4, 
          labels = c("b)","",""), label_fontface = "plain",
          nrow=1, align="hv", axis = "bt")
dev.off()

#### Best fit totals ####
result.format %>% 
  count(best)
