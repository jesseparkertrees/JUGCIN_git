##################################################################
#######  BUTTERNUT SPATIAL ANALYSES  #############################
##################################################################
# LOAD PACKAGES
library(vegan)
library(usdm)
library(ggplot2)
library(grid)
library(maps)
library(dplyr)
library(terra)
library(tidyr)
library(geodata)
##################################################################
## LOAD DATA AND FILTER
pops_health6<-read.csv("pops_health6.csv")
pops_health6[is.na(pops_health6$JC_struc), ] # should be 0
summary(as.factor(pops_health6$NAME_1))
pops_health6[pops_health6$NAME_1=='',]
pops_health6[pops_health6$NAME_2=='',]
pops_health6[pops_health6$HASC_2=='',]
colnames(pops_health6)
pops_health6$X<-NULL
##################################################################
# SUMMARIZE NEWHYBIDS CATEGORIES
## filter out trees that there were not NewHybrids assignments for
dim(pops_health7)
pops_health7<-pops_health6[pops_health6$NewHyb_FinalAssignment!="",]
sum(pops_health6$NewHyb_FinalAssignment=="")
dim(pops_health6)
dim(pops_health7)
## summarize categories
summary(as.factor(pops_health7$NewHyb_FinalAssignment))
histplot<-ggplot(pops_health7, aes(x = NewHyb_FinalAssignment)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Hybrid Category",
    y = "Count")
histplot
newhyb_table <- pops_health7 %>%
  count(NewHyb_FinalAssignment)
write.csv(newhyb_table,"./outputs/summaries/NewHybridsCategoriesSummaryTableALL.csv", row.names = FALSE)

## summarize categories by admin unit
summary_table <- pops_health7 %>%
  count(HASC_2, NewHyb_FinalAssignment) %>%
  pivot_wider(
    names_from = NewHyb_FinalAssignment,
    values_from = n,
    values_fill = 0
  )

summary_table
summary_table <- pops_health7 %>%
  count(HASC_2, NewHyb_FinalAssignment) %>%
  pivot_wider(
    names_from = NewHyb_FinalAssignment,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(
    pops_health7 %>%
      group_by(HASC_2) %>%
      summarize(mean_JC_struc = mean(JC_struc, na.rm = TRUE)),
    by = "HASC_2"
  )
summary_table <- summary_table %>% mutate(Hybrids = rowSums(across(-c(HASC_2, JC, JA, mean_JC_struc))))
summary_table <- summary_table %>% mutate(Hybrid_perc = Hybrids/rowSums(across(-c(HASC_2, mean_JC_struc, Hybrids))))
write.csv(summary_table,"./outputs/summaries/NewHybridsCategoriesSummaryTableCounties.csv", row.names = FALSE)
summary_table2 <- pops_health7 %>%
  count(NAME_1, NewHyb_FinalAssignment) %>%
  pivot_wider(
    names_from = NewHyb_FinalAssignment,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(
    pops_health7 %>%
      group_by(NAME_1) %>%
      summarize(mean_JC_struc = mean(JC_struc, na.rm = TRUE)),
    by = "NAME_1"
  )
summary_table2 <- summary_table2 %>% mutate(Hybrids = rowSums(across(-c(NAME_1, JC, JA, mean_JC_struc))))
summary_table2 <- summary_table2 %>% mutate(Hybrid_perc = Hybrids/rowSums(across(-c(NAME_1, mean_JC_struc, Hybrids))))
summary_table2
write.csv(summary_table2,"./outputs/summaries/NewHybridsCategoriesSummaryTableStates.csv", row.names = FALSE)

## plot on map
usa <- gadm(country = "USA", level = 2, path = "../JUGCIN_git_externalFiles")
can <- gadm(country = "CAN", level = 2, path = "../JUGCIN_git_externalFiles")
## combine USA and CAN
admin <- rbind(usa, can)
unique(admin$NAME_1)
unique(summary_table2$NAME_1)
states_df <- as.data.frame(admin)
states_df <- states_df %>%
  left_join(summary_table2, by = "NAME_1")
states <- merge(admin, states_df[, c("NAME_1", "Hybrids", "Hybrid_perc", "mean_JC_struc")],
                by = "NAME_1", all.x = TRUE)
library(terra)
states_simple <- simplifyGeom(states, tolerance = 0.01)
plot(states_simple, "Hybrids")
library(sf)
library(ggplot2)
states_sf <- st_as_sf(states_simple)
ggplot(states_sf) +
  geom_sf(aes(fill = Hybrid_perc), color = "grey40", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Hybrid\nproportion", na.value = "white") +
  theme_void()
ggplot(states_sf) +
  geom_sf(aes(fill = Hybrids), color = "grey40", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Hybrids", na.value = "white") +
  theme_void()

##################################################################
# spatial analysis
pops_health5<-read.csv("pops_health5_04272026.csv")
pops_health5$X<-NULL
library(maps)
# Extract state (returns "state_name:region")
states <- map.where(database = "state", pops_health5$x_new, pops_health5$y_new)
states <- sapply(strsplit(states, ":"), "[", 1)
print(states) 
pops_health5$state<-states

pops_health5$HybridYN<-ifelse(pops_health5$NH_final_assignment=="NewHyb_Butternut" | pops_health5$NH_final_assignment=="NewHyb_Heartnut", 0,1 )
summary(as.factor(pops_health5$NH_final_assignment))

summary(pops_health5$CLUSTER_ID)
pops_health5$CLUSTER_ID<-as.character(pops_health5$CLUSTER_ID)

sum_clust<-pops_health5 %>%
  group_by(CLUSTER_ID, NH_final_assignment, state) %>%
  summarise(count= n(), .groups="drop")

sum_state1 <- pops_health5 %>%
  group_by(state, NH_final_assignment) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = NH_final_assignment,
    values_from = n,
    values_fill = 0  )
colnames(sum_state1)
sum_state1$hyb<-rowSums(sum_state1[,c(3:9,11:14)])
sum_state1$JC<-rowSums(sum_state1[,2])
sum_state1$JA<-rowSums(sum_state1[,10])
sum_state1$tot<-rowSums(sum_state1[,2:14])


sum_clust2 <- sum_clust %>%
  tidyr::pivot_wider(
    names_from = NH_final_assignment,
    values_from = count,
    values_fill = 0  )
colnames(sum_clust2)
sum_clust2$hyb<-rowSums(sum_clust2[,c(2:6,8:9,11:14)])
sum_clust2$JC<-rowSums(sum_clust2[,7])
sum_clust2$JA<-rowSums(sum_clust2[,10])
sum_clust2$tot<-rowSums(sum_clust2[,2:14])

dim(sum_clust2[sum_clust2$CLUSTER_ID!=0 & sum_clust2$hyb>0,])
(58+45)/353

dim(sum_clust2[sum_clust2$CLUSTER_ID!=0 & sum_clust2$JA>0,])
7/353

length(sum_clust2$hyb[sum_clust2$hyb==0])
length(sum_clust2$hyb[sum_clust2$hyb>0])
length(unique(sum_clust2$CLUSTER_ID))
58/159

ILM<-pops_health5[pops_health5$CLUSTER_ID=="3",]
ILM$height_bin <- cut(
  ILM$PlantHeight_ft,
  breaks = c(-Inf, 5, 15, 30, Inf),
  labels = c("<5 ft", "5–15 ft", "15–30 ft", ">30 ft"),
  right = FALSE
)

ILM2 <- ILM[!is.na(ILM$height_bin), ]

length(ILM2$height_bin)
length(ILM2$NH_final_assignment)
height_hyb_table <- table(
  ILM2$height_bin,
  ILM2$NH_final_assignment) %>% as.data.frame() %>%
  pivot_wider(names_from=Var2, values_from = Freq)
#write.csv(height_hyb_table,"height_hyb_table.csv")
summary(as.factor(pops_health5$CLUSTER_ID))
CP<-pops_health5[pops_health5$CLUSTER_ID=="7",]
CPfil<-CP[CP$PlantHeight_ft>0,]
plot(CP$x_new, CP$y_new)
plot(CPfil$x_new, CPfil$y_new)

CP$height_bin <- cut(
  ILM$PlantHeight_ft,
  breaks = c(-Inf, 5, 15, 30, Inf),
  labels = c("<5 ft", "5–15 ft", "15–30 ft", ">30 ft"),
  right = FALSE
)


#check correlation
library(tidyr)
pops_health5_var<-pops_health5 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_1,wc2.1_30s_bio_11,wc2.1_30s_bio_12,wc2.1_30s_bio_3) %>% drop_na() %>% as.data.frame()
cor(pops_health5_var)
library(usdm)
usdm::vif(pops_health5_var)
pops_health5_var<-pops_health5 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_12, wc2.1_30s_bio_15) %>% drop_na() %>% as.data.frame()
usdm::vif(pops_health5_var)

#rescale data
pops_health5 <- pops_health5 |>
  dplyr::mutate(
    PlantHeight_ft = as.numeric(PlantHeight_ft),
    human_footprint_z = scale(wildareas.v3.2009.human.footprint),
    ForestEdge_30m_z = scale(ForestEdge_30m),
    ForestEdge_NorAmer_custom_z = scale(ForestEdge_NorAmer_custom),
    bio6_z = scale(wc2.1_30s_bio_6),
    bio12_z = scale(wc2.1_30s_bio_12),
    bio15_z = scale(wc2.1_30s_bio_15),
    ph_z = scale(phh2o_0.5cm),
    nit_z = scale(nitrogen_0.5cm),
    CLUSTER_ID=as.factor(CLUSTER_ID),
    NLCD_1km=as.factor(NLCD_1km),
    NLCD_5km=as.factor(NLCD_5km),
    NLCD=as.factor(NLCD),
    TCC_z=scale(TCC),
    TCC_1km_z=scale(TCC_1km),
    TCC_5km_z=scale(TCC_5km)
    )

##################################################################
# spatial analysis
colnames(pops_health5)

library("glmmTMB")
library(performance)
library("ggeffects")

mod_beta <- glmmTMB(HybridYN ~ ForestEdge_30m_z*TCC_5km_z*bio12_z*nit_z+(1|CLUSTER_ID),family = binomial(),data = pops_health5)

mod_beta <- glmmTMB(HybridYN ~ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health5)

mod_beta <- glmmTMB(HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health5)

summary(mod_beta)
r2(mod_beta)
plot(ggpredict(mod_beta, terms = "TCC_5km_z"))
plot(ggpredict(mod_beta, terms = "human_footprint_z"))
plot(ggpredict(mod_beta, terms = "NLCD_5km"))
plot(ggpredict(mod_beta, terms = c("ForestEdge_30m_z","bio6_z")))
plot(ggpredict(mod_beta, terms = "bio12_z"))
plot(ggpredict(mod_beta, terms = "bio15_z"))

par(mar=c(8,8,8,8))
boxplot(pops_health5$human_footprint_z~pops_health5$NH_final_assignment, las=2)
boxplot(pops_health5$TCC_1km~pops_health5$NH_final_assignment, las=2)

library("DHARMa")
sim <- simulateResiduals(mod_beta)
plot(sim)
testDispersion(sim) # no dispersion problem

#removed ph_Z, bio6_z, human_footprint_z, ForestEdge_NorAmer_custom_z (not significant) and bio_12_z (sig. effect with TCC_5km and nit_z)
# odds of being a hybrid increases with increasing nitrogen, annual precip, and on/near developed land () and with decreasing distance from forest edge (within forest), and decreasing canopy cover (averaged over 5km)
# R2 for Mixed Models: Conditional R2: 0.468/ Marginal R2: 0.124


# remove sites with no hybrids or JA


# health analysis
health_vars <- pops_health5 %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))
colnames(pops_health5)
#ancestry <- select(pops_health5, "HybridYN","CLUSTER_ID")
ancestry <- select(pops_health5, "HybridYN")
ancestry$HybridYN<-as.factor(ancestry$HybridYN)
#ancestry$CLUSTER_ID<-as.factor(ancestry$CLUSTER_ID)
#ancestry <- health_merge$NH_final_assignment %>% as.factor()
pops_health5$TCC_5km_z
envdata<-select(pops_health5,"bio6_z","ForestEdge_30m_z","ph_z","PlantHeight_ft", "bio15_z", "TCC_z","TCC_5km_z")
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-pops_health5[complete_rows,]
dim(ancestry_clean)
health_vars_scaled <- scale(health_vars_clean)
envdata_scaled<-scale(envdata_clean[,]) %>% as.data.frame
colnames(envdata_scaled)
predictors<-cbind(envdata_scaled,ancestry_clean)
library(vegan)
hyb_rda <- rda(health_vars_scaled ~ ., data = predictors)
summary(hyb_rda)
RsquareAdj(hyb_rda)
scores(hyb_rda, display = "species", scaling = 2)
anova(hyb_rda)
anova(hyb_rda, by = "term")

ind_scores <- scores(hyb_rda, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
#ind_scores_df$ID <- edge@ind.names
#ind_scores_df$pop <- edge@pop
ind_scores_df$HybridIndex <- ids$HybridYN %>% as.factor()

env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 2  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_manual(values=c('1' = "purple4", '0' = "yellow")) +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"%"), y = paste0("RDA2"," ", percexpl_rda2,"%"), color = "Hybrid?") +
  theme_minimal()
plot

ind_scores_df$Infec <- ids$InfectionSeverity
plot1<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = Infec), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained)"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Infection Severity") +
  theme_minimal()
plot1

ind_scores_df$cankersTrunk <- ids$AreaInfectedByCanker_Trunk_perc
plot2<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = cankersTrunk), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained)"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Area infected by cankers (trunk)") +
  theme_minimal()
plot2


#run without seedlings
# health analysis
colnames(pops_health5)
pops_health5_adults<-pops_health5 %>% filter(PlantHeight_ft>10)
dim(pops_health5_adults)
health_vars <- pops_health5_adults %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))
ancestry <- select(pops_health5_adults, "HybridYN")
ancestry$HybridYN<-as.factor(ancestry$HybridYN)
#ancestry$CLUSTER_ID<-as.factor(ancestry$CLUSTER_ID)
#ancestry <- health_merge$NH_final_assignment %>% as.factor()
envdata<-select(pops_health5_adults,"bio6_z","ForestEdge_30m_z","ph_z","PlantHeight_ft", "bio15_z", "TCC_z","TCC_5km_z")
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-pops_health5_adults[complete_rows,]
dim(ancestry_clean)
health_vars_scaled <- scale(health_vars_clean)
envdata_scaled<-scale(envdata_clean[,]) %>% as.data.frame
colnames(envdata_scaled)
predictors<-cbind(envdata_scaled,ancestry_clean)
library(vegan)
hyb_rda <- rda(health_vars_scaled ~ ., data = predictors)
summary(hyb_rda)
RsquareAdj(hyb_rda)
scores(hyb_rda, display = "species", scaling = 2)
anova(hyb_rda)
anova(hyb_rda, by = "term")

ind_scores <- scores(hyb_rda, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
#ind_scores_df$ID <- edge@ind.names
#ind_scores_df$pop <- edge@pop
ind_scores_df$HybridIndex <- ids$HybridYN %>% as.factor()

env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 2  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_manual(values=c('1' = "purple4", '0' = "yellow")) +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"%"), y = paste0("RDA2"," ", percexpl_rda2,"%"), color = "Hybrid?") +
  theme_minimal()
plot

ind_scores_df$Infec <- ids$InfectionSeverity
plot1<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = Infec), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained)"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Infection Severity") +
  theme_minimal()
plot1

ind_scores_df$cankersTrunk <- ids$AreaInfectedByCanker_Trunk_perc
plot2<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = cankersTrunk), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained)"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Area infected by cankers (trunk)") +
  theme_minimal()
plot2


