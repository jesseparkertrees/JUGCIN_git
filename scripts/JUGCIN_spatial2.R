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
write.csv(newhyb_table,"./summaries/NewHybridsCategoriesSummaryTableALL.csv", row.names = FALSE)

## summarize categories by admin unit
summary_table <- pops_health7 %>%
  count(HASC_2, NewHyb_FinalAssignment) %>%
  pivot_wider(names_from = NewHyb_FinalAssignment,values_from = n,values_fill = 0)
summary_table
summary_table <- pops_health7 %>%
  count(HASC_2, NAME_2, NAME_1, NewHyb_FinalAssignment) %>%
  pivot_wider(
    names_from = NewHyb_FinalAssignment,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(
    pops_health7 %>%
      group_by(HASC_2, NAME_2, NAME_1) %>%
      summarize(mean_JC_struc = mean(JC_struc, na.rm = TRUE), .groups="drop"),
    by = c("HASC_2","NAME_2", 'NAME_1'))
summary_table <- summary_table %>% mutate(Hybrids = rowSums(across(-c(HASC_2, NAME_2, NAME_1, JC, JA, mean_JC_struc))))
summary_table <- summary_table %>% mutate(Total = (rowSums(across(-c(HASC_2, NAME_2, NAME_1, mean_JC_struc, Hybrids)))))
summary_table <- summary_table %>% mutate(Hybrid_perc = (Hybrids/Total)*100)
write.csv(summary_table,"./summaries/NewHybridsCategoriesSummaryTableCounties.csv", row.names = FALSE)

## plot on map
library(geodata)
library(sf)
can <- gadm("CAN",level = 2,path = "../JUGCIN_git_externalFiles",resolution = 2)
usa <- gadm("USA",level = 2,path = "../JUGCIN_git_externalFiles",resolution = 2)
admin1 <- rbind(usa, can)
admin1_sf<-st_as_sf(admin1)
sf::sf_use_s2(FALSE)
admin1_east <- st_crop(admin1_sf,xmin = -100,xmax = -62, ymin = 25,ymax = 51)
sf::sf_use_s2(TRUE)
admin1_east <- admin1_east %>%
  group_by(HASC_2) %>%
  summarize()
admin1_east <- admin1_east %>%
  left_join(summary_table, by = "HASC_2")
sum(admin1_east$Total, na.rm = TRUE)
sum(as.numeric(summary_table$Total[!is.na(summary_table$Total)]))

label_pts <- st_point_on_surface(admin1_east)
label_coords <- cbind(
  st_drop_geometry(label_pts),
  st_coordinates(label_pts))

plot<-ggplot(admin1_east) +
  geom_sf(aes(fill = Hybrid_perc), color = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "white") +
  theme_void() +
  labs(fill = "Percent of\ntotal that are\nhybrids") +
  geom_sf_text(data = subset(label_pts, !is.na(Total)),aes(label = Total),size = 0.5,color = "black",fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per county", size=2) +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5),
        legend.position=c(0.85,0.35))
plot
ggsave(filename="./summaries/HybridSummariesByCounty_percentHybrid.png",plot=plot, width=5, height=4, units='in', bg='white',dpi=600)

plot2<-ggplot(admin1_east) +
  geom_sf(aes(fill = mean_JC_struc), color = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "white", direction=-1) +
  theme_void() + 
  labs(fill="Average\nAncestry\nProportion") +
  geom_sf_text(data=subset(label_pts, !is.na(Total)), aes(label=Total), size=0.5, color='black', fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per county", size=2) +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5),
        legend.position=c(0.85,0.35))
plot2
ggsave(filename="./summaries/HybridSummariesByCounty_ancestryProportion.png",plot=plot2, width=5, height=4, units='in', bg='white', dpi=600)

## summarize by state
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
summary_table2 <- summary_table2 %>% mutate(Total = (rowSums(across(-c(NAME_1, mean_JC_struc, Hybrids)))))
summary_table2 <- summary_table2 %>% mutate(Hybrid_perc = (Hybrids/Total)*100)
write.csv(summary_table2,"./summaries/NewHybridsCategoriesSummaryTableStates.csv", row.names = FALSE)

## plot on map
us <- ne_states(country = "United States of America", returnclass = "sf")
ca <- ne_states(country = "Canada", returnclass = "sf")
admin <- bind_rows(us, ca)
admin <- admin %>%
  left_join(summary_table2, by = c("name" = "NAME_1"))
admin_east <- st_crop(admin, xmin = -100, xmax = -62, ymin = 25, ymax = 51)
library(sf)
label_pts <- st_point_on_surface(admin_east)
label_coords <- cbind(
  st_drop_geometry(label_pts),
  st_coordinates(label_pts))
plot3<-ggplot(admin_east) +
  geom_sf(aes(fill = Hybrid_perc), color = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "white") +
  theme_void() +
  labs(fill = "Percent of\ntotal that are\nhybrids") +
  geom_sf_text(
    data = subset(label_pts, !is.na(Total)),
    aes(label = Total),
    size = 2,
    color = "black",
    fontface='bold') +
    labs(caption = "*Numbers indicate total number of trees sampled per state")+
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5),
        legend.position=c(0.89,0.35))
plot3
ggsave(filename="./summaries/HybridSummariesByState_percentHybrid.png",plot=plot3, width=5, height=4, units='in', bg='white' )

colnames(admin_east)
plot4<-ggplot(admin_east) +
  geom_sf(aes(fill = mean_JC_struc), color = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "white", direction=-1) +
  theme_void() + 
  labs(fill="Average\nAncestry\nProportion") +
  geom_sf_text(data=subset(admin_east, !is.na(Total)), aes(label=Total), size=2, color='black', fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per state") +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5),
        legend.position=c(0.89,0.35))
plot4
ggsave(filename="./summaries/HybridSummariesByState_ancestryProportion.png",plot=plot4, width=5, height=4, units='in', bg='white' )

##################################################################

# spatial analysis
colnames(pops_health7)
pops_health7$NAME_1
summary(as.factor(pops_health7$NewHyb_FinalAssignment))
pops_health7$HybridYN<-ifelse(pops_health7$NewHyb_FinalAssignment=="JC" | pops_health7$NewHyb_FinalAssignment=="JA", 0,1 )
summary(as.factor(pops_health7$HybridYN))
summary(as.factor(pops_health7$CLUSTER_ID))


################delete?############################
#ILM<-pops_health5[pops_health5$CLUSTER_ID=="3",]
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
################delete?############################

##################################################################
#check correlation
library(tidyr)
colnames(pops_health7)
pops_health7_var<-pops_health7 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_1,wc2.1_30s_bio_11,wc2.1_30s_bio_12,wc2.1_30s_bio_3) %>% drop_na() %>% as.data.frame()
cor(pops_health7_var)
library(usdm)
usdm::vif(pops_health7_var)
pops_health7_var<-pops_health7 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_12, wc2.1_30s_bio_15) %>% drop_na() %>% as.data.frame()
usdm::vif(pops_health7_var)

#rescale data
pops_health8 <- pops_health7 |>
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
    NAME_1=as.factor(NAME_1),
    NAME_2=as.factor(NAME_2),
    HASC_2=as.factor(HASC_2),
    NLCD_1km=as.factor(NLCD_1km),
    NLCD_5km=as.factor(NLCD_5km),
    NLCD=as.factor(NLCD),
    TCC_z=scale(TCC),
    TCC_1km_z=scale(TCC_1km),
    TCC_5km_z=scale(TCC_5km)
    )

##################################################################
# spatial analysis
colnames(pops_health8)

library("glmmTMB")
library(performance)
library("ggeffects")

mod_beta <- glmmTMB(HybridYN ~ ForestEdge_30m_z*TCC_5km_z*bio12_z*nit_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta)

mod_beta2 <- glmmTMB(HybridYN ~ ForestEdge_NorAmer_custom_z*TCC_5km_z*bio12_z*nit_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta2)

mod_beta3 <- glmmTMB(HybridYN ~ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta3)

mod_beta4 <- glmmTMB(HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta4)

mod_beta5 <- glmmTMB(HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta5)

mod_beta6 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta6)

mod_beta7 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta7)

mod_beta8 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|NAME_1),family = binomial(),data = pops_health8)
summary(mod_beta8)

r2(mod_beta)

AIC(mod_beta, mod_beta2, mod_beta3, mod_beta4, mod_beta5, mod_beta6, mod_beta7)

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


