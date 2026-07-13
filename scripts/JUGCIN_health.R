#############################################################################################
################### BUTTERNUT HEALTH DATA MODELING ##########################################
#############################################################################################
# load and clean data
library(dplyr)
list.files()
health_merge<- read.csv("./pops_health6.csv")
dim(health_merge)
colnames(health_merge)
health_merge$X<-NULL
health_merge$X.1<-NULL
health_merge$x_old<-NULL
health_merge$y_old<-NULL
health_merge$NewHyb_FinalAssignment.x<-NULL
health_merge$NewHyb_FinalAssignment.y<-NULL
dim(health_merge)

##############################################################
########### START HERE #######################################
##############################################################

health_merge$CLUSTER_ID<-as.factor(health_merge$CLUSTER_ID)
summary(health_merge$CLUSTER_ID)
health_merge$NewHyb_FinalAssignment<-as.factor(health_merge$NewHyb_FinalAssignment)
summary(health_merge$NewHyb_FinalAssignment)


health_merge$HybridNH<-!health_merge$NH_final_assignment %in% c("NewHyb_Butternut","NewHyb_Heartnut")
summary(health_merge$HybridNH)
summary(health_merge$HybridNH %>% as.factor())
sum(health_merge$JC_total>0.95)
#write.csv(health_best,"health_merge_clean_040826_dupsREMOVED.csv",na="")
#health_merge<-read.csv("health_merge_clean_040826_dupsREMOVED.csv")
health_merge$PercentLiveCanopy_edit<-as.numeric(health_merge$PercentLiveCanopy)*0.01  
health_merge$CLUSTER_ID<-as.factor(health_merge$CLUSTER_ID)
health_merge$InfectionSeverity<-as.numeric(health_merge$InfectionSeverity)
health_merge$DBH_cm<-as.numeric(health_merge$DBH_cm)
health_merge$PlantHeight_ft<-as.numeric(health_merge$PlantHeight_ft)

summary(health_merge$PlantHeight_ft)
summary(health_merge$InfectionSeverity)
colnames(health_merge)
summary(health_merge$FIPS)

YesSeedling<-health_merge[health_merge$Seedling_YN=="Y",]
dim(YesSeedling)
summary(YesSeedling$PlantHeight_ft)

#health_merge$hybridYN<-ifelse(health_merge$NH_cat=="NewHyb_Butternut",0,1) %>% as.factor()
#write.csv(health_merge,"health_merge_clean_040926_dupsREMOVED.csv")
hybrids<-health_merge[health_merge$HybridNH==1,]
dim(hybrids)
pures<-health_merge[health_merge$HybridNH==0,]
dim(pures)

t.test(pures$InfectionSeverity, hybrids$InfectionSeverity, alternative = "greater")
t.test(pures$AreaInfectedByCanker_Trunk_perc, hybrids$AreaInfectedByCanker_Trunk_perc, alternative = "greater")
t.test(pures$AreaInfectedByCanker_RootFlare_prec, hybrids$AreaInfectedByCanker_RootFlare_prec, alternative = "greater")
t.test(pures$PercentLiveCanopy_edit, hybrids$PercentLiveCanopy_edit, alternative = "less")

#health_merge_VT<-health_merge_VT[health_merge_VT$CLUSTER_ID=="161",]
#write.csv(health_merge_VT, "health_merge_VT_04012926.csv")
library(lme4)
library(performance)
health_merge$fastStruc_logit <- qlogis(pmin(pmax(health_merge$JC_JP_fastStrucK2,0.001),0.999))

mod <- lm(InfectionSeverity ~ fastStruc_logit+DBH_cm, data=health_merge)
summary(mod)
plot(mod)
r2(mod)


mod2 <- lmer(InfectionSeverity ~ fastStruc_logit+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge)
summary(mod2)
plot(mod2)
r2(mod2)
plot(mod2)

mod3 <- lmer(InfectionSeverity ~ HybridNH+(1|CLUSTER_ID), data=health_merge)
summary(mod3)
plot(mod3)
r2(mod3)
plot(mod3)

colnames(health_merge)

mod3 <- lmer(AreaInfectedByCanker_Trunk_perc ~ HybridNH+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge)
summary(mod3)
plot(mod3)

mod3 <- lmer(AreaInfectedByCanker_Trunk_perc ~ NH_final_assignment+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge)
summary(mod3)
plot(mod3)
r2(mod3)

mod3 <- lmer(AreaInfectedByCanker_RootFlare_prec ~ HybridNH+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge)
summary(mod3)
plot(mod3)
r2(mod3)

colnames(health_merge)
mod3 <- lmer(InfectionSeverity ~ HybridNH+scale(PlantHeight_ft)+scale(ForestEdge_30m)+(1|CLUSTER_ID), data=health_merge)
summary(mod3)
anova(mod3)
r2(mod3)

plot(mod3)
colnames(health_merge)
mod4 <- lmer(InfectionSeverity ~ HybridNH+PlantHeight_ft+scale(wc2.1_30s_bio_6)+(1|CLUSTER_ID), data=health_merge)
summary(mod4)
anova(mod4)
r2(mod4)
plot(mod4)

colnames(health_merge)


library(betareg)
plot(x=health_merge$PlantHeight_ft,y=health_merge$InfectionSeverity)
plot(x=health_merge$PlantHeight_ft,y=health_merge$AreaInfectedByCanker_Trunk_perc)
plot(x=health_merge$PlantHeight_ft,y=health_merge$AreaInfectedByCanker_RootFlare_prec)
plot(x=health_merge$DBH_cm,y=health_merge$AreaInfectedByCanker_RootFlare_prec)
plot(x=health_merge$DBH_cm,y=health_merge$AreaInfectedByCanker_Trunk_perc)
plot(x=health_merge$DBH_cm,y=health_merge$InfectionSeverity)

######################## test using a filter (only plants over 3ft) ################################
summary(health_merge$CLUSTER_ID)
ilm<-health_merge[health_merge$CLUSTER_ID==3,]
write.csv(ilm,"JBP_IsleLaMotte_butternut.csv")
noHT<-ilm[is.na(ilm$PlantHeight_ft),]
dim(noHT)


health_merge1 <- health_merge[!is.na(health_merge$PlantHeight_ft) & 
                                health_merge$PlantHeight_ft > 3, ]

dim(health_merge1)

#mod <- lm(InfectionSeverity ~ fastStruc_logit+DBH_cm, data=health_merge1)
#summary(mod)
#plot(mod)

#mod2 <- lmer(InfectionSeverity ~ fastStruc_logit+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge1)
#summary(mod2)
#r2(mod2)
#plot(mod2)


mod3 <- lmer(AreaInfectedByCanker_RootFlare_prec ~ hybridYN+PlantHeight_ft+(1|CLUSTER_ID), data=health_merge1)
summary(mod3)
plot(mod3)

#health_vars <- health_merge1 %>%
#  select(InfectionSeverity,
#         AreaInfectedByCanker_Trunk_perc,
#         AreaInfectedByCanker_RootFlare_prec, PercentLiveCanopy_edit) %>%
#  mutate(across(everything(), as.numeric))
health_vars <- health_merge1 %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))

ancestry <- select(health_merge1, "HybridNH","CLUSTER_ID")
ancestry$HybridNH<-as.factor(ancestry$HybridNH)
ancestry$CLUSTER_ID<-as.factor(ancestry$CLUSTER_ID)
#ancestry <- health_merge$NH_final_assignment %>% as.factor()

envdata<-select(health_merge1,"bio6","ForestEdge_30m","ph","PlantHeight_ft")
colnames(health_merge1)
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-health_merge1[complete_rows,]
length(ancestry_clean)
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
ind_scores_df$HybridIndex <- ids$HybridNH

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
  scale_color_manual(values=c('TRUE' = "purple4", 'FALSE' = "yellow")) +
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

######################## test multivariate  ################################

library(dplyr)
sum(!is.na(health_merge$PercentLiveCanopy_edit))
sum(!is.na(health_merge$EpicormicBranches_FromTrunk))
sum(!is.na(health_merge$AreaInfectedByCanker_RootFlare_prec))
sum(!is.na(health_merge$AreaInfectedByCanker_Trunk_perc))
sum(!is.na(health_merge$InfectionSeverity))

sum(complete.cases(health_merge[, c(
  "PercentLiveCanopy_edit",
  "AreaInfectedByCanker_RootFlare_prec",
  "AreaInfectedByCanker_Trunk_perc",
  "InfectionSeverity"
)]))

sum(complete.cases(health_merge[,
"PercentLiveCanopy_edit"]))

sum(complete.cases(health_merge[,"InfectionSeverity"]))

sum(complete.cases(health_merge[, c(
  "PercentLiveCanopy_edit",
  "InfectionSeverity"
)]))

sum(complete.cases(health_merge1[, c(
  "AreaInfectedByCanker_RootFlare_prec",
  "AreaInfectedByCanker_Trunk_perc",
  "InfectionSeverity"
)]))


health_vars <- health_merge %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))

#health_vars <- health_merge %>%
#  select(PercentLiveCanopy,
#         InfectionSeverity,
#         AreaInfectedByCanker_Trunk_perc,
#         AreaInfectedByCanker_RootFlare_prec) %>%
#  mutate(across(everything(), as.numeric))
ancestry <- select(health_merge, "HybridNH","CLUSTER_ID")
ancestry$HybridNH<-as.factor(ancestry$HybridNH)
ancestry$CLUSTER_ID<-as.factor(ancestry$CLUSTER_ID)
#ancestry <- health_merge$NH_final_assignment %>% as.factor()

colnames(health_merge)
envdata<-select(health_merge,"bio6","ForestEdge_30m","ph","PlantHeight_ft")
colnames(health_merge)
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-health_merge[complete_rows,]
length(ancestry_clean)
health_vars_scaled <- scale(health_vars_clean)
envdata_scaled<-scale(envdata_clean[,]) %>% as.data.frame
colnames(envdata_scaled)
predictors<-cbind(envdata_scaled,ancestry_clean)
predictors$CLUSTER_ID <- droplevels(predictors$CLUSTER_ID)
library(vegan)
hyb_rda <- rda(health_vars_scaled ~ ., data = predictors)
summary(hyb_rda)
plot(hyb_rda, type="n", scaling=2)
points(hyb_rda, display="sites", pch=16)
text(hyb_rda, display="bp", col="red", scaling=2)
RsquareAdj(hyb_rda)
scores(hyb_rda, display = "species", scaling = 2)
anova(hyb_rda)
anova(hyb_rda, by = "term")

ind_scores <- scores(hyb_rda, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
#ind_scores_df$ID <- edge@ind.names
#ind_scores_df$pop <- edge@pop
ind_scores_df$HybridIndex <- ids$HybridNH

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
  scale_color_manual(values=c('TRUE' = "purple4", 'FALSE' = "yellow")) +
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

plot(x=health_merge$ph,y=health_merge$AreaInfectedByCanker_RootFlare_prec)
plot(x=health_merge$ph,y=health_merge$InfectionSeverity)
plot(x=health_merge$ph,y=health_merge$AreaInfectedByCanker_Trunk_perc)
plot(x=health_merge$bio6,y=health_merge$AreaInfectedByCanker_RootFlare_prec)
plot(x=health_merge$bio6,y=health_merge$InfectionSeverity)
plot(x=health_merge$bio6,y=health_merge$AreaInfectedByCanker_Trunk_perc)
plot(x=ids$ph,y=ids$AreaInfectedByCanker_RootFlare_prec)
plot(x=ids$ph,y=ids$InfectionSeverity)
plot(x=ids$ph,y=ids$AreaInfectedByCanker_Trunk_perc)
plot(x=ids$bio6,y=ids$AreaInfectedByCanker_RootFlare_prec)
plot(x=ids$bio6,y=ids$InfectionSeverity)
plot(x=ids$bio6,y=ids$AreaInfectedByCanker_Trunk_perc)



data1<-cbind(health_vars_scaled,envdata_scaled,ancestry_clean)
library(brms)
# transform: scale 0 and 1 slightly
model <- brm(InfectionSeverity ~ HybridNH+PlantHeight_ft+bio6+ph+ForestEdge_30m+(1|CLUSTER_ID) , data = data1, family = gaussian(), chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.99))
summary(model)

model2 <- brm(
  JC_total ~ human_footprint+bio6_z*nit_z,
  data = avg_cm,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000
)

summary(model)


################################################################
########## look into pop dynamics at different populations #####
################################################################
#health_merge<-read.csv("health_merge_clean_040826_dupsREMOVED.csv")
health_merge$PercentLiveCanopy_edit<-as.numeric(health_merge$PercentLiveCanopy)*0.01  
health_merge$CLUSTER_ID<-as.factor(health_merge$CLUSTER_ID)
health_merge$InfectionSeverity<-as.numeric(health_merge$InfectionSeverity)
health_merge$FIPS<-as.factor(health_merge$FIPS)
health_merge$DBH_cm<-as.numeric(health_merge$DBH_cm)
health_merge$PlantHeight_ft<-as.numeric(health_merge$PlantHeight_ft)

health_inds<-health_merge[!is.na(health_merge$LabID),]
sum(is.na(health_inds$LabID))
health_inds$CLUSTER_ID<-as.factor(health_inds$CLUSTER_ID)
summary(health_inds$CLUSTER_ID)
dim(health_inds)

plot(x=health_inds$JC_JP_fastStrucK2,y=health_inds$DBH_cm)
colnames(health_inds)
ggplot(health_inds, aes(x = NH_final_assignment, y = DBH_cm)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "DBH", title=paste0("All Sites; DBH by NewHybrids hybrid class (n = ", sum(!is.na(health_inds$DBH_cm)), ")"  ))
ggplot(health_inds, aes(x = NH_final_assignment, y = PlantHeight_ft)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Plant Height (ft)", title=paste0("All Sites; Plant height by NewHybrids hybrid class (n = ", sum(!is.na(health_inds$PlantHeight_ft)), ")"  ))
ggplot(health_inds, aes(x = NH_final_assignment, y = InfectionSeverity)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Infection Severity (1-5)", title=paste0("All Sites; Infection Severity by NewHybrids hybrid class (n = ", sum(!is.na(health_inds$InfectionSeverity)), ")"  ))
windows()
ggplot(health_inds, aes(x = NH_final_assignment, y = InfectionSeverity)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16)
  ) +
  labs(
    x = "Hybrid Class",
    y = "Infection Severity",
    title = paste0(
      "All Sites; Infection severity by NewHybrids hybrid class (n = ",
      sum(!is.na(health_inds$InfectionSeverity)),
      ")"
    )
  )


ggplot(health_inds, aes(x = NH_final_assignment, y = AreaInfectedByCanker_RootFlare_prec)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14)
  ) +
  labs(
    x = "Hybrid Class",
    y = "percent area infected by canker, root flare",
    title = paste0(
      "All Sites; Percent area infected by canker (root flare) by NewHybrids hybrid class (n = ",
      sum(!is.na(health_inds$AreaInfectedByCanker_RootFlare_prec)),
      ")"
    )
  )

ggplot(health_inds, aes(x = NH_final_assignment, y = AreaInfectedByCanker_Trunk_perc)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14)
  ) +
  labs(
    x = "Hybrid Class",
    y = "percent area infected by canker, trunk",
    title = paste0(
      "All Sites; Percent area infected by canker (trunk) by NewHybrids hybrid class (n = ",
      sum(!is.na(health_inds$AreaInfectedByCanker_Trunk_perc)),
      ")"
    )
  )
colnames(health_inds)
ggplot(health_inds, aes(x = NH_final_assignment, y = as.numeric(PercentLiveCanopy))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18)
  ) +
  labs(
    x = "Hybrid Class",
    y = "percent area infected by canker, root flare",
    title = paste0(
      "All Sites; Percent live canopy by NewHybrids hybrid class (n = ",
      sum(!is.na(health_inds$PercentLiveCanopy)),
      ")"
    )
  )



#isle la motte
ILM_inds<-health_inds[health_inds$CLUSTER_ID=="102",]
summary(ILM_inds$DBH_cm)
plot(x=ILM_inds$JP_fastStruc_K2_JC,y=ILM_inds$DBH_cm)
library(ggplot2)
ILM_inds
ggplot(ILM_inds, aes(x = majority_class_totalProb, y = DBH_cm)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "DBH", title=paste0("Isle La Motte; DBH by NewHybrids hybrid class (n = ", sum(!is.na(ILM_inds$DBH_cm)), ")"  ))
ggplot(ILM_inds, aes(x = majority_class_totalProb, y = PlantHeight_ft)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Plant Height (ft)", title=paste0("Isle La Motte; Plant height by NewHybrids hybrid class (n = ", sum(!is.na(ILM_inds$PlantHeight_ft)), ")"  ))

colnames(ILM_inds)
ggplot(ILM_inds, aes(x = majority_class_totalProb, y = AreaInfectedByCanker_Trunk_perc)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Area Infected by Canker", title=paste0("Isle La Motte; Area infected by canker by NewHybrids hybrid class (n = ", sum(!is.na(ILM_inds$PlantHeight_ft)), ")"  ))


#charlotte park
CP_inds<-health_inds[health_inds$CLUSTER_ID=="73",]
dim(CP_inds)
plot(x=CP_inds$JP_fastStruc_K2_JC,y=CP_inds$DBH_cm)
ggplot(CP_inds, aes(x = majority_class_totalProb, y = DBH_cm)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "DBH", title=paste0("Charlotte Park DBH by NewHybrids hybrid class (n = ", sum(!is.na(CP_inds$DBH_cm)), ")"  ))
ggplot(CP_inds, aes(x = majority_class_totalProb, y = PlantHeight_ft)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Plant Height (ft)", title=paste0("Charlotte Park; Plant height by NewHybrids hybrid class (n = ", sum(!is.na(CP_inds$PlantHeight_ft)), ")"  ))


#west chicago prairies
WCP_inds<-health_inds[health_inds$CLUSTER_ID=="4",]
dim(WCP_inds)
sum(is.na(WCP_inds$LabID))
plot(x=WCP_inds$JP_fastStruc_K2_JC,y=WCP_inds$DBH_cm)
ggplot(WCP_inds, aes(x = majority_class_totalProb, y = DBH_cm)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "DBH", title=paste0("West Chicago Prairies DBH by NewHybrids hybrid class (n = ", sum(!is.na(WCP_inds$DBH_cm)), ")"  ))
ggplot(WCP_inds, aes(x = majority_class_totalProb, y = PlantHeight_ft)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic() +
  theme(axis.text=element_text(size =14), axis.title = element_text(size = 16),plot.title=element_text(size =18)) +
  labs(x = "Hybrid Class", y = "Plant Height (ft)", title=paste0("West Chicago Prairies; Plant height by NewHybrids hybrid class (n = ", sum(!is.na(WCP_inds$PlantHeight_ft)), ")"  ))





library("glmmTMB")
library(performance)

avg_cm<-read.csv("health_merge_clean_040926_dupsREMOVED.csv")
avg_cm$hybridYN<-as.factor(avg_cm$hybridYN)
hist(avg_cm$human_footprint)
hist(log10(avg_cm$POP_SQMI_2))
hist(avg_cm$TCC)
colnames(avg_cm)
summary(avg_cm$hybridYN)
mod_beta6 <- lm(InfectionSeverity ~ hybridYN+PlantHeight_ft,data = avg_cm)
plot(ggpredict(mod_beta6, terms = "hybridYN"))
summary(mod_beta6)
r2(mod_beta6)
plot(ggpredict(mod_beta6, terms = "bio6"))
plot(ggpredict(mod_beta6, terms = "TCC"))
plot(ggpredict(mod_beta6, terms = "human_footprint"))


