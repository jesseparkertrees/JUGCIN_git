
################################################################################################ initial plotting #########################################################################################
library(dplyr)
library(terra)
library(geodata)
library(sf)
library(ggplot2)
library(ggnewscale)
bio6<-rast("E:/env_data/data/wc2.1_30s_bio/raw_rasters_touse/wc2.1_30s_bio_6.tif")

pts<-dplyr::select(pops_filt, c("LabID","JC_total","x","y"))
pts_vect<-vect(pts, crs='epsg:4326', c('x','y'))
pts_sf <- st_as_sf(pts_vect)

#mapping ranges2
juci<-vect("E:/BUTTERNUT/data/juglcine.shp", crs='EPSG:4267')

#download admin shapefiles
#us_states <- gadm("USA", level=1, path="./data")  # Level 1 = states
#can<-gadm("Canada", level=1, path="./data")

#load admin shapefiles
us_states  <- readRDS("./data/gadm/gadm41_USA_1_pk.RDS")
can_states <- readRDS("./data/gadm/gadm41_CAN_1_pk.RDS")
north_america<-rbind(us_states, can_states)
us_sf <- st_as_sf(us_states)
can_sf <- st_as_sf(can_states)

# Clip US states to raster extent
raster_extent <- st_as_sfc(st_bbox(c(xmin=-100, xmax=-65, ymin=24, ymax=50), crs=st_crs(us_sf)))
can_sf_clipped <- st_crop(can_sf, raster_extent)
us_sf_clipped <- st_crop(us_sf, raster_extent)
single_sf <- dplyr::bind_rows(list(us_sf_clipped,can_sf_clipped))
dissolve_sf <- st_union(single_sf)
pts_df <- pts_sf %>%
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2])
bio6_crop<-crop(bio6, dissolve_sf)
df <- as.data.frame(bio6_crop, xy=TRUE)
colnames(df) <- c("x", "y", "presence")

# Add species and source info as data frames for legend mapping
juci_sf<-st_as_sf(juci)
juci_sf$species <- "J. cinerea"
juci_sf$range_label <- "Little's J. cinerea Range"

pts_df$source <- "DNA Samples"
library(ggnewscale)  # Make sure this is loaded
# Add a new column for separate fill legend
pts_df$source_fill <- "DNA Samples"
# Add a variable for legend mapping
windows()
map <- ggplot() +
  # Raster layer
  geom_raster(data = df, aes(x = x, y = y, fill = presence)) +
  scale_fill_gradientn(
    colours = c("transparent", "#440154FF", "#440154FF", "#440154FF", "#21908CFF", "#FDE725FF"),
    values = scales::rescale(c(0, 0.05, 0.15, 0.3, 0.5, 0.98)),
    na.value = "white",
    name = expression("Min. Temp. of the Coldest Month"),
    guide = guide_colorbar(
      label.theme = element_text(size = 20, family = "raleway"),
      title.theme = element_text(size = 20, family = "raleway")
    )
  ) +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() +
  # Base map
  geom_sf(data = us_sf_clipped, fill = NA, color = "black", linewidth = 0.25) +
  
  # QUBI range as purple line with its own color legend
  geom_sf(data = juci_sf, aes(color = "Little's J. cinerea Range"), fill = NA, linewidth = 1.25) +
  scale_color_manual(
    name = "",
    values = c("Little's J. cinerea Range" = "#FDE725FF"),
    labels=c(expression("Little's "*italic("J. cinerea  ")*"Range"))
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linetype = "solid", size = 1.25),
      order = 2
    )
  ) +
  # Points (pink triangle with black border)
  geom_point(
    data = pts_df,
    aes(x = X, y = Y, shape = "DNA Samples", fill = "DNA Samples"),
    size = 2,
    stroke = 0.5,
    color = "black"
  ) +
  scale_fill_manual(
    name = "",
    values = c("DNA Samples" = "magenta3")
  ) +
  scale_shape_manual(
    name = "",
    values = c("DNA Samples" = 23)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 23, fill = "magenta3", color = "black", stroke = 0.5, size = 3),
      order = 1
    ),
    shape = "none"
  ) +
  
  coord_sf(clip = "off") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.7, "cm"),
    legend.key = element_rect(fill = "transparent", colour = NA),
    text = element_text(family = "raleway"),
    legend.text = element_text(family = "raleway", size = 24),
    legend.title = element_text(family = "raleway", size = 24)
  )
map
#ggsave("butternut_DNAandRange.png",plot=map, width = 8, height = 8, bg='white')

#########################################################################################
############ spatial analysis set up and data acquisition #########################################################################################

#put in other R file

#########################################################################################
######### summary and linear modelling
#########################################################################################
# remove "pure" individuals
census_matched_hybs<-census_matched[census_matched$JC_total<0.95,] 
colnames(census_matched_hybs)
hist(census_matched_hybs$JC_total, breaks=15)

#average per site (many genotypes have exact same coordinate)
str(census_matched_hybs$FIPS)
str(census_matched_hybs$CLUSTER_ID)
str(census_matched_hybs$NLCD)
str(census_matched_hybs$TCC)

census_matched_hybs$FIPS<-as.factor(census_matched_hybs$FIPS)
census_matched_hybs$CLUSTER_ID<-as.factor(census_matched_hybs$CLUSTER_ID)


avg_cm <- census_matched_hybs %>%
  group_by(x, y) %>%
  summarise(
    n_ind = n(),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ !is.numeric(.x)), ~ first(.x)),
    .groups = "drop"
  )
dim(avg_cm)

#check correlation
library(tidyr)

colnames(avg_cm)
avg_cm_var<-avg_cm %>% select(forest_edge_distance, ForestEdge_30m, human_footprint, POP_SQMI_2, bio6, nitrogen, ph) %>% drop_na() %>% as.data.frame()
cor(avg_cm_var)
library(usdm)
usdm::vif(avg_cm_var)

#rescale data
avg_cm <- avg_cm |>
  dplyr::mutate(
    human_footprint_z = scale(human_footprint),
    ForestEdge_30m_z = scale(ForestEdge_30m),
    bio6_z = scale(bio6),
    bio1_z = scale(bio1),
    bio12_z = scale(bio12),
    bio11_z = scale(bio11),
    bio3_z = scale(bio3),
    ph_z = scale(ph),
    nit_z = scale(nitrogen),
    CLUSTER_ID=as.factor(CLUSTER_ID),
    FIPS=as.factor(FIPS),
    POP_SQMI_2_z = scale(POP_SQMI_2),
    NLCD=as.factor(NLCD),
    TCC_z=scale(TCC))
#avg_cm <- avg_cm |>
#  dplyr::mutate(
#    human_footprint_z = scale(human_footprint),
#    ForestEdge_30m_z = scale(ForestEdge_30m),
#    bio6_z = scale(bio6),
#    bio1_z = scale(bio1),
#    bio12_z = scale(bio12),
#    bio11_z = scale(bio11),
#    bio3_z = scale(bio3),
#    ph_z = scale(ph),
#    nit_z = scale(nitrogen),
#    CLUSTER_ID=as.factor(CLUSTER_ID),
#    FIPS=as.factor(FIPS),
#    POP_SQMI_2_z = scale(POP_SQMI_2))
dim(avg_cm)
hist(avg_cm$GAdmix, breaks=30)

avg_cm$GAdmix_transform <- (avg_cm$GAdmix*(nrow(avg_cm)-1)+0.5)/nrow(avg_cm)
library(lme4)
library(lmerTest)   # optional, for p-values
mod <- lmer(
  JC_total ~ human_footprint_z + bio6_z +TCC+NLCD+ (1|CLUSTER_ID),
  data = avg_cm
)
summary(mod)
plot(mod)

library("ggeffects")
plot(ggpredict(mod, terms = "human_footprint_z"))
plot(ggpredict(mod, terms = "NLCD"))
plot(ggpredict(mod, terms = "bio6_z"))
library("performance")
performance::r2(mod)


################### try beta regression #####################################



avg_cm$POP_log <- log10(avg_cm$POP_SQMI_2)
mod_beta8 <- glmmTMB(JC_total ~ scale(POP_log) + scale(TCC),family = beta_family(),data = avg_cm)
summary(mod_beta8)
r2(mod_beta8)
plot(ggpredict(mod_beta8, terms = "POP_log"))
plot(ggpredict(mod_beta8, terms = "TCC"))

#plot all together
library(ggeffects)
library(dplyr)
library(ggeffects)

pred_POP_log  <- ggpredict(mod_beta8, terms = "POP_log")
pred_TCC <- ggpredict(mod_beta8, terms = "TCC")

pred_POP_log$panel <- "Population Density (Log scale)"
pred_TCC$panel <- "Tree Canopy Cover"

pred_all <- bind_rows(pred_POP_log, pred_TCC)

coefs <- summary(mod_beta8)$coefficients$cond
rownames(coefs)<-c("Intercept","Population Density (Log scale)","Tree Canopy Cover")
pvals <- data.frame(
  panel = c("Population Density (Log scale)", "Tree Canopy Cover"),
  label = c(
    paste0("p = ", signif(coefs["Population Density (Log scale)", "Pr(>|z|)"], 3)),
    paste0("p = ", signif(coefs["Tree Canopy Cover", "Pr(>|z|)"], 3))),
  x = c(2.5,60),  # put in upper-right corner
  y = Inf)

library(showtext)
library(sysfonts)
font_add_google("Raleway", "raleway")
showtext_auto()

windows()
ggplot(pred_all, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  facet_wrap(~panel, scales = "free_x", nrow=1) +
  labs(
    x = "Predictor Value",
    y = expression(italic("J. cinerea")*" Admixture Proportion"),
    color = "Nitrogen",
    fill = "Nitrogen"
  ) +
  geom_text(
    data = pvals,
    aes(x = x, y = y, label = label),
    hjust = 1.1,
    vjust = 1.1,
    inherit.aes = FALSE,
    size = 4,
    family="raleway"
  ) +
  theme_classic(base_family = "raleway") +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 13),
    legend.position = "none"
  )



#BEST? model when filtering to JC_total<0.95
#mod_beta8 <- glmmTMB(JC_total ~ NLCD+POP_SQMI_2_z,family = beta_family(),  data = avg_cm)
#summary(mod_beta8)
#r2(mod_beta8)

#mod_beta <- glmmTMB(GAdmix ~ human_footprint_z + bio6_z +nit_z,family = beta_family(),  data = avg_cm)
#summary(mod_beta)
#r2(mod_beta)

library("DHARMa")
sim <- simulateResiduals(mod_beta8)
plot(sim)

#test if NLCD is meaningful [it is not]
#mod_full  <- glmmTMB(GAdmix ~ human_footprint_z + bio6_z + TCC_z + NLCD + (1|FIPS),  family=beta_family(), data=avg_cm)

#mod_reduced <- glmmTMB(GAdmix ~ human_footprint_z + bio6_z + TCC_z + (1|FIPS), family=beta_family(), data=avg_cm)

#anova(mod_reduced, mod_full)


plot(ggpredict(mod_beta8, terms = "POP_SQMI_2_z"))
plot(ggpredict(mod_beta8, terms = "NLCD"))

#write.csv(avg_cm,"avg_cm_031026.csv")



############################ try making admixture categorical...introgressed...or not?################
library(pROC)
library(ggplot2)
library(dplyr)
library(ggeffects)

census_matched_hybs<-census_matched[census_matched$JC_total<1,] 
colnames(census_matched_hybs)
hist(census_matched_hybs$JC_total, breaks=15)

str(census_matched_hybs$FIPS)
str(census_matched_hybs$CLUSTER_ID)
str(census_matched_hybs$NLCD)
str(census_matched_hybs$TCC)

census_matched_hybs$FIPS<-as.factor(census_matched_hybs$FIPS)
census_matched_hybs$CLUSTER_ID<-as.factor(census_matched_hybs$CLUSTER_ID)

avg_cm <- census_matched_hybs %>%
  group_by(x, y) %>%
  summarise(
    n_ind = n(),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    across(where(~ !is.numeric(.x)), ~ first(.x)),
    .groups = "drop"
  )
dim(avg_cm)

#rescale data
avg_cm <- avg_cm |>
  dplyr::mutate(
    human_footprint_z = scale(human_footprint),
    ForestEdge_30m_z = scale(ForestEdge_30m),
    bio6_z = scale(bio6),
    bio1_z = scale(bio1),
    bio12_z = scale(bio12),
    bio11_z = scale(bio11),
    bio3_z = scale(bio3),
    ph_z = scale(ph),
    nit_z = scale(nitrogen),
    CLUSTER_ID=as.factor(CLUSTER_ID),
    FIPS=as.factor(FIPS),
    POP_SQMI_2_z = scale(POP_SQMI_2),
    NLCD=as.factor(NLCD),
    TCC_z=scale(TCC))

avg_cm$hybrid_cat<-ifelse(avg_cm$JC_total>0.95,0 ,1) %>% as.factor()#hybrids are 1, pures are 0
str(avg_cm$hybrid_cat)
sum(avg_cm$hybrid_cat==0)
sum(avg_cm$hybrid_cat==1)

latGLM<-glm(hybrid_cat ~ y, family = binomial, data=avg_cm)
summary(latGLM)
plot(ggeffects::ggpredict(latGLM))

bin_mod <- glm(hybrid_cat ~ human_footprint_z + ForestEdge_30m_z+POP_SQMI_2_z+NLCD+TCC_z,
           data = avg_cm,
           family = binomial)
summary(bin_mod)
1 - (bin_mod$deviance / bin_mod$null.deviance)
#forest edge and TCC and NLCD don't matter

#add climate variables
avg_cm <- avg_cm[!is.na(avg_cm$nit_z), ]
dim(avg_cm)
avg_cm$human_footprint_z <- as.numeric(avg_cm$human_footprint_z)
avg_cm$bio6_z <- as.numeric(avg_cm$bio6_z)
avg_cm$nit_z <- as.numeric(avg_cm$nit_z)

bin_mod2 <- glm(hybrid_cat ~ human_footprint+bio6_z*nit_z+y,
               data = avg_cm,
               family = binomial)
summary(bin_mod2)
1 - (bin_mod2$deviance / bin_mod2$null.deviance)
pred_hf <- ggpredict(bin_mod2, terms = "human_footprint [all]")
plot(pred_hf)
pred_y <- ggpredict(bin_mod2, terms = "y [all]")
plot(pred_y)
plot(ggeffects::ggpredict(bin_mod2))

bin_mod_realval <- glm(hybrid_cat ~ human_footprint + bio6*nitrogen,
                data = avg_cm,
                family = binomial)

roc_obj <- roc(avg_cm$hybrid_cat, fitted(bin_mod2))
auc(roc_obj)


library(showtext)
library(sysfonts)
font_add_google("Raleway", "raleway")
showtext_auto()
ggplot(pred_hf, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    x = "Human Footprint",
    y = "Probability of Hybrid"
  ) +
  theme_classic() +
  theme(text = element_text(family = "raleway"),axis.title = element_text(size = 14),axis.text = element_text(size = 12))

par(mfrow=c(2,2))
plot(bin_mod2)
avg_cm[c("523","474","7","478"),"LabID"]

#check for over dispersion
deviance(bin_mod2) / df.residual(bin_mod2) #if around 1, no overdispersion

pred_int <- ggpredict(bin_mod_realval, terms = c("bio6", "nitrogen"))
plot(pred_int)

#plot all together

pred_hf  <- ggpredict(bin_mod_realval, terms = "human_footprint [all]")
pred_bio6 <- ggpredict(bin_mod_realval, terms = "bio6 [all]")
pred_nit<- ggpredict(bin_mod_realval, terms = "nitrogen [all]")
pred_int <- ggpredict(
  bin_mod_realval,
  terms = c("bio6 [all]", "nitrogen [3.97,7.04,10.11]")
)

pred_hf$panel <- "Human Footprint"
pred_bio6$panel <- "Min Winter Temp"
pred_nit$panel <- "Soil Nitrogen"

pred_int$panel <- "Temp. × Nit."
pred_int$group <- factor(pred_int$group, labels = c("Low N", "Average N", "High N"))
pred_hf$group <- NA
pred_bio6$group <- NA
pred_nit$group <- NA

pred_all <- bind_rows(pred_hf, pred_bio6, pred_nit, pred_int)

coefs <- summary(bin_mod2)$coefficients
rownames(coefs)<-c("Intercept","Human Footprint","Min Temp. of Coldest Month", "Soil Nitrogen","Temp. × Nit.")
pvals <- data.frame(
  panel = c("Human Footprint", "Min Winter Temp", "Soil Nitrogen","Temp. × Nit."),
  label = c(
    paste0("p = ", signif(coefs["Human Footprint", "Pr(>|z|)"], 3)),
    paste0("p = ", signif(coefs["Min Temp. of Coldest Month", "Pr(>|z|)"], 3)),
    paste0("p = ", signif(coefs["Soil Nitrogen", "Pr(>|z|)"], 3)),
    paste0("p (Temp. × Nit.) = ", signif(coefs["Temp. × Nit.", "Pr(>|z|)"], 3))
  ),
  x = c(35,-5,15,0),  # put in upper-right corner
  y = Inf)
pred_all$group_plot <- pred_all$group
pred_all$group_plot[is.na(pred_all$group_plot)] <- "Single"  # placeholder for coloring

windows()
plot<-ggplot(pred_all, aes(x = x, y = predicted)) +
  # single-predictor lines: red, no legend
  geom_line(
    data = subset(pred_all, is.na(group)),
    aes(y = predicted),
    color = "red",
    size = 0.75
  ) +
  geom_ribbon(
    data = subset(pred_all, is.na(group)),
    aes(ymin = conf.low, ymax = conf.high),
    fill = "red",
    alpha = 0.15
  ) +
  # interaction lines: colored, appear in legend
  geom_line(
    data = subset(pred_all, !is.na(group)),
    aes(color = group),
    size = 0.75
  ) +
  geom_ribbon(
    data = subset(pred_all, !is.na(group)),
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15
  ) +
  facet_wrap(~panel, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("Low N" = "#1b9e77", "Average N" = "#d95f02", "High N" = "#7570b3")) +
  geom_text(
    data = pvals,
    aes(x = x, y = y, label = label),
    hjust = 1.1,
    vjust = 1.1,
    inherit.aes = FALSE,
    size = 4,
    family="raleway"
  ) +
  scale_fill_manual(values = c("Low N" = "#1b9e77", "Average N" = "#d95f02", "High N" = "#7570b3")) +
  labs(
    x = "Predictor Value",
    y = "Predicted Probability of Hybrid",
    color = "Nitrogen",
    fill = "Nitrogen"
  ) +
  theme_classic(base_family = "raleway")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.7, "cm"),
    legend.key = element_rect(fill = "transparent", colour = NA),
    text = element_text(family = "raleway"),
    legend.text = element_text(family = "raleway", size = 16),
    legend.title = element_text(family = "raleway", size = 16)
  )
plot

#ggsave("binom_glm_preds.png",plot=plot, width=6, height=3, bg="white")

avg_cm$hybrid_cat
#try a Bayesian approach
library(brms)
model <- brm(JC_total ~ human_footprint+bio6_z*nit_z + (1|CLUSTER_ID), data = avg_cm, family = Beta(), chains = 4, cores = 4, iter = 4000)

model2 <- brm(
  JC_total ~ human_footprint+bio6_z*nit_z,
  data = avg_cm,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000
)

summary(model)
library(posterior)
draws <- as_draws_df(model)
mean(draws$b_human_footprint < 0)
quantile(draws$b_human_footprint, c(0.025, 0.5, 0.975))
library(bayesplot)
mcmc_areas(draws, pars = "b_human_footprint")
plot(model)


############################ try GAM #######################################
library(mgcv)
mod_gam <- gam(
  JC_total ~ s(human_footprint_z, k=6) + s(bio6_z) + s(nit_z) + s(TCC_z)  +  s(CLUSTER_ID, bs = "re") + s(FIPS, bs = "re"),,
  family = betar(link = "logit"),
  data = avg_cm
)
summary(mod_gam)

#control for spatial autocorrelation
mod_gam2 <- gam(
  GAdmix ~ s(human_footprint_z) +
    s(bio6_z) +
    s(nit_z) +
    s(x, y),
  family = betar(link="logit"),
  data = avg_cm
)
summary(mod_gam2)

mod_gam <- gam(
  GAdmix ~ 
    s(human_footprint_z) +
    s(bio6_z) +
    s(nit_z) +
    s(CLUSTER_ID, bs = "re") +
    s(FIPS, bs = "re"),
  family = betar(link="logit"),
  data = avg_cm
)

summary(mod_gam)

plot(mod_gam, select=1, shade=TRUE)



plot(y=census_matched_hybs$JC1, x=census_matched_hybs$JA.x)
plot(y=census_matched_hybs$JC2, x=census_matched_hybs$JA.x)
model_jc1 <- lm(JC1 ~ JA.x, data = census_matched_hybs)
summary(model_jc1)

model_jc2 <- lm(JC2 ~ JA.x, data = census_matched_hybs)
summary(model_jc2)

#add ancestry index
census_matched_hybs$ancestry_index<-census_matched_hybs$JC1-census_matched_hybs$JC2
plot(x=census_matched_hybs$ancestry_index, y=census_matched_hybs$JA.x)
model_AncInd <- lm(JA.x ~ ancestry_index, data = census_matched_hybs)
summary(model_AncInd)


census_matched_hybs$JA_adj <- pmin(pmax(census_matched_hybs$JA.x, 0.001), 0.999)
census_matched_hybs$JA_logit <- qlogis(census_matched_hybs$JA_adj)

model <- lm(JA_logit ~ ancestry_index, data = census_matched_hybs)
summary(model)
ggplot(census_matched_hybs, aes(ancestry_index, JA.x)) +
  geom_point() +
  geom_smooth(method = "lm")

#ternary plot
library(ggtern)
ggtern(data = census_matched_hybs,
       aes(JC1, JC2, JA.x, color = JA.x)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  theme_bw()

#pca
census_matched$anc_prop<-census_matched$JC1-census_matched$JC2
qmat<-select(census_matched, c(LabID, JA.x, anc_prop, GAdmix))
#qmat<-select(census_matched, c(LabID,JC1, JC2, JA.x))
qmat<-na.omit(qmat)

colnames(census_matched)
predis<-census_matched[,c(1,5:28,31,41,48,50)]
predis<-na.omit(predis)
predis<-predis[predis$LabID %in% qmat$LabID,]
qmat<-qmat[qmat$LabID %in% predis$LabID,]
dim(predis)
dim(qmat)
predis$LabID<-NULL
qmat$LabID<-NULL

rad1<-rda(X=qmat$GAdmix,Y=predis)
summary(rad1)

rad2<-rda(qmat~ ., data=predis)
summary(rad2)

rad<-rda(qmat$GAdmix~POP_SQMI_2+ForestEdge_30m+ph+bio6+human_footprint+nitrogen, data=predis)
summary(rad)


vif.cca(rad)
RsquareAdj(rad)
anov<-anova.cca(rad, by="axis", permutations=99)
anov
#only first 2 RDA sig
# Test each environmental variable separately
anov2<-anova(rad, by = "term", permutations = 999)
anov2
#

summary(rad)$concont
screeplot(rad)

plot(rad, scaling=3)  ## d.fault is axes 1 and 2
#select how many axes to keep
axes<-2

library(reshape2)
library(ggplot2)

# Extract SNP (species) loadings on RDA1
loadings <- scores(rad, display = "species", scaling = 3, choices=2)
# Optional: histogram of loadings to visualize distribution
# Choose threshold z (e.g., 3 standard deviations)
screeplot(rad)
load.rda <- summary(rad)$species[,1:2]

env_scores <- scores(rad, display = "bp", scaling = 3)
env_loadings <- summary(rad)$biplot
round(env_loadings, 3)

ind_scores <- scores(rad, display = "sites", scaling = 3, choices=1:axes)
ind_scores_df <- as.data.frame(ind_scores[, 1:axes])
#ind_scores_df$ID <- g@ind.names

hi <- qmat$JA.x
ind_scores_df$HybridIndex <- hi

env_arrows <- scores(rad, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 10  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(rad)
perc_cons<-round((rad$CCA$tot.chi/rad$tot.chi)*100, 2)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
env_arrows_df$Variable
env_arrows_df$Variable2<-c("Pop per sq mi","Forest edge 30m","ph","bio6","human footprint","nitrogen")
env_arrows_df$xpos<-env_arrows_df$RDA1*1
env_arrows_df$ypos<-env_arrows_df$RDA2*1

plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_gradient(low = "yellow", high = "purple4") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey60") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = xpos, y = ypos, label = Variable2),
            color = "black", vjust = 1, size = 3, family="mono", fontface="bold") +
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained) [total constrained = ",perc_cons,"%]"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Hybrid Index") +
  theme_minimal(base_family = "mono") 

plot






















#GAM model
library(mgcv)
mod_gam <- gam(
  GAdmix ~ 
    s(x, y) + 
    bio12 + human_footprint + nitrogen + ph + ForestEdge_30m,
  data = avg_cm,
  family = gaussian(), method = 
    "REML"
)
summary(mod_gam)
vis.gam(mod_gam, view = c("x","y"), plot.type = "contour")

#do beta regression of variables and hybrid index
library(betareg)
mod<-betareg(GAdmix ~ human_footprint+ForestEdge_30m+bio6+bio1+P1920, data = avg_cm)
s<-summary(mod)
s$coefficients$mean[2,4]
# Extract R-squared and p-value
r2_value <- summary(mod)$pseudo.r.squared
p_value <- summary(mod)$coefficients$mean[2,4]  # P-value for V4
# Create text label for plot
lm_text <- paste0("Pseudo R² = ", round(r2_value, 3), 
                  "\nP-value = ", format.pval(p_value, digits = 3, eps = 0.001))
# Create scatter plot with regression line and stats annotation
plot<-ggplot(census_matched, aes(x = var_x, y = var_y)) +
  geom_point() +  # Color points by "pop"
  geom_line(color='black',show.legend=FALSE,size=0.5,aes(y = predict(mod, census_matched), linetype = "logit")) +
  xlab("Predictor") +
  ylab("Butternut Ancestry Proportion") +
  theme_minimal() +
  labs(title = "Ancestry proportion as function of predictor", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))+  
  annotate("text", 
           x = max(var_x) * 0.4, 
           y = max(var_y) * 0.6, 
           label = lm_text, 
           size = 3, family = "mono", color="black")
plot


#do beta regression of variables and hybrid index


library(betareg)
var_y<-hybs_avg$JC_total
var_x<-hybs_avg$x
mod <- betareg(var_y ~ var_x)
s<-summary(mod)
s$coefficients$mean[2,4]
# Extract R-squared and p-value
r2_value <- summary(mod)$pseudo.r.squared
p_value <- summary(mod)$coefficients$mean[2,4]  # P-value for V4
# Create text label for plot
lm_text <- paste0("Pseudo R² = ", round(r2_value, 3), 
                  "\nP-value = ", format.pval(p_value, digits = 3, eps = 0.001))
# Create scatter plot with regression line and stats annotation
plot<-ggplot(hybs_avg, aes(x = var_x, y = var_y)) +
  geom_point() +  # Color points by "pop"
  geom_line(color='black',show.legend=FALSE,size=0.5,aes(y = predict(mod, hybs_avg), linetype = "logit")) +
  xlab("Distance from Forest Edge") +
  ylab("Butternut Ancestry Proportion") +
  theme_minimal() +
  labs(title = "Ancestry Proportion as function of distance to forest edge", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))+  
  annotate("text", 
           x = max(var_x) * 0.4, 
           y = max(var_y) * 0.6, 
           label = lm_text, 
           size = 3, family = "mono", color="black")
plot


mod <- glm(JC_total ~ x, family = quasibinomial, data = ints2)
summary(mod)


#################################################################################################

##summarize by state
dat<-vect(pops_filt, crs='epsg:4326', c('x','y'))
#load admin shapefiles
us_states  <- readRDS("./data/gadm/gadm41_USA_1_pk.RDS")
can_states <- readRDS("./data/gadm/gadm41_CAN_1_pk.RDS")
north_america<-rbind(us_states, can_states)
us_sf <- st_as_sf(us_states)
can_sf <- st_as_sf(can_states)

par(mfrow=c(1,1))
plot(north_america)
points(dat)
ints<-terra::intersect(dat, north_america)
ints2<-as.data.frame(ints, geom='XY')
colnames(ints2)
summary_by_state <- ints2 %>%
  group_by(NAME_1) %>%
  summarise(
    n = n(),
    mean_JC = mean(JC_total, na.rm = TRUE),
    sd_JC = sd(JC_total, na.rm = TRUE),
    min_JC = min(JC_total, na.rm = TRUE),
    max_JC = max(JC_total, na.rm = TRUE),
    mean_LAT = mean(y, na.rm = TRUE)
  )
summary_by_state
write.csv(summary_by_state, file="summary_by_states_03052026.csv")

#get different fonts
#fonts<-print(systemfonts::system_fonts())
#View(subset(systemfonts::system_fonts(), grepl("Avenir", family)))
#library('showtext')
#font_add("Perpetua", "C:/WINDOWS/Fonts/PER_____.TTF")
#showtext_auto()

library(showtext)
library(sysfonts)
font_add_google("Raleway", "raleway")
showtext_auto()

states_plot<-ggplot(summary_by_state, aes(x = reorder(NAME_1, mean_LAT), y = mean_JC)) +
  geom_col() +
  geom_text(aes(label = paste0("n = ", n)), hjust = -0.2, size = 3, family='montserrat') +
  coord_flip() +
  labs(
    x = "State",
    y = "Mean Butternut Ancestry",
    title = "Average Butternut Ancestry by State (ordered by Latitude)"
  ) +
  theme_minimal(base_family = "montserrat")
#ggsave(filename="states_summary_plot.png",plot=states_plot, bg='white',width = 3.5, height = 1.7)


#summarize by hybrid category
summary_by_state_binomial <- ints2 %>%
  mutate(
    JC_quarter = cut(
      JC_total,
      breaks = c(0, 0.95, 1),
      include.lowest = TRUE,
      labels = c("0–0.95", "0.95–1")
    )
  ) %>%
  group_by(NAME_1) %>%
  summarise(
    n = n(),
    mean_JC = mean(JC_total, na.rm = TRUE),
    sd_JC = sd(JC_total, na.rm = TRUE),
    min_JC = min(JC_total, na.rm = TRUE),
    max_JC = max(JC_total, na.rm = TRUE),
    mean_LAT = mean(y, na.rm = TRUE),
    n_0_95 = sum(JC_quarter == "0–0.95", na.rm = TRUE),
    n_95_1 = sum(JC_quarter == "0.95–1", na.rm = TRUE),
    )

#perform Fisher's exact test on hybrid vs non hybrid columns
summary_by_state_binomial_filt<-summary_by_state_binomial[summary_by_state_binomial$n>5,]
binom_mat<-cbind(summary_by_state_binomial_filt$n_0_95, summary_by_state_binomial_filt$n_95_1) %>% as.matrix()
rownames(binom_mat)<-summary_by_state_binomial_filt$NAME_1
fisher.test(binom_mat, workspace=200000000, simulate.p.value = TRUE, B=200000)

#binom GLM to see effect of latitude
latGLM<-glm(hybrid_cat ~ y, family = binomial, data=avg_cm)
summary(latGLM)
plot(ggeffects::ggpredict(latGLM))
1 - (latGLM$deviance / latGLM$null.deviance)
roc_obj_y <- roc(avg_cm$hybrid_cat, fitted(latGLM))
auc(roc_obj_y)


#summarize by tenths
summary_by_state_tenths <- ints2 %>%
  mutate(
    JC_quarter = cut(
      JC_total,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      labels = c(
        "0–0.1","0.1–0.2","0.2–0.3","0.3–0.4","0.4–0.5",
        "0.5–0.6","0.6–0.7","0.7–0.8","0.8–0.9","0.9–1.0"
      )
    )
  ) %>%
  group_by(NAME_1) %>%
  summarise(
    n = n(),
    mean_JC = mean(JC_total, na.rm = TRUE),
    sd_JC = sd(JC_total, na.rm = TRUE),
    min_JC = min(JC_total, na.rm = TRUE),
    max_JC = max(JC_total, na.rm = TRUE),
    mean_LAT = mean(y, na.rm = TRUE),
    
    n_0_to_10 = sum(JC_quarter == "0–0.1", na.rm = TRUE),
    n_10_to_20 = sum(JC_quarter == "0.1–0.2", na.rm = TRUE),
    n_20_to_30 = sum(JC_quarter == "0.2–0.3", na.rm = TRUE),
    n_30_to_40 = sum(JC_quarter == "0.3–0.4", na.rm = TRUE),
    n_40_to_50 = sum(JC_quarter == "0.4–0.5", na.rm = TRUE),
    n_50_to_60 = sum(JC_quarter == "0.5–0.6", na.rm = TRUE),
    n_60_to_70 = sum(JC_quarter == "0.6–0.7", na.rm = TRUE),
    n_70_to_80 = sum(JC_quarter == "0.7–0.8", na.rm = TRUE),
    n_80_to_90 = sum(JC_quarter == "0.8–0.9", na.rm = TRUE),
    n_90_to_100 = sum(JC_quarter == "0.9–1.0", na.rm = TRUE)
  )
#perform Fisher's exact test on hybrid vs non hybrid columns
summary_by_state_filt<-summary_by_state_tenths[summary_by_state_tenths$n>10,]
fish_mat<-summary_by_state_filt[,8:17] %>% as.matrix()
rownames(fish_mat)<-summary_by_state_filt$NAME_1
fisher.test(fish_mat, workspace=200000000, simulate.p.value = TRUE, B=200000)


write.csv(summary_by_state_tenths, file="summary_by_states_tenths_03052026.csv")


ggplot(ints2, aes(x = JC_total)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.1),
    color = "black",
    fill = "steelblue"
  ) +
  facet_wrap(~ NAME_1, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Butternut Ancestry (Hybrid Index)",
    y = "Count",
    title = "Distribution of Butternut Ancestry by State"
  ) +
  theme_minimal()










allex<-na.omit(allex)
colnames(allex)
allext<-allex[3:31]
#remove correlated variables
library("usdm")
vifcor(allext, th=0.85)
vstep<-vifstep(allext, th=10, keep=c("bio10"))
selpreds<-preds[[vstep@results$Variables]]
extract_locs2<-extract(selpreds, locs, xy=TRUE)
allex2<-cbind(coords$ID, extract_locs2)
allex2<-na.omit(allex2)
colnames(allex2)





#####################################################################################################
preds<-rast("E:/env_data/data/predictors_30s_easternUS.tif")
preds<-preds[[names(preds) != "cec" & names(preds) != "elevation" & names(preds) != "sand"]]
names(preds)
plot(preds$bio1)
points(pts_vect)


locs<-vect(coords, crs=crs(preds), c('x','y'))
extract_locs<-extract(preds, locs, xy=TRUE)
allex<-cbind(coords$ID, extract_locs)
allex<-na.omit(allex)
colnames(allex)
allext<-allex[3:31]
#remove correlated variables
library("usdm")
vifcor(allext, th=0.85)
vstep<-vifstep(allext, th=10, keep=c("bio10"))
selpreds<-preds[[vstep@results$Variables]]
extract_locs2<-extract(selpreds, locs, xy=TRUE)
allex2<-cbind(coords$ID, extract_locs2)
allex2<-na.omit(allex2)
colnames(allex2)
#allext2<-allex2[3:8]
allext<-allex2[3:14]
allext
length(allext$bio10)
allex$`coords$ID`==g@ind.names
length(g@ind.names)






#####################################################################################################################
############## get rasters ##########################################################################################
#####################################################################################################################
setwd("E:/env_data/data/wc2.1_30s_bio/raw_rasters_touse/")

bio3<-rast("./wc2.1_30s_bio_3.tif")
bio10<-rast("./wc2.1_30s_bio_10.tif")
bio15<-rast("./wc2.1_30s_bio_15.tif")
bio18<-rast("./wc2.1_30s_bio_18.tif")
bio1<-rast("./wc2.1_30s_bio_1.tif")
bio2<-rast("./wc2.1_30s_bio_2.tif")
bio4<-rast("./wc2.1_30s_bio_4.tif")
bio5<-rast("./wc2.1_30s_bio_5.tif")
bio6<-rast("./wc2.1_30s_bio_6.tif")
bio7<-rast("./wc2.1_30s_bio_7.tif")
bio8<-rast("./wc2.1_30s_bio_8.tif")
bio9<-rast("./wc2.1_30s_bio_9.tif")
bio10<-rast("./wc2.1_30s_bio_10.tif")
bio11<-rast("./wc2.1_30s_bio_11.tif")
bio12<-rast("./wc2.1_30s_bio_12.tif")
bio13<-rast("./wc2.1_30s_bio_13.tif")
bio14<-rast("./wc2.1_30s_bio_14.tif")
bio16<-rast("./wc2.1_30s_bio_16.tif")
bio17<-rast("./wc2.1_30s_bio_17.tif")
bio19<-rast("./wc2.1_30s_bio_19.tif")



#####################################################################################################################
#####################################################################################################################
#####################################################################################################################




#hyb_rda_bio11 <- rda(g2_matrix ~ ., data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio11+bio7+slope+nit, data = allext)
#hyb_rda <- rda(g2_matrix ~ bio11+bio7+bio5+bio2+slope+nit, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio15+bio16+bio6+bio5+slope, data = allext)
#hyb_rda <- rda(g2_matrix ~ bio15+bio16+bio6+bio5+slope+nit, data = allext)
#hyb_rda <- rda(g2_matrix ~ bio13+bio15+bio3+bio7+bio8+slope+ph+ocd+TPI+TRI, data = allext)
#hyb_rda <- rda(g2_matrix_imp ~ bio6+bio5+bio15+bio16, data = allext)
#hyb_rda <- rda(g2_matrix ~ bio6, data = allext)
#hi
#hyb_rda <- rda(g2_matrix_imp ~ bio6 + bio5 + Condition(hi), data = allext)

#removed variables with vifstep th10 ( 
#Variables      VIF
#1      bio10 9.520761
#2      bio15 5.837862
#3      bio18 7.392154
#4       bio2 8.623756
#5       bio8 3.541612
#6      slope 7.143252
#7        nit 7.496956
#8        soc 7.390529
#9       sand 8.867468
#10       ocd 2.978820
#11      cfvo 2.752021
#12       TPI 2.644676
#13       TRI 6.394199)
#########################################################################################

hyb_rda <- rda(g2_matrix ~ ., data = allext)

#hyb_rda <- rda(gsub_matrix ~ ., data = allext)
#hyb_rda<-rda(formula = g2_matrix ~ bio5 + bio6 + bio18 + bio2 + slope+ ocd +cfvo + TPI + TRI, data = allext)
vif.cca(hyb_rda)
RsquareAdj(hyb_rda)
#17.307 r squared, 8.197 adj r squ
anov<-anova.cca(hyb_rda, by="axis", permutations=99)
#all 12 RDA axes signficant, first 4 are greater than 2
# Test each environmental variable separately
anov2<-anova(hyb_rda, by = "term", permutations = 999)
#all 12 variables are sig.

summary(hyb_rda)$concont
screeplot(hyb_rda)

plot(hyb_rda, scaling=3)  ## d.fault is axes 1 and 2
#select how many axes to keep
axes<-4

library(reshape2)
library(ggplot2)
po <- read.csv("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/141_2spec_qs.csv", header=TRUE)
pop<-cbind(po$Individual,po$lon, po$lat, po$pop5) %>% as.data.frame()
pop<-pop[(pop$V1 %in% edge@ind.names),]
po<-cbind(po$Individual,po$lon, po$lat) %>% as.data.frame()
po<-po[(po$V1 %in% edge@ind.names),]

###########################################################
############# get candidate SNPs and plot biplot
##########################################################

# Extract SNP (species) loadings on RDA1
snp_loadings <- scores(hyb_rda, display = "species", scaling = 3, choices=1:axes)
# Optional: histogram of loadings to visualize distribution
# Choose threshold z (e.g., 3 standard deviations)
screeplot(hyb_rda)
load.rda <- summary(hyb_rda)$species[,1:2]

# Optional: histogram of loadings to visualize distribution for each axis
par(mfrow = c(2, 2))
for (i in 1:axes) {
  hist(snp_loadings[, i],
       main = paste("SNP loadings on RDA", i),
       xlab = "Loading", breaks=100, xlim=c(-0.15,0.15))
}

par(mfrow = c(1, 1))
# Identify outliers for each of the first 6 RDA axes
# Identify outliers based on loading extremes
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]  # return names of outlier loci
}
axes<-1
cand_list <- list()
for (i in 1:axes) {
  cand_list[[i]] <- outliers(snp_loadings[, i], 3.5)
}

# Combine candidate loci from all 4 axes and remove duplicates
hyb_rda.cand <- unique(unlist(lapply(cand_list, names)))
# Check number of candidate loci
length(hyb_rda.cand)
#2449


env_scores <- scores(hyb_rda, display = "bp", scaling = 3)
env_loadings <- summary(hyb_rda)$biplot
round(env_loadings, 3)



ind_scores <- scores(hyb_rda, display = "sites", scaling = 3, choices=1:axes)
ind_scores_df <- as.data.frame(ind_scores[, 1:axes])
ind_scores_df$ID <- g@ind.names

est2<-read.csv("est2_5.csv")
est<-est2[(est2$hi.INDLABEL %in% g@ind.names),]
hi <- est$hi.h_posterior_mode[match(g@ind.names, est$hi.INDLABEL)]
ind_scores_df$HybridIndex <- hi
# OR use sNMF based hybrid indices
est2<-read.csv("q_scores_2spec_snmf.csv")
est<-est2[(est2$Individual %in% g@ind.names),]

hi <- est$oc[match(g@ind.names, est$Individual)] %>% as.data.frame()
hi$Individual<-g@ind.names
colnames(hi)<-c("hi","Individual")
ind_scores_df$HybridIndex <- hi$hi

env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 10  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda)
perc_cons<-round((hyb_rda$CCA$tot.chi/hyb_rda$tot.chi)*100, 2)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
env_arrows_df$Variable
env_arrows_df$Variable2<-c("Mean Temp. Warmest Quarter","Precip. \nSeasonality","Precip. of Warmest Quarter","Mean Diurnal Range","Mean Temp. \nWettest \nQuarter","Slope","Total \nNitrogen","Soil \nOrganic \nCarbon","Organic \nCarbon \nDensity","Course \nFragments","Topo. \nPosition \nIndex","Topo. \nRoughness \nIndex")
env_arrows_df$xpos<-env_arrows_df$RDA1*1
env_arrows_df$ypos<-env_arrows_df$RDA2*1
env_arrows_df$xpos[5]<-env_arrows_df$xpos[5]-0.2
env_arrows_df$xpos[11]<-env_arrows_df$xpos[11]+0.3
env_arrows_df$ypos[11]<-env_arrows_df$ypos[11]+0.3
env_arrows_df$xpos[12]<-env_arrows_df$xpos[12]+0.3
env_arrows_df$ypos[12]<-env_arrows_df$ypos[12]+0.3
env_arrows_df$xpos[10]<-env_arrows_df$xpos[10]
env_arrows_df$ypos[10]<-env_arrows_df$ypos[10]+0.2

plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_gradient(low = "yellow", high = "purple4") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey60") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = xpos, y = ypos, label = Variable2),
            color = "black", vjust = 1, size = 3, family="mono", fontface="bold") +
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained) [total constrained = ",perc_cons,"%]"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Hybrid Index") +
  theme_minimal(base_family = "mono") 

plot
ggsave("10_14_25_RDAplot_vifstep10_white_sNMFancestry.png",plot=plot, width = 8, height = 8, bg='white')


###########################################################
############# run IBD on subset of loci (gghybrid and RDA genes)
##########################################################
g<-gl.compliance.check(g)
g@ind.names
#g<-gl.drop.ind(g,c("qbicolor","qlyrata"))
#rownames(g@other$latlon) <- indNames(g)
g@loc.names
gsub<-gl.keep.loc(g,c("NC_044904.1_15750202", "NC_044904.1_36123896", "NC_044906.1_59613866", "NC_044909.1_44793175" ,"NC_044909.1_44857771","NC_044909.1_51361006", "NC_044910.1_34266153", "NC_044910.1_35898948", "NC_044911.1_64048033", "NC_044912.1_11153726","NC_044913.1_38831166")) %>% gl.compliance.check() %>% gl.filter.monomorphs()

tRNA pseudouridine synthase A
NC_044905.1_13384374
mitogen actived 3 kinase
NC_044904.1_36125034
gsub<-gl.keep.loc(g,c("NC_044905.1_13384374", "NC_044904.1_36125034")) %>% gl.compliance.check() %>% gl.filter.monomorphs()



#gsub_matrix<-as.matrix(gsub)
gsub1<-gl.drop.pop(gsub, c("CORE-IL1", "CORE-IN1", "CORE-IN2", "CORE-IN3", "CORE-IN4", "CORE-MI1", "CORE-OH1", "CORE-OH2","P0")) %>% gl.filter.monomorphs()
gsub1@other$latlon$lat<-as.numeric(gsub1@other$latlon$lat)
gsub1@other$latlon$lon<-as.numeric(gsub1@other$latlon$lon)
length(gsub1@other$latlon$lon)
gsub1@other$latlon <- gsub1@other$latlon[indNames(gsub1), ]
gsub_matrix<-as.matrix(gsub1)


ibd<-gl.ibd(gsub1, coordinates=gsub1@other$latlon)


###########################################################
############# COOL PLOTS OF ALLELE STATE CHANGE WITH ENV
##########################################################
prot_func_OC<-read.csv("prot_func_OC_hybrda.csv")

allext$Individual<-gsub@ind.names
colnames(allext)
#est2<-read.csv("est2_5.csv")
#est<-est2[(est2$hi.INDLABEL %in% g@ind.names),]
#hi <- est$hi.h_posterior_mode[match(g@ind.names, est$hi.INDLABEL)]
#hi<-as.data.frame(hi)
est2<-read.csv("q_scores_2spec_snmf.csv")
est<-est2[(est2$Individual %in% g@ind.names),]
hi <- est$oc[match(g@ind.names, est$Individual)] %>% as.data.frame()
hi$Individual<-g@ind.names
colnames(hi)<-c("hi","Individual")
# 1. Melt gsub_matrix into long format for ggplot
library(reshape2)
gsub_long <- reshape2::melt(as.matrix(gsub_matrix))
colnames(gsub_long) <- c("Individual", "SNP", "Allele")

# 2. Add hybrid index and environment data
gsub_long$HybridIndex <- hi$hi[match(gsub_long$Individual, hi$Individual)]
gsub_long$EnvVar <- allext$ph[match(gsub_long$Individual, allext$Individual)]  # example with bio5

# 3. Merge SNP functional annotations (optional, for better facet labels)
gsub_long <- merge(gsub_long,
                   prot_func_OC[, c("name", "product")],
                   by.x = "SNP", by.y = "name", all.x = TRUE)
gsub_long$samplenum<-substr(gsub_long$Individual,8,10) %>% as.numeric()
gsub_long<-gsub_long[!(gsub_long$samplenum>76 & gsub_long$samplenum<142),] 
unique(gsub_long$Individual)
# 4. Faceted genomic cline plots
library(dplyr)
indiv_dat <- gsub_long %>%
  group_by(samplenum, HybridIndex, EnvVar) %>%
  summarise(.groups = "drop")  # ensures unique rows per individual
mod <- lm(EnvVar ~ HybridIndex, data = indiv_dat)
summary(mod)
rsq <- summary(mod)$r.squared
pval <- round(summary(mod)$coefficients[2,4], 5)
plot<-ggplot(gsub_long, aes(x = HybridIndex, y = EnvVar, color = factor(Allele))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.7) +
  geom_smooth(aes(group=1), method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ product, scales = "free_y", labeller = label_wrap_gen(width = 27))+
  scale_color_manual(values = c("0" = "blue", "1" = "purple", "2" = "red"),
                     labels = c("AA", "AB", "BB")) +
  labs(x = bquote(bold("Hybrid Index (0 = " * italic("Q. bicolor")*", 1 = "*italic("Q. lyrata")*")")),
       y = "pH", color = "Allele state")+
  theme_minimal(base_family = "mono") +
  theme(strip.text = element_text(size=10, face="bold"),
        axis.title.x = element_text(size = 17, face="bold"),
        axis.title.y = element_text(size = 17,face="bold"),
        legend.position = c(0.95, 0.05),
        legend.title = element_text(size=17, face="bold"),
        legend.justification = c("right", "bottom"))+
  labs(caption = paste0("Linear model: R² = ", round(rsq, 3),
                        ", p = ",pval)
  ) +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 12, face = "bold")
  )
plot  
ggsave("twoloci_white.png",plot=plot, width = 12, height = 8, bg='white')

#plot all env variables in one plot for one protein product
library(tidyr)
gsub_long2<-gsub_long[gsub_long$product=="auxin transporter-like protein 2",]
allext_long <- allext %>%
  pivot_longer(cols = -Individual, names_to = "EnvVarName", values_to = "EnvVar")

gsub_long2 <- gsub_long2 %>%
  left_join(allext_long, by = "Individual", relationship = "many-to-many")
gsub_long2 <- gsub_long2[gsub_long2$EnvVarName %in% 
                           c("cfvo","TRI","slope","bio10","bio15","soc"), ]
library(dplyr)
stats_df <- gsub_long2 %>%
  group_by(EnvVarName, product) %>%
  summarise(
    mod = list(lm(EnvVar.y ~ HybridIndex, data = cur_data())),
    ymax = max(EnvVar.y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    rsq = summary(mod)$r.squared,
    pval = summary(mod)$coefficients[2, 4],
    label = paste0("R² = ", round(rsq, 3), 
                   ",\np = ", signif(pval, 3))
  )

plot<-ggplot(gsub_long2, aes(x = HybridIndex, y = EnvVar.y, color = factor(Allele))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ EnvVarName,
             scales = "free_y",
             labeller = label_wrap_gen(width = 27),
             ncol = 3) +
  scale_color_manual(values = c("0" = "blue", "1" = "purple", "2" = "red"),
                     labels = c("AA", "AB", "BB")) +
  labs(x = "Hybrid Index",
       y = "Environmental variable",
       color = "Allele state for auxin transporter-like protein 2: ") +
  geom_text(
    data = stats_df,
    aes(x = 0.5, y = ymax, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1.5,
    size = 4, family="mono", fontface="bold"
  )+
  theme_minimal(base_family = "mono") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(t = -10)  # move up (negative top margin)
  )+
  theme(strip.text = element_text(size=19, face="bold"),
        axis.title.x = element_text(size = 17, face="bold"),
        axis.title.y = element_text(size = 17,face="bold"),
        legend.title = element_text(size=15),
        theme(plot.caption = element_text(hjust = 0.5, size = 20, face = "bold")
        ))

plot
ggsave("8_20_25_justTN_allENV_RDAselLoci_white.png",plot=plot, width = 12, height = 8, bg='white')





#plot just one plot

library(reshape2)
gsub_long <- reshape2::melt(as.matrix(gsub_matrix))
colnames(gsub_long) <- c("Individual", "SNP", "Allele")

# 2. Add hybrid index and environment data
gsub_long$HybridIndex <- hi$hi[match(gsub_long$Individual, hi$Individual)]
gsub_long$EnvVar <- allext$bio10[match(gsub_long$Individual, allext$Individual)]

# 3. Merge SNP functional annotations (optional, for better facet labels)
gsub_long <- merge(gsub_long,
                   prot_func_OC[, c("name", "product")],
                   by.x = "SNP", by.y = "name", all.x = TRUE)
gsub_long$samplenum<-substr(gsub_long$Individual,8,10) %>% as.numeric()
gsub_long<-gsub_long[!(gsub_long$samplenum>76 & gsub_long$samplenum<142),] 
unique(gsub_long$Individual)
gsub_long<-gsub_long[gsub_long$product=='auxin transporter-like protein 2',]


# 4. Faceted genomic cline plots
library(dplyr)
indiv_dat <- gsub_long %>%
  group_by(samplenum, HybridIndex, EnvVar) %>%
  summarise(.groups = "drop")  # ensures unique rows per individual
mod <- lm(EnvVar ~ HybridIndex, data = indiv_dat)
summary(mod)
rsq <- summary(mod)$r.squared
pval <- round(summary(mod)$coefficients[2,4], 5)
plot<-ggplot(gsub_long, aes(x = HybridIndex, y = EnvVar, color = factor(Allele))) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.7) +
  geom_smooth(aes(group=1), method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ product, scales = "free_y", labeller = label_wrap_gen(width = 27))+
  scale_color_manual(values = c("0" = "blue", "1" = "purple", "2" = "red"),
                     labels = c("AA", "AB", "BB")) +
  labs(x = bquote(bold("Hybrid Index (0 = " * italic("Q. bicolor")*", 1 = "*italic("Q. lyrata")*")")),
       y = "Mean Temp. of Warmest Quarter", color = "Allele state")+
  theme_minimal(base_family = "mono") +
  theme(strip.text = element_text(size=10, face="bold"),
        axis.title.x = element_text(size = 17, face="bold"),
        axis.title.y = element_text(size = 17,face="bold"),
        legend.position = c(0.95, 0.05),
        legend.title = element_text(size=12, face="bold"),
        legend.justification = c("right", "bottom"))+
  labs(caption = paste0("Linear model: R² = ", round(rsq, 3),
                        ", p = ",pval)
  ) +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 12, face = "bold")
  )
plot  
ggsave("10_13_25_bio10_singleplot_white.png",plot=plot, width = 7, height = 6, bg='white')

###########################################################
############# predict changes to ENV on allele states
##########################################################
install.packages('AlleleShift')
library(AlleleShift)
library(BiodiversityR) # also loads vegan
library(poppr) # also loads adegenet
library(ggplot2)
library(ggsci)

install.packages("gganimate")
library(ggforce)
library(dplyr)
library(ggrepel)
library(patchwork)
library(GGally)
library(mgcv)
library(ggmap)
library(gggibbous)
library(gganimate)

library(AlleleShift)
genepo<-gl2genepop(gsub1,"genepo.gen")
colnames(gsub1@tab)[1:10]

fut<-rast("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/modeling_sept2024/102224_models/cropped_futureclimate/ecearth3_2070_585.tif")
preds











############################################
############################################
############################################
############################################

# Choose a candidate locus to plot
locus_name <- hyb_rda.cand[20]
# Extract genotype calls (0 = hom ref, 1 = het, 2 = hom alt)
geno_values <- as.numeric(g2_matrix_imp[,locus_name])

# Assign shapes for allele state
allele_shapes <- rep(21, length(geno_values)) # default
allele_shapes[geno_values == 0] <- 21  # hom ref. circle
allele_shapes[geno_values == 1] <- 22  # het. square
allele_shapes[geno_values == 2] <- 24  # hom alt. triangle

# Extract hybrid index for each individual (replace with your actual vector)
# Example: hybrid_index <- your_hybrid_index_vector
# Must be same length as number of individuals

hybrid_index <- hi  # replace with actual variable
# Color palette for hybrid index
hybrid_colors <- colorRampPalette(c("blue", "white", "red"))(100)
color_by_hi <- hybrid_colors[as.numeric(cut(hybrid_index, breaks=100))]
# Plot with color = hybrid index, shape = allele state
plot(hyb_rda, type="n", scaling=1, xlim=c(-2,6), ylim=c(-3,1),
     main=paste("Hybrid Index & Allele State for", locus_name))
points(hyb_rda, display="sites", pch=allele_shapes, cex=1.2,
       col="black", bg=color_by_hi, scaling=3)
text(hyb_rda, scaling=3, display="bp", col="#0868ac", cex=1)
####

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)

# Example: your SNP list
snp_strings <- hyb_rda.cand
snp_strings <- both$locus
# Split into chromosome and position
snp_df <- data.frame(
  chr = sub("_[0-9]+$", "", snp_strings),                 # drop last "_number"
  pos = as.integer(sub("^.*_", "", snp_strings))          # take only last number after "_"
)
# Check
print(snp_df)
# Function to get protein products (as before)
get_protein_products <- function(gff_file, snp_df) {
  gff_gr <- import(gff_file)
  snp_gr <- GRanges(seqnames = snp_df$chr,
                    ranges = IRanges(start = snp_df$pos, end = snp_df$pos))
  hits <- findOverlaps(snp_gr, gff_gr)
  overlap_df <- data.frame(
    snp_chr   = seqnames(snp_gr)[queryHits(hits)],
    snp_pos   = start(snp_gr)[queryHits(hits)],
    gene_name = mcols(gff_gr)$gene[subjectHits(hits)],
    type      = mcols(gff_gr)$type[subjectHits(hits)],
    product   = mcols(gff_gr)$product[subjectHits(hits)],
    protein_id= mcols(gff_gr)$protein_id[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  return(overlap_df)
}

# Run
gff_path <- "./snpEff/data/genomic.sorted.gff"
protein_info <- get_protein_products(gff_path, snp_df)
write.csv(protein_info,"adaptiveSNPsAndRDA_protein_products.csv")
protein_info$product
protein_info$name<-paste0(protein_info$snp_chr, "_", protein_info$snp_pos)
length(unique(protein_info$name))
(unique(protein_info$gene_name))

# Save
write.csv(protein_info, "protein_products_for_snps.csv", row.names = FALSE)














































edge<-g[g@strata$.!="CORE"]
edge<-gl.compliance.check(edge)
edge<-gl.drop.ind(edge,c("qlyrata","SRR5284357_quly_GA", "SRR5284358_quly_MO", "SRR5632350_quly_FL","SRR5632442_quly_IL", "sample_120", "sample_121", "sample_122", "sample_123", "sample_124", "sample_125",
                         "sample_126"))
edge<-gl.filter.monomorphs(edge)
edge<-gl.filter.callrate(edge, threshold=1)
#saveRDS(edge,"edge_vc3")

##############################################################################
###### start here ############################################################
library(vcfR)
library(dartR)
library(dplyr)

setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline")
edge<-readRDS("edge_vc3")
#perform RDA


#loci<-read.csv("impeded_and_envSelected_SNPs.csv")
#loci_ada<-read.csv("./711_gcom_adaptive_loci_0.70_hybridTN.csv") 
#edge<-gl.keep.loc(edge,loci_ada$gcom_adaptive_0.80_TN.locus)


df1<-cbind(edge@ind.names,as.character(edge@pop)) %>% as.data.frame()

po <- read.csv("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/141_2spec_qs.csv", header=TRUE)
pop<-cbind(po$Individual,po$lon, po$lat, po$pop5) %>% as.data.frame()
pop<-pop[(pop$V1 %in% edge@ind.names),]
po<-cbind(po$Individual,po$lon, po$lat) %>% as.data.frame()
po<-po[(po$V1 %in% edge@ind.names),]

coords <- po[match(edge@ind.names, po$V1), ] %>% na.omit()
colnames(coords)<-c("ID","x","y")
coords$x<-as.numeric(coords$x)
coords$y<-as.numeric(coords$y)

library(terra)
preds<-rast("D:/data/predictors_30s_easternUS.tif")
preds<-preds[[names(preds) != "cec"]]
#preds<-c(preds$bio1,preds$bio10, preds$bio11,preds$bio12,preds$bio13,preds$bio14,preds$bio15,preds$bio16,preds$bio17,preds$bio18,preds$bio19,preds$bio2,preds$bio3,preds$bio4,preds$bio5,preds$bio6,preds$bio7,preds$bio8,preds$bio9)  
names(preds)

locs<-vect(coords, crs=crs(preds), c('x','y'))

extract_locs<-extract(preds, locs, xy=TRUE)
allex<-cbind(coords$ID, extract_locs)
allex<-na.omit(allex)
colnames(allex)
slx<-allex[match(po$V1, allex$`coords$ID`), ]
slx$pop<-pop$V4
slx<-slx[order(slx$pop), ]
barplot(names.arg=slx$pop, height=slx$TRI, las=2)
#higher slope generally, #higher slope generally, #higher slope generally, 
#allext<-allex[3:21]
allext<-allex[3:32]

#remove correlated variables
library("usdm")
vstep<-vifstep(allext, th=10)
selpreds<-preds[[vstep@results$Variables]]
extract_locs2<-extract(selpreds, locs, xy=TRUE)
allex2<-cbind(coords$ID, extract_locs2)
allex2<-na.omit(allex2)
colnames(allex2)
#allext2<-allex2[3:8]
allext2<-allex2[3:13]

library(vegan)
library(ggplot2)
library(grid)
g2_matrix<-as.matrix(edge)

hyb_rda <- rda(g2_matrix ~ .,data=allext2, scale=T)
RsquareAdj(hyb_rda)
summary(hyb_rda)$concont
plot(hyb_rda, scaling=3)  ## d.fault is axes 1 and 2

ind_scores <- scores(hyb_rda, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
ind_scores_df$ID <- edge@ind.names
ind_scores_df$pop <- edge@pop

est2<-read.csv("est2_5.csv")
est<-est2[(est2$hi.INDLABEL %in% edge@ind.names),]
hi <- est$hi.h_posterior_mode[match(edge@ind.names, est$hi.INDLABEL)]
ind_scores_df$HybridIndex <- hi

env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 10  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
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
       x = paste0("RDA1"," ", percexpl_rda1,"%"), y = paste0("RDA2"," ", percexpl_rda2,"%"), color = "Hybrid Index") +
  theme_minimal()
plot
#ggsave("715_RDA Plot with Individuals and Environmental Vectors label.png", plot=plot, width = 10, height = 10, bg='grey')
#ggsave("715_RDA Plot with Individuals and Environmental Vectors WHITE.png", plot=plot, width = 8, height = 8, bg='white')

# SNP matrix (response) and group factor (e.g., Hybrid vs Non-hybrid)
allext2$hi<-hi
adonis_result <- adonis2(g2_matrix ~ hi + ., data=allext2,method = "euclidean")
print(adonis_result)
cor.test(allext2$hi, allext2$bio11, method = "pearson")


snp_loadings <- scores(hyb_rda, display = "species", scaling = 3)[, 1:2]
top_snps <- order(abs(snp_loadings[,1]), decreasing = TRUE)[1:20]
head(snp_loadings[top_snps, ])

screeplot(hyb_rda)
load.rda <- summary(hyb_rda)$species[,1:2]
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) 
cand2 <- outliers(load.rda[,2], 3) 

hyb_rda.cand <- c(names(cand1), names(cand2)) ## j.st the names of the candidates
length(hyb_rda.cand)
hyb_rda.cand<-hyb_rda.cand[!duplicated(hyb_rda.cand)] ## 7.duplicate detections (detected on multiple RDA axes)

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(g2_matrix) %in% hyb_rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(g2_matrix) %in% hyb_rda.cand, 'red', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
#png("rda_snps_selected_plot.png",height=1000,width=1000)
#plot(hyb_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="hyb_RDA, axes 1 and 2")
#points(hyb_rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
#points(hyb_rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
#text(hyb_rda, scaling=3, display="bp", col="#0868ac", cex=1)
#dev.off()


hyb_rda.cand
gc<-readRDS("711_gc.rds")
hist(gc$gc$S1.prop_1)
snp_genotypes <- g2_matrix[, hyb_rda.cand]
est2<-readRDS("711_est2.rds")
esthi<-est2$hi 
hi<-est2$hi$beta_mean
length(hi)
hist(hi, breaks=100, labels=esthi$INDLABEL)
esthi
hist_data <- data.frame(ID = esthi$INDLABEL, HybridIndex = hi, pop=esthi$POPID)
levels(hist_data$ID)<-c(levels(hist_data$ID),"SRR26194539-Qlyrata_ref","SRR26194547-Qbicolor_ref")
hist_data[6,1]<-"SRR26194539-Qlyrata_ref"	
hist_data[5,1]<-"SRR26194547-Qbicolor_ref"
hist_data <- hist_data[order(hist_data$HybridIndex), ]
hist_data$ID <- factor(hist_data$ID, levels = hist_data$ID) 
hist_data
hIndex<-ggplot(hist_data, aes(x = ID, y = HybridIndex)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Individual", y = "Hybrid Index")
#ggsave("hybrid_indeces_allSamples.png",hIndex, width = 16, height = 8, bg='white')

gc<-readRDS("711_gc.rds")
gcc<-gc$gc %>% as.data.frame()
length(gcc$locus)
hist(gcc$exp_mean_log_v, breaks=100)
dim(gc)

length(hyb_rda.cand)
#2057 SNPs environmentally selected according to RDA
#1055 at sd3
gcc_sub<-gcc[gcc$locus %in% hyb_rda.cand,]
dim(gcc_sub)
#608 of these SNPs have greater than 0.7 min.diff (parental diff in gghybrid)
hist(hist(gcc$exp_mean_log_v, breaks=1000))
hist(gcc_sub$exp_mean_log_v, breaks=1000)
hist(gcc$invlogit_mean_logit_centre, breaks=20)
hist(gcc_sub$invlogit_mean_logit_centre, breaks=20)

g2_matrix_sub<-g2_matrix[,colnames(g2_matrix) %in% gcc$locus]
dim(g2_matrix_sub)
dim(allext2)






loci<-read.csv("impeded_and_envSelected_SNPs.csv")
loci_ada<-read.csv("./711_gcom_adaptive_loci_0.70_hybridTN.csv") 
edge_sel<-gl.keep.loc(edge,loci_ada$gcom_adaptive_0.80_TN.locus)




# Run RDA with bio11 as the only predictor
#check VIF
test<-c(preds$bio5,preds$bio11,preds$bio7,preds$bio2,preds$slope,preds$nit)  
test.locs<-vect(coords, crs=crs(test), c('x','y'))

test.extract_locs<-extract(test, locs, xy=TRUE)
test.allex<-cbind(coords$ID, test.extract_locs)
test.allex<-na.omit(test.allex)
colnames(test.allex)
test.allext<-test.allex[3:8]
vstep<-vif(test.allext)
vstep



#######################################################################################################################
############# START HERE ###########################################################################################
library(vcfR)
library(dartR)
library(dplyr)
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline")
edge<-readRDS("edge_vc3")
gl.filter.allna(edge)
gl.filter.monomorphs(edge)
gl.filter.callrate(edge)
po <- read.csv("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/141_2spec_qs.csv", header=TRUE)
pop<-cbind(po$Individual,po$lon, po$lat, po$pop5) %>% as.data.frame()
pop<-pop[(pop$V1 %in% edge@ind.names),]
po<-cbind(po$Individual,po$lon, po$lat) %>% as.data.frame()
po<-po[(po$V1 %in% edge@ind.names),]
coords <- po[match(edge@ind.names, po$V1), ] %>% na.omit()
colnames(coords)<-c("ID","x","y")
coords$x<-as.numeric(coords$x)
coords$y<-as.numeric(coords$y)
library(vegan)
library(usdm)
library(terra)
library(ggplot2)
library(grid)
edge@ind.names
g2_matrix<-as.matrix(edge)
preds<-rast("D:/data/predictors_30s_easternUS.tif")
preds<-preds[[names(preds) != "cec"]]
locs<-vect(coords, crs=crs(preds), c('x','y'))
extract_locs<-extract(preds, locs, xy=TRUE)
allex<-cbind(coords$ID, extract_locs)
allex<-na.omit(allex)
colnames(allex)
allext<-allex[3:32]
vstep<-vifstep(allext, th=10)

#hyb_rda_bio11 <- rda(g2_matrix ~ ., data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio11+bio7+slope+nit, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio11+bio7+bio5+bio2+slope+nit, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio15+bio16+bio6+bio5, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio15+bio16+bio6+bio5+slope, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio15+bio16+bio6+bio5+slope+nit, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio13+bio15+bio3+bio7+bio8+slope+ph+ocd+TPI+TRI, data = allext)

hyb_rda_bio11 <- rda(g2_matrix ~ bio3+bio15+bio6+bio7+nit+slope+ph+ocd+cfvo+TPI+TRI, data = allext)

#hyb_rda_bio11 <- rda(g2_matrix ~ bio15+bio16, data = allext)
#hyb_rda_bio11 <- rda(g2_matrix ~ bio5 +bio6, data = allext)
# View summary to inspect loadings
summary(hyb_rda_bio11)
RsquareAdj(hyb_rda_bio11)
screeplot(hyb_rda_bio11)
anov<-anova.cca(hyb_rda_bio11, by="axis")
saveRDS(anov, "8_14_25_anova_hybRDA_vifstep10")
#seven axes are signficiant, only keep four: signficant p value and F value greater than 2

# Extract SNP (species) loadings on RDA1
snp_loadings_bio11 <- scores(hyb_rda_bio11, display = "species", scaling = 3, choices=1:4)
# Optional: histogram of loadings to visualize distribution
# Choose threshold z (e.g., 3 standard deviations)
screeplot(hyb_rda_bio11)
load.rda <- summary(hyb_rda_bio11)$species[,1:2]

# Optional: histogram of loadings to visualize distribution for each axis
par(mfrow = c(2, 2))
for (i in 1:4) {
  hist(snp_loadings_bio11[, i],
       main = paste("SNP loadings on RDA", i),
       xlab = "Loading")
}
par(mfrow = c(1, 1))
# Identify outliers for each of the first 6 RDA axes
# Identify outliers based on loading extremes
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]  # return names of outlier loci
}

cand_list <- list()
for (i in 1:4) {
  cand_list[[i]] <- outliers(snp_loadings_bio11[, i], 3.5)
}
# Combine candidate loci from all 4 axes and remove duplicates
hyb_rda.cand <- unique(unlist(lapply(cand_list, names)))
# Check number of candidate loci
length(hyb_rda.cand)
#6610 candidates at SD 3, 3586 at SD 3.5
# Set up the color scheme for plotting:
dim(g2_matrix)
bgcol  <- ifelse(colnames(g2_matrix) %in% hyb_rda.cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(g2_matrix) %in% hyb_rda.cand, 'red', '#00000000')


## a.es 1 & 2 - zooming in to just the SNPs here...
plot(hyb_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="hyb_RDA, axes 1 and 2")
points(hyb_rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(hyb_rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(hyb_rda, scaling=3, display="bp", col="#0868ac", cex=1)

#get just the negative RDA1 values
#snp_loadings_bio11 <- as.numeric(snp_loadings_bio11)
#names(snp_loadings_bio11) <- rownames(scores(hyb_rda_bio11, display = "species", scaling = 3))
#str(snp_loadings_bio11)
#negative_rda1_snps <- snp_loadings_bio11[snp_loadings_bio11 < 0]
#length(negative_rda1_snps)
#length(snp_loadings_bio11)
#negative_snp_names <- names(negative_rda1_snps)
#cand_neg_names<-candidate_snps_bio11[candidate_snps_bio11 %in% negative_snp_names]

#hi_locNEG<-gccc[gccc$locus %in% cand_neg_names, ]
#length(unique(gc$locus))
#length(unique(cand_neg_names))
#length(hi_locNEG$locus)
#hist(hi_locNEG$invlogit_mean_logit_centre, breaks=100)
#length(hi_locNEG)
#length(gccc$locus)


ind_scores <- scores(hyb_rda_bio11, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
ind_scores_df$ID <- edge@ind.names
ind_scores_df$pop <- edge@pop

est2<-read.csv("est2_5.csv")
est<-est2[(est2$hi.INDLABEL %in% edge@ind.names),]
hi <- est$hi.h_posterior_mode[match(edge@ind.names, est$hi.INDLABEL)]
ind_scores_df$HybridIndex <- hi

env_arrows <- scores(hyb_rda_bio11, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 10  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda_bio11)
perc_cons<-round((hyb_rda_bio11$CCA$tot.chi/hyb_rda_bio11$tot.chi)*100, 2)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  # Add custom text in bottom right
  #add sample names
  #geom_text(data = ind_scores_df,
  #          aes(x = RDA1, y = RDA2, label = ID),
  #          color = "black", vjust = 1, size = 1) +
  labs(title = "",
       x = paste0("RDA1"," ", percexpl_rda1,"% (constrained) [total constrained = ",perc_cons,"%]"), y = paste0("RDA2"," ", percexpl_rda2,"% (constrained)"), color = "Hybrid Index") +
  theme_minimal()
plot
ggsave("8_14_25_RDAplot_vifstep10_grey.png",plot=plot, width = 8, height = 8, bg='grey')
ggsave("8_14_25_RDAplot_vifstep10_white.png",plot=plot, width = 8, height = 8, bg='white')

plot <- ggplot() +
  geom_point(data = ind_scores_df,
             aes(x = RDA1, y = RDA2, color = HybridIndex, shape = pop),
             size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  scale_shape_manual(values = 0:25) +  # Use up to 26 different shapes if needed
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = 1, size = 5) +
  labs(title = "",
       x = paste0("RDA1 ", percexpl_rda1, "%"),
       y = paste0("RDA2 ", percexpl_rda2, "%"),
       color = "Hybrid Index",
       shape = "Population") +
  theme_minimal()
plot




gc<-readRDS("711_gc.rds")
gcc<-gc$gc %>% as.data.frame()
length(gcc$locus)
hist(gcc$exp_mean_log_v, breaks=100)
hist(gcc$mean_log_v, breaks=100)
dim(gcc)
length(hyb_rda.cand)
#2449
gcc_sub<-gcc[gcc$locus %in% hyb_rda.cand,]
dim(gcc_sub)
#474 of these SNPs have greater than 0.7 min.diff (parental diff in gghybrid)
hist(gcc$mean_log_v, breaks=100, xlim=c(-2,2), xlab="mean log(v)")
hist(gcc_sub$mean_log_v, breaks=100,xlim=c(-2,2),xlab="mean log(v)")

hist(gcc$invlogit_mean_logit_centre, breaks=100)
hist(gcc_sub$invlogit_mean_logit_centre, breaks=100)
g2_matrix_sub<-g2_matrix[,colnames(g2_matrix) %in% gcc$locus]
dim(g2_matrix_sub)

dim(allext2)
saveRDS(hyb_rda.cand, "8_14_25_hyb_rda_cand_3point5SD")
ada_loci<-read.csv("./adaptive_liberal_loci_0.70_hybridTN.csv")
imp_loci<-read.csv("./impeded_liberal_loci_0.70_hybridTN.csv")
length(ada_loci$gcom_ada_0.80_TN_liberal.locus)
#109
ada_loci$gcom_ada_0.80_TN_liberal.locus
ada_both<-ada_loci[(ada_loci$gcom_ada_0.80_TN_liberal.locus %in% hyb_rda.cand),]
length(ada_both)
#0
length(imp_loci$gcom_impeded_0.80_TN.locus)
#531
imp_loci$gcom_impeded_0.80_TN.locus
imp_both<-imp_loci[(imp_loci$gcom_impeded_0.80_TN.locus %in% hyb_rda.cand),]
length(imp_both)
#85


dim(gcc_sub)
length(imp_both)
length(ada_both)
length(hyb_rda.cand)

ada_both_ann<-ada.ann[ada.ann$name %in% ada_both,]
length(unique(ada_both_ann$name))
protfunc<-read.csv("./snpEff/Extracted_adaptive_protein_functions.csv")
protfunc$name<-paste0(protfunc$snp_chr, "_", protfunc$snp_pos)
prot_func_ada<-protfunc[protfunc$name %in% ada_both,]
length(unique(prot_func_ada$gene_name))
unique(prot_func_ada$product)
write.csv(unique(prot_func_imp$product),"8_14_25_vifstep10_RDA3point5SD_AdaptiveAllproteinfunctions.csv")


fromOvercup<-gcc_sub[gcc_sub$invlogit_mean_logit_centre<0.4 & gcc_sub$locus %in% ada_both, ] 
fromOvercup<-fromOvercup$locus
ada.ann<-read.delim("./snpEff/Extracted_adaptive_loci_annotations_long.csv", header=TRUE, sep=",")
ada.ann$name <- paste0(ada.ann$CHROM, "_", ada.ann$POS)
ada_both_ann<-ada.ann[ada.ann$name %in% fromOvercup,]
length(unique(ada_both_ann$name))
protfunc<-read.csv("./snpEff/Extracted_adaptive_protein_functions.csv")
protfunc$name<-paste0(protfunc$snp_chr, "_", protfunc$snp_pos)
prot_func_ada<-protfunc[protfunc$name %in% fromOvercup,]
length(prot_func_ada$snp_chr)
length(unique(prot_func_ada$gene_name))
unique(prot_func_ada$product)
write.csv(unique(prot_func_ada$product),"8_14_25_vifstep10_RDA3point5SD_AdaptiveGeneFuncsFromOvercup.csv")

fromSWO<-gcc_sub[gcc_sub$invlogit_mean_logit_centre>0.6 & gcc_sub$locus %in% ada_both, ] 
fromSWO<-fromSWO$locus
ada.ann<-read.delim("./snpEff/Extracted_adaptive_loci_annotations_long.csv", header=TRUE, sep=",")
ada.ann$name <- paste0(ada.ann$CHROM, "_", ada.ann$POS)
ada_both_ann<-ada.ann[ada.ann$name %in% fromSWO,]
length(unique(ada_both_ann$name))
protfunc<-read.csv("./snpEff/Extracted_adaptive_protein_functions.csv")
protfunc$name<-paste0(protfunc$snp_chr, "_", protfunc$snp_pos)
prot_func_ada<-protfunc[protfunc$name %in% fromSWO,]
length(prot_func_ada$snp_chr)
length(unique(prot_func_ada$gene_name))
unique(prot_func_ada$product)
write.csv(unique(prot_func_ada$product),"8_14_25_vifstep10_RDA3point5SD_AdaptiveGeneFuncsFromSWO.csv")

#write.csv(imp_both, "impeded_and_envSelected_SNPs_4variBIOCLIM.csv")
imp.ann<-read.delim("./snpEff/Extracted_impeded_loci_annotations_long.csv", header=TRUE, sep=",")
imp.ann$name <- paste0(imp.ann$CHROM, "_", imp.ann$POS)
imp_both_ann<-imp.ann[imp.ann$name %in% imp_both,]
length(unique(imp_both_ann$name))
protfunc<-read.csv("./snpEff/Extracted_impeded_protein_functions.csv")
protfunc$name<-paste0(protfunc$snp_chr, "_", protfunc$snp_pos)
prot_func_imp<-protfunc[protfunc$name %in% imp_both,]
length(unique(prot_func_imp$gene_name))
unique(prot_func_imp$product)
write.csv(unique(prot_func_imp$product),"8_11_25_impeded_andenvselected_vifstep10_proteinfunctions.csv")

#write.csv(ada_both, "adaptive_and_envSelected_SNPs_4variBIOCLIM.csv")
ada.ann<-read.delim("./snpEff/Extracted_adaptive_loci_annotations_long.csv", header=TRUE, sep=",")
protfunc<-read.csv("./snpEff/Extracted_adaptive_protein_functions.csv")
length(protfunc$snp_chr)
unique(protfunc$product)

ada.ann$name <- paste0(ada.ann$CHROM, "_", ada.ann$POS)
ada_both_ann<-ada.ann[ada.ann$name %in% ada_both,]
length(unique(ada_both_ann$name))

protfunc$name<-paste0(protfunc$snp_chr, "_", protfunc$snp_pos)
prot_func_ada<-protfunc[protfunc$name %in% ada_both,]
length(prot_func_ada$snp_chr)
length(unique(prot_func_ada$gene_name))
unique(prot_func_ada$product)
write.csv(unique(prot_func_ada$product),"8_11_25_adaptive_andenvselected_vifstep10_proteinfunctions.csv")































# g2_matrix: samples x SNPs
# cline_centers: named vector (names = SNPs), values = cline centers (e.g., ancestry value at midpoint)

# Example: match SNPs in g2_matrix with cline center values
cline_centers <- gc$invlogit_mean_logit_centre  # or similar
names(cline_centers) <- gc$locus

# Ensure cline centers match g2_matrix columns
cline_centers <- cline_centers[colnames(g2_matrix)]

# Confirm alignment
all(names(cline_centers) == colnames(g2_matrix))  # should be TRUE
length(cline_centers)








####################################################
###########       GGHYBRIDS           ##############
####################################################
library(vcfR)
library(dartR)
library(SNPfiltR)

setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline")
list.files()
vc<-read.vcfR("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline_ready_snps.recode.vcf")
g<-vcfR2genlight(vc)
names<-read.csv("bamlist.csv", header=FALSE)
gl.filter.callrate(g, method='ind')
g@ind.names<-names$V1
g@pop<-names$V4 %>% as.factor()
g@strata<-names[,3:4] %>% as.data.frame()
g<-gl.compliance.check(g)
g<-g[g@ind.names!="qbicolor"]

edge<-g[g@strata$V3!="CORE"]
edge@pop
edge<-gl.filter.callrate(edge, method='ind',threshold=0.8)
edge<-gl.filter.callrate(edge, threshold=1)
edge<-gl.filter.monomorphs(edge)
edge@pop
#67 genotypes, 5769 SNPs 0% missing data
gl2structure(g, outpath=getwd(), outfile = "gl3.str", addcolumns = g@pop)

#gghybrids
#install.packages("devtools"); devtools::install_github("ribailey/gghybrid")
library(gghybrid)
dat2<-gghybrid::read.data("gl3.str", precol.headers = 0 , nprecol = 2, markername.dup = 0, INDLABEL = 1, MISSINGVAL= -9, NUMINDS= 135, NUMLOCI=43339,ONEROW=0,POPID=1)
length(dat2$loci)

prepdata2<-gghybrid::data.prep(data=dat2$data, loci=dat2$loci, sourceAbsent = FALSE, S0="P0", S1="P1", alleles=dat2$alleles, precols=dat2$precols, return.genotype.table = TRUE, return.locus.table = TRUE, INDLABEL.name="INDLABEL",POPID.name = "POPID", min.diff=0.5)
str(prepdata2)
length(prepdata2$locus.data$locus)
#5532 loci left!!

est2<-esth(prepdata2$data.prep, read.data.precols = dat2$precols, nitt=10000, burnin=5000, init.var=0.002)
#write.csv(est2, "est2_4.csv")
write.csv(est2,"est2_5.csv")
hist(est2$hi$h_posterior_mode)

#saved est2.csv based on running with all core and edge samples and both ref samples but with CORE-IN1 as the P0 and qlyrata as P1. this was with ~9500 SNPs, callrate of 1, filtered a few individuals (down to 0.9 i think?). shows that the ref sample is 0.39 quercus lyrata. 
hist(est2$hi$h_posterior_mode, breaks=100)

table(prepdata2$geno.data$Source)
prepdata2$geno.data
table(prepdata2$geno.data$POPID, prepdata2$geno.data$Source)

cline<-ggcline(data.prep.object=prepdata2$data.prep, esth.object=est2, nitt=10000, burnin=5000, include.Source=TRUE, return.likmeans = TRUE, plot.test.subject = c("NC_044904.1_10198924","NC_044904.1_1023825"), plot.col=c("orange","cyan"), plot.ylim=c(-3,5), plot.pch.v.centre = c(1, 3), read.data.precols = dat2$precols)
getwd()
#saveRDS(cline,"ggcline_object_062325")
#saveRDS(est2,"esth_object_062325")
#saveRDS(prepdata2,"prepdata_object_062325")

getwd()
write.table(cline$gc,"gc_out",quote = F,row.names = F)
cline<-readRDS("ggcline_object_062325")
cline$gc
setdiff(c("NW_022154797.1_53141","NC_044906.1_29230604","NC_044911.1_31159311","NC_044906.1_7813205","NC_044906.1_55367622"), 
        cline$gc$locus)
cline$gc$locus
#png(filename="CLINE_OUT.png",width=1200,height=500)
plot_clinecurve(ggcline.object=cline$gc,cline.locus=c("NW_022154797.1_53141","NC_044906.1_29230604","NC_044911.1_31159311","NC_044906.1_7813205","NC_044906.1_55367622"),locus.column="locus",cline.col=c("orange","blue","green","red","magenta"),null.line.locus=c("NW_022154797.1_53141","NC_044906.1_29230604","NC_044911.1_31159311","NC_044906.1_7813205","NC_044906.1_55367622"),null.line.col=c("orange","blue","green","red","magenta"),cline.centre.line=c("NW_022154797.1_53141","NC_044906.1_29230604","NC_044911.1_31159311","NC_044906.1_7813205","NC_044906.1_55367622"),cline.centre.col=c("orange","blue","green","red","magenta")
)
plot_clinecurve(ggcline.object=cline$gc,cline.locus=cline$gc$locus,locus.column="locus")
)

#Add a title and axis labels:
#title(main = "NW_022154797.1_53141, NC_044906.1_29230604, NC_044911.1_31159311, NC_044906.1_7813205, NC_044906.1_55367622", xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1.5,cex.lab=1.5)
#dev.off()
length(cline$locus)
subset<-cline[cline$exp_mean_log_v>1 & cline$v_pvalue<0.05,]
length(subset$locus)

hyb_rda.cand
#get candidate snps in both cline and RDA analysis
both<-subset[(subset$locus %in% hyb_rda.cand),]
length(both$locus)
both_loc<-both$locus
both_loc<-gsub(".1_", ".1: ", both_loc)
lapply(both_loc, write, "06102025_rda_and_gghybrids_loci_sd2point5.txt", append=TRUE, ncolumns=1)

gc_out<-read.table("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline/gc_out", header=TRUE)
sel<-gc_out[gc_out$exp_mean_log_v>1 & gc_out$v_pvalue<0.05,]
dim(sel)

length(sel$locus)/length(gc_out$locus)

inds <- indNames(edge)
meta_ordered <- est2[match(inds, est2$hi.INDLABEL), ]
sorted_meta <- meta_ordered[order(meta_ordered$hi.h_posterior_mode), ]
sort_order <- match(sorted_meta$hi.INDLABEL, inds)
edge_sorted <- edge[sort_order, ]
edge_sorted@ind.names

sub<-gl.keep.loc(edge_sorted, both$locus)
sub@ind.names
"C:\Users\jpark107\genetics\plink-1.07-x86_64.zip"
setwd("C:/Users/jpark107/genetics/")
gl2vcf(x=sub,plink_path="C:/Users/jpark107/Downloads/plink_win64_20241022/" , outfile="GGHYBRID_LOCI_20", outpath=getwd())
sub2<-gl.keep.loc(edge_sorted, subset$locus)
sub2@ind.names
gl2vcf(x=sub2,plink_path="C:/Users/jpark107/Downloads/plink_win64_20241022/" , outfile="GGHYBRID_LOCI_50", outpath=getwd())

est2


#sub<-gl.keep.loc(ed, sel$locus)
#sub@ind.names
"C:\Users\jpark107\genetics\plink-1.07-x86_64.zip"
#gl2vcf(x=sub,plink_path="C:/Users/jpark107/Downloads/plink_win64_20241022/" , outfile="GGHYBRID_LOCI_50", outpath=getwd())







install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
# Connect to Ensembl Plants and the Quercus lobata dataset
mart <- useMart(biomart = "plants_mart",
                dataset = "qlobata_eg_gene",
                host = "https://plants.ensembl.org")
# Paste your gene list here
genes <- c("QL01p023136", "QL01p023736", "QL01p023741", "QL01p025468", "QL01p025993",
           "QL01p026245", "QL01p026257", "QL01p026260", "QL01p026685", "QL01p031287",
           "QL01p031553", "QL01p032655", "QL01p032659", "QL01p032663", "QL02p003969",
           "QL02p003978", "QL02p009188", "QL02p009195", "QL02p017067", "QL02p029261",
           "QL02p029270", "QL02p031474", "QL02p036925", "QL02p036937", "QL02p042449",
           "QL02p042453", "QL02p042456", "QL02p051635", "QL02p057252", "QL02p088845",
           "QL02p104378", "QL03p007813", "QL03p035915", "QL03p047190", "QL03p054439",
           "QL03p055365", "QL03p055373", "QL04p056643", "QL04p080658", "QL05p020644",
           "QL05p020652", "QL05p025724", "QL05p026667", "QL05p043168", "QL05p047742",
           "QL05p047746", "QL06p011235", "QL06p016165", "QL08p000831", "QL08p031153",
           "QL08p031156", "QL08p031162", "QL08p042224", "QL08p061123", "QL08p061129",
           "QL08p064889", "QL08p064899", "QL11p016263", "QL12p012367", "QL93p0093_0048",
           "QL93p0093_0054")
# Query gene annotations (add more attributes if desired)
results <- getBM(attributes = c("ensembl_gene_id", 
                                "external_gene_name", 
                                "description", 
                                "go_id", 
                                "name_1006"),
                 filters = "ensembl_gene_id",
                 values = genes,
                 mart = mart)
# Preview the results
head(results, 20)
# Save to CSV
write.csv(results, "quercus_lobata_gene_functions.csv", row.names = FALSE)





# Define color vector with names
locus_colors <- c("orange", "blue", "green", "purple", "red", "magenta","yellow","darkgreen","brown")
names(locus_colors) <- both$locus

# Plot with named color vectors
plot_clinecurve(
  ggcline.object = cline$gc,
  cline.locus = both$locus,
  locus.column = "locus",
  cline.col = locus_colors,
  null.line.locus = both$locus,
  null.line.col = locus_colors,
  cline.centre.line = both$locus,
  cline.centre.col = locus_colors
)
legend(
  "topright",                        # Position (or use "bottomleft", etc.)
  legend = names(locus_colors),     # Locus names
  col = locus_colors,               # Corresponding colors
  lty = 1,                           # Line type (match plot)
  lwd = 2,                           # Line width
  title = "Locus",                  # Legend title
  cex = 0.8                         # Text size
)
library('viridis')
color_palette <- viridis(length(subset$locus), option='turbo')

# Plot with named color vectors
plot_clinecurve(
  ggcline.object = cline$gc,
  cline.locus = subset$locus,
  locus.column = "locus",
  cline.col = color_palette,
  null.line.locus = subset$locus,
  null.line.col = color_palette,
  cline.centre.line = subset$locus,
  cline.centre.col = color_palette
)
legend(
  "topright",                        # Position (or use "bottomleft", etc.)
  legend = subset$locus,     # Locus names
  col = color_palette,               # Corresponding colors
  lty = 1,                           # Line type (match plot)
  lwd = 2,                           # Line width
  title = "Locus",                  # Legend title
  cex = 0.8                         # Text size
)
















######################################################################################
######################################################################################


#library(terra)
#output_folder<-"D:/data/wc2.1_30s_bio/raw_rasters_touse/QUBI_and_QULY/150km/ready/done"
#files <- list.files(output_folder, pattern = "\\.tif$", full.names = TRUE)
#preds<-rast(files)

#output_folder<-"D:/data/wc2.1_30s_bio/raw_rasters_touse/"
#files <- (list.files(output_folder, pattern = ".tif$", full.names = TRUE))[3:21]
#ds<-rast(files)

#ds<-crop(ds,preds, mask=TRUE)
#preds<-c(ds,preds)
#na_mask <- app(preds, fun = function(x) any(is.na(x)))
#r_masked <- mask(preds, na_mask, maskvalues=1, updatevalue=NA)
#summary(r_masked)
#names(preds)[21]<-"slope"
#names(preds)[20]<-"elevation"
#names(preds)[1:19]<-c("bio1","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9")

preds<-rast("D:/data/predictors_30s_easternUS.tif")
preds<-preds[[names(preds) != "cec"]]

setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/modeling_sept2024/Qlyrata")
qubi<-read.csv("D:/sdm stuff 2025/data_wgs84.csv")
qubi<-qubi[qubi$pa==1,]
qubi<-vect(qubi, crs='EPSG:4326', c('x','y'))

quly<-read.csv("Qlyrata_occurences.csv")
colnames(quly)
quly<-vect(quly, crs='EPSG:4326', c('decimalLongitude','decimalLatitude'))

hybr<-read.csv("C:/Users/jpark107/OneDrive - University of Tennessee/134_coords.csv",header=TRUE)
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline/")
est2<-read.csv("est2_5.csv")
est<-est2[(est2$hi.INDLABEL %in% hybr$ID),]
hi <- est$hi.h_posterior_mode[match(hybr$ID, est$hi.INDLABEL)]
hybr$hi<-hi
#hybr<-hybr[hybr$SWO_prop<0.9,]
hybs<-vect(hybr, crs='EPSG:4326', c('LON','LAT'))

exhyb<-extract(preds, hybs, xy=TRUE)
exquly<-extract(preds, quly, xy=TRUE)
exqubi<-extract(preds, qubi, xy=TRUE)
exquly$spec<-"Q. lyrata"
exqubi$spec<-"Q. bicolor"
exhyb$spec<-"samples"
exquly$prop<-1
exqubi$prop<-0
exhyb$prop<-hybr$hi

allex<-rbind(exqubi, exquly, exhyb)
allex<-na.omit(allex)
allext<-allex[2:31]

#remove correlated variables
library("usdm")
vstep<-vifstep(allext, th=5)
selpreds<-preds[[vstep@results$Variables]]
exhyb<-extract(selpreds, hybs, xy=TRUE)
exquly<-extract(selpreds, quly, xy=TRUE)
exqubi<-extract(selpreds, qubi, xy=TRUE)
exquly$spec<-"Q. lyrata"
exqubi$spec<-"Q. bicolor"
exhyb$spec<-"samples"
exquly$prop<-1
exqubi$prop<-0
exhyb$prop<-hybr$hi
allex<-rbind(exqubi, exquly, exhyb)
allex<-na.omit(allex)
colnames(allex)
allext<-allex[2:16]
library(vegan)
library(ggplot2)
library(grid)
pred.pca <- rda(allext, scale = TRUE)
pca_scores <- scores(pred.pca, display = "sites", scaling = 1)
pca_df <- as.data.frame(pca_scores)
pca_df$spec <- allex$spec %>% as.factor() # your species factor
pca_df$prop <- allex$prop
loadings <- scores(pred.pca, display = "species", scaling = 1)
loadings_df <- as.data.frame(loadings)
loadings_df$var <- rownames(loadings_df)
# Calculate percent variance explained by each PC
eig_vals <- pred.pca$CA$eig
var_exp <- eig_vals / sum(eig_vals) * 100

# Create axis labels with % variance
x_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
y_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")
# Your ggplot with updated axis labels
envpca<-ggplot(pca_df, aes(x = PC1, y = PC2, color = prop, shape = spec)) +
  geom_point(aes(size = spec), alpha = 0.5) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = PC1 * 0.15, yend = PC2 * 0.15),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings_df,
            aes(x = PC1 * 0.12, y = PC2 * 0.12, label = var),
            size = 4, hjust = 0.5, color = "black", inherit.aes = FALSE) +
  labs(title = "",
       x = x_lab,  # updated here
       y = y_lab,  # updated here
       color = "Hybrid Index",
       shape = "Species",
       size = "Species") +
  theme_minimal() +
  scale_shape_manual(values = c("Q. bicolor" = 0, "Q. lyrata" = 1, "samples" = 17)) +
  scale_size_manual(values = c("Q. bicolor" = 1, "Q. lyrata" = 1, "samples" = 5)) +
  scale_color_gradient2(high = "yellow", mid="#D4B65B", low = "purple4", midpoint=.5)
envpca
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/cline/")
ggsave("ENV PCA Qbicolor adn Qlyrata_colorsreverseHI.png", plot=envpca, width = 8, height = 8, bg='grey')


hyb_rda <- rda(allext ~ allex$spec,data=extracts, scale=T)
RsquareAdj(hyb_rda)
summary(hyb_rda)$concont
plot(hyb_rda, scaling=3)  ## d.fault is axes 1 and 2

ind_scores <- scores(hyb_rda, display = "sites", scaling = 3)
ind_scores_df <- as.data.frame(ind_scores[, 1:2])
ind_scores_df$ID <- ed@ind.names
hi<-est2$hi.h_posterior_mode[2:67]
ind_scores_df$HybridIndex <- hi

env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# Optional: rescale arrows to fit the RDA plot scale
arrow_multiplier <- 10  # Adjust if arrows are too big/small
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier

# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = HybridIndex), size = 3) +
  scale_color_gradient(low = "purple4", high = "yellow") +
  # Add environmental vectors as arrows
  geom_segment(data = env_arrows_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Add labels to environmental vectors
  geom_text(data = env_arrows_df,
            aes(x = RDA1, y = RDA2, label = Variable),
            color = "black", vjust = -0.5) +
  
  labs(title = "RDA Plot with Individuals and Environmental Vectors",
       x = "RDA1", y = "RDA2", color = "Hybrid Index") +
  theme_minimal()





################################################
#LEA and Tess
#genome scan for selection
library(LEA)
g@ind.names
gl2geno(g, outfile= "sub_geno",outpath=getwd()) 
popus<-g@pop
#detect admixture
project = snmf("sub_geno.geno", K = 2, entropy = TRUE, repetitions = 100, seed=42,project = "new", iterations=1000000)
project = load.snmfProject("sub_geno.snmfProject")
plot(project, col = "blue", pch = 19, cex = 1.2)
# select the best run for K = 4 clusters 
best = which.min(cross.entropy(project, K = 2)) 
my.colors <- c("tomato", "lightblue", "olivedrab", "gold","pink") 
windows()
barchart(project, K = 2, run = best, border = NA,sort.by.Q=FALSE, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(bp$order), labels = popus, las=2, cex.axis = .75)
q_scores<-Q(project, K = 2, run = best)
q_scores<-as.data.frame(q_scores)
q_scores$Individual<-gsub@ind.names
colnames(q_scores)<-c("oc","swo","Individual")
write.csv(q_scores,"q_scores_2spec_snmf.csv")

p = snmf.pvalues(project, entropy = TRUE, ploidy = 2, K = 4) 
pvalues = p$pvalues 
par(mfrow = c(2,1)) 
hist(pvalues, col = "orange") 
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)

#impute missing values
#impute(project, "gl_geno.lfmm", method = 'mode', K = 4, run = best) 
colnames(ex)
env<-ex[2:17]
coor<-ex[18:19] %>% as.matrix()
length(coor)
getwd()
write.env(env, "sub_env.env")
project=NULL
project = lfmm("sub_geno.geno", "sub_env.env", K = 3, repetitions = 5, project = "new")
load.lfmmProject("sub_geno_sub_env.lfmmProject")
p = lfmm.pvalues(project, K = 3, d=5) 
pvalues = p$pvalues
par(mfrow = c(2,1)) 
hist(pvalues, col = "lightblue") 
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)

project= load.lfmmProject("sub_geno_sub_env.lfmmProject")
lfmm<-read.lfmm("sub_geno.lfmm")
tess3.obj <- tess3(X = lfmm, coord = coor, K = 1:14, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 4, max.iteration = 1000000, rep=10) 

plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix <- qmatrix(tess3.obj, K = 2)
# STRUCTURE-like barplot for the Q-matrix 
pal<-CreatePalette(color.vector = c("tomato", "chartreuse","blue", "gold",
                                    "violet", "olivedrab", "purple", "brown","orange","pink"), palette.length = 4)
barplot(q.matrix, border = NA, sort.by.Q=FALSE,space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), las = 3, cex.axis = .4) 
install.packages("rworldmap")
par(mfrow=c(1,1))
plot(q.matrix, coor, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = pal)
map<-map(database= "county", add = T, interior = T, col="grey",lwd=0.2)
map<-map(database= "state", add = T, interior = T, col="black",lwd=1)
































edge@loc.names
plot_clinecurve(cline$gc, cline.locus=c(edge@loc.names), locus.column = "locus", cline.col="#E495A5", cline.centre.line="CLTA",cline.centre.col="black")
calc_AIC(data.prep.object = prepdata$data.prep,
         esth.object = est2,
         esth.colname = "h_posterior_mode",
         ggcline.object = cline,
         ggcline.pooled = FALSE,
         test.subject)


edge@strata
oakP0<-edge[edge@pop=='CORE-IN5']
oakP1<-edge[edge@pop=='qlyrata_ref']
oakPutHyb<-edge[edge@strata$V3=='EDGE']
oakPutHyb<-gl.drop.ind(oakPutHyb, ind.list=c("qlyrata"))
library(adegenet)
P0_matrix <- as.matrix(oakP0)
P1_matrix<-as.matrix(oakP1)
hyb_matrix<-as.matrix(oakPutHyb)
names<-indNames(oakPutHyb)
#BGC
## install and load devtools
library(devtools)
## install bgc-hm
## this will take a bit (on the order of an hour, depending on what dependencies you already have installed)  as it requires compiling a substantial amount of C++ code
devtools::install_github("zgompert/bgc-hm")
## load the Bgc-hm package
library(bgchm)
?est_genocl()
## estimate parental allele frequencies, uses analytical solution 
p_out<-est_p(G0=P0_matrix,G1=P1_matrix,model="genotype",ploidy="diploid",HMC=FALSE)
g
## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
h_out<-est_hi(Gx=hyb_matrix,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/bgc")
## plot hybrid index estimates with 90% equal-tail probability intervals
## sorted by hybrid index, just a nice way to visualize that in this example we have
## few hybrids with intermediate hybrid indexes

plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=hyb_matrix,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="diploid",hier=TRUE,n_iters=4000)

## how variable is introgression among loci? Lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## examine a plot of the joint posterior distribution for the SDs
pp_plot(objs=gc_out,param1="sdv",param2="sdc",probs=c(0.5,0.75,0.95),colors="black",addPoints=TRUE,palpha=0.1,pdf=FALSE,pch=19)

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)

## summarize loci with credible deviations from genome-average gradients, here the focus is
## specifically on steep clines indicative of loci introgressing less than the average
locinum<-which(sz_out$gradient[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$gradient[,2] > 1) ## number of loci with credibly steep clines

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=hyb_matrix,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

## plot the results
tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)
## note that some individuals appear to be likely backcrosses (close to the outer lines of the triangles)
## but the individals with intermediate hybrid indexes are clearly not F1s but rather late generation hybrids
write.csv(q_out, "bgchm_q_out_run1.csv")
View(q_out)
locn<-g@loc.names %>% as.data.frame()
?seq
locn$num<-seq(1,length(locn$.),1)
loci<-locn[locn$num %in% locinum,]
snp_fil<-gl.keep.loc(g, loci$.)
snp_fil@chromosome
install.packages("snpStats")
gl.ld.haplotype(snp_fil)
# Extract chromosome and position

library(ggplot2)
snp_data <- data.frame(
  SNP = locNames(snp_fil),
  Chromosome = snp_fil@chromosome,    # or gl@other$loc.metrics$chrom
  Position = snp_fil@position         # or gl@other$loc.metrics$position
)

# Make sure Chromosome is a factor for plotting
snp_data$Chromosome <- as.factor(snp_fil$chromosome)

ggplot(snp_data, aes(x = Position, y = Chromosome)) +
  geom_point(alpha = 0.6, size = 0.8) +
  theme_minimal() +
  labs(title = "SNP Positions Across Chromosomes",
       x = "Genomic Position",
       y = "Chromosome")

plot(p_out$p0[1, ], type = "l", main = "Traceplot of alpha (SNP 1)", ylab = "alpha")
plot(p_out$p1[1, ], type = "l", main = "Traceplot of beta (SNP 1)", ylab = "beta")
hist(p_out$p0[1, ], breaks = 50, main = "Posterior of alpha (SNP 1)")
hist(p_out$p1[1, ], breaks = 50, main = "Posterior of beta (SNP 1)")
library(coda)
effectiveSize(mcmc(p_out$p0[1, ]))  # alpha ESS for SNP 1
effectiveSize(mcmc(p_out$p1[1, ]))  # beta ESS for SNP 1
acf(p_out$p0[1, ], main = "Autocorrelation of alpha (SNP 1)")
acf(p_out$p1[1, ], main = "Autocorrelation of beta (SNP 1)")



#function for assessing convergence (From chatgpt)
assess_bgc_convergence <- function(p_out, ess_threshold = 100, acf_threshold = 0.3) {
  nsnps <- nrow(p_out$p0)
  results <- data.frame(
    SNP = 1:nsnps,
    alpha_mean = NA, alpha_sd = NA, alpha_ess = NA, alpha_acf1 = NA,
    beta_mean = NA, beta_sd = NA, beta_ess = NA, beta_acf1 = NA,
    alpha_flag = FALSE, beta_flag = FALSE
  )
  
  for (i in 1:nsnps) {
    alpha_chain <- mcmc(p_out$p0[i, ])
    beta_chain  <- mcmc(p_out$p1[i, ])
    
    # Alpha stats
    results$alpha_mean[i] <- mean(alpha_chain)
    results$alpha_sd[i]   <- sd(alpha_chain)
    results$alpha_ess[i]  <- effectiveSize(alpha_chain)
    results$alpha_acf1[i] <- acf(alpha_chain, plot = FALSE)$acf[2]
    
    # Beta stats
    results$beta_mean[i] <- mean(beta_chain)
    results$beta_sd[i]   <- sd(beta_chain)
    results$beta_ess[i]  <- effectiveSize(beta_chain)
    results$beta_acf1[i] <- acf(beta_chain, plot = FALSE)$acf[2]
    
    # Flag if below ESS or above ACF threshold
    results$alpha_flag[i] <- results$alpha_ess[i] < ess_threshold || results$alpha_acf1[i] > acf_threshold
    results$beta_flag[i]  <- results$beta_ess[i] < ess_threshold || results$beta_acf1[i] > acf_threshold
  }
  
  return(results)
}

summary_df <- assess_bgc_convergence(p_out)

# View worst-converging SNPs
subset(summary_df, alpha_flag | beta_flag)
length(summary_df$alpha_flag=='TRUE')
length(summary_df$beta_flag=='TRUE')
