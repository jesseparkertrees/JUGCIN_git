# r script to compile all data (health, Aziz DNA/coords, extracted env variables, JP DNA data (NH and fastStruc))

# inputs are: (1) CSV of health data in ./data/ folder (2) CSV of Aziz's DNA and coordinate data in ./data/ folder (3) environmental variables in ../EnvPreds/ folder (directory outside of project directory)

#############################################################################################
library(dplyr)
# READ IN AZIZ"S DATA AND FILTER
pops <- read.csv("./data/ButternutLandscapeGPS_Aziz07022026.csv", header=TRUE)

colnames(pops)
pops$JC_total<-pops$JP_fastStruc_K2_JC
pops<-rename(pops,"y_old"="latitude_dd", "x_old"="Longitude_dd")
colnames(pops)

pops$x_old <- as.numeric(pops$x_old)
pops$y_old <- as.numeric(pops$y_old)
sum(is.na(pops$y_old))
sum(is.na(pops$x_old))
sum(!is.na(pops$y_old))
sum(!is.na(pops$x_old))
dim(pops)

#############################################################################################
# READ IN HEALTH DATA AND FILTER, ALIGN SAMPLE NAMES WITH AZIZ'S DATA
health<-read.csv("./data/JUGCIN_Emma_Collected_dupsDel.csv")
health$Latitude<-as.numeric(health$Latitude)
health$Longitude<-as.numeric(health$Longitude)
dim(health)

lookup <- unique(pops$Description)
lookup <- lookup[!is.na(lookup) & lookup != ""]

heal <- health[health$Aziz_descrip %in% lookup, ]
noheal <- health[!(health$Aziz_descrip %in% lookup), ]

dim(heal)
dim(noheal)
heal[duplicated(heal$Aziz_descrip),]
sum(duplicated(heal$Aziz_descrip))

heal$LabID<-NULL
heal<-rename(heal, "Description"="Aziz_descrip")
dim(heal)
length(unique(heal$Description))

pops_health<-left_join(pops,heal,by="Description")
dim(pops_health)
colnames(pops_health)

## update coordinates with those from health data 
pops_health$y_new<-ifelse(is.na(pops_health$Latitude),as.numeric(pops_health$y_old),as.numeric(pops_health$Latitude)) %>% as.numeric()
pops_health$x_new<-ifelse(is.na(pops_health$Longitude),as.numeric(pops_health$x_old),pops_health$Longitude) %>% as.numeric()
dim(pops_health)

sum(is.na(pops_health$LabID))
sum(is.na(pops_health$y_new))
sum(is.na(pops_health$x_new))

## remove locations from Asia and pacific ocean and those without coordinates
pops_health2 <- pops_health[!is.na(pops_health$x_new) & !is.na(pops_health$y_new) & pops_health$x_new < 0 & pops_health$x_new > -125,]
dim(pops_health2)
sum(is.na(pops_health2$y_new))
sum(is.na(pops_health2$x_new))

## write those lacking coordinates to another file
noCoords <- pops_health[is.na(pops_health$x_new) & is.na(pops_health$y_new),]
dim(noCoords)
#write.csv(noCoords, "butternut_no_coords.csv")

## remove locations from plantations (coordinates not correct) 
summary(as.factor(pops_health2$Natural.Forest.Plantation))
summary(as.factor(pops_health2$Situation))

dim(pops_health2)

pops_health3<- pops_health2[pops_health2$Natural.Forest.Plantation!="Plantation" & pops_health2$Natural.Forest.Plantation!="Plantation or private property",]
dim(pops_health3)
summary(as.factor(pops_health3$Natural.Forest.Plantation))
summary(as.factor(pops_health3$Situation))

##################################################################
# JOIN MICROSATELLITE DATA FROM HOBAN TO THIS DATA
## first run "JUGCIN_processMsatData.R"
## then load and filter msat data that includes newhybrids and STRUCTURE results
list.files("./outputs/msatNEWHYBRIDS/")
msat<-read.csv("./outputs/msatNEWHYBRIDS/nhq_all.csv")
dim(msat)
colnames(msat)
msat$GPS.1<-as.numeric(msat$GPS.1)
msat$GPS.2<-as.numeric(msat$GPS.2)
## remove non-natural trees
msat2<-msat %>% filter(Source== "Wild", !is.na(GPS.1))
dim(msat2)
## remove out of range coordinates
summary(msat2$GPS.1)
summary(msat2$GPS.2)
msat2$GPS.2<-ifelse(msat2$GPS.2>0,0-msat2$GPS.2,msat2$GPS.2)
dim(msat2)
summary(msat2$GPS.1)
summary(msat2$GPS.2)
msat2$LabID<-paste0("HobanMSAT_",msat2$SAMPLE.ID )

## remove microsatellite trees that fall within 25m of a GBS sampled tree (assessed using the "Near" tool in ArcGIS Pro)
library(terra)
msat_vect<-vect(msat2, crs='epsg:4326', c('GPS.2','GPS.1'))
aziz_vect<-vect(pops_health3, crs='epsg:4326', c('x_new','y_new'))
# Pairwise distance matrix
d <- distance(msat_vect, aziz_vect)
# TRUE if a point in msat_vect is within 100 m of any point in aziz_vect
within100 <- apply(d <= 100, 1, any)
sum(within100)
# Subset
msat_within100 <- msat_vect[within100, ]
msat_within100_df<-as.data.frame(msat_within100)
msat_over100 <- msat_vect[!within100, ]
msat_over100_df<-as.data.frame(msat_within100)
# also calculated using Near function in arcGIS, can SKIP
#toRemove<-read.csv("./data/msatIDs_toRemove.csv", header=FALSE)
#length(toRemove$V1)
#sum(toRemove$V1 %in% msat_within100_df$LabID)
#sum(!toRemove$V1 %in% msat2$LabID)
#msat3<-msat2[!msat2$LabID %in% toRemove$V1,]

msat3<-msat2[!msat2$LabID %in% msat_over100_df$LabID,]
dim(msat2)
dim(msat3)

## rename columns to match
colnames(msat3)
colnames(pops_health3)
msat3<-rename(msat3, y_new=GPS.1, x_new=GPS.2, DateSampled=DATE.COLLECTED, location=PLACE.COLLECTED, JC_struc= k1andk2, JA_struc=k3_3, DBH_in=DBH..in., PlantHeight_ft=HEIGHT..ft.)

## convert inches to centimeters
msat3$DBH_in<-as.numeric(msat3$DBH_in)
msat3$DBH_cm<-msat3$DBH_in*2.54

## 
msat3<-select(msat3, c(LabID,JA_struc, JC_struc, PlantHeight_ft, DBH_cm, x_new, y_new, NewHyb_FinalAssignment))
pops_health3<-rename(pops_health3, "JC_struc"="JP_fastStruc_K2_JC","JA_struc"="JP_fastStruc_K2_JA")

pops_health3$DBH_cm<-as.numeric(pops_health3$DBH_cm)
pops_health4<-bind_rows(pops_health3,msat3)
dim(pops_health4)

################################################################################
####### OBTAIN AND EXTRACT ENVIRONMENTAL VARIABLES AND ADD TO TABLE  ###########
################################################################################
library(terra)
library(usdm)
library(ggplot2)
library(grid)
library(geodata)
## create 50km buffer around points
coords<-pops_health4 %>% select("LabID","y_new","x_new")
dat<-vect(coords, crs='epsg:4326', c('x_new','y_new'))
chull<-hull(dat, type="concave_ratio", param=0.1, allowHoles=FALSE)
buffer<-buffer(chull, width=50000)
plot(buffer)
points(dat)
extent<-ext(buffer)
extent[1:2]

## download worldclim data using geodata package
#worldclim_global(var='bio', path="../JUGCIN_git_externalFiles/", version="2.1", res=0.5)
#unzip
rasterpath<-"../JUGCIN_git_externalFiles/climate/wc2.1_30s/"
files <- list.files(path=rasterpath, pattern = "^wc2\\.1.*\\.tif$", full.names = TRUE)
bio_stack <- rast(files)

## get "human footprint" data from geodata package: https://www.nature.com/articles/sdata201667
#footprint(year=2009, path="../JUGCIN_git_externalFiles/")
hum_foot<-rast("../JUGCIN_git_externalFiles/landuse/wildareas-v3-2009-human-footprint_geo.tif")

## get soils data
#soil_world(var=c("phh2o", "nitrogen","ocd"),depth=5, path="../JUGCIN_git_externalFiles/")
files<-list.files("../JUGCIN_git_externalFiles/soil_world", full.names = TRUE)
soil<-rast(files)

## crop and combine bioclim, human footprint, and soils data 
stack<-c(bio_stack, hum_foot)
stack1_buf<-crop(stack,buffer,mask=TRUE)
soil_buf<-crop(soil, buffer, mask=TRUE)
stack<-c(stack1_buf,soil_buf)
par(mfrow=c(1,1))
plot(stack$`wildareas-v3-2009-human-footprint`)
points(dat)
plot(stack$`phh2o_0-5cm`)
points(dat)

## extract values
extract_locs<-terra::extract(stack, dat)
dim(extract_locs)
colnames(extract_locs)
length(pops_health4$LabID)
sum(pops_health4$LabID==dat$LabID)
extract_locs$LabID<-pops_health4$LabID
extract_locs$ID<-NULL
colnames(extract_locs)
extract_locs<-extract_locs[,c(24,1:23)]

# get NLCD for North America from <https://www.cec.org/north-american-environmental-atlas/land-cover-30m-2020/>
nlcd<-rast("../JUGCIN_git_externalFiles/land_cover_2020v2_30m_mappackage/nlcd30m_NA.tif")
dat_proj2<-project(dat,nlcd)
chull_proj<-hull(dat_proj2, type="concave_ratio", param=0.1, allowHoles=FALSE)
buffer_proj<-buffer(chull_proj, width=50000)
nlcd_crop<-crop(nlcd, buffer_proj, mask=TRUE)
plot(nlcd_crop)
points(dat_proj2)
levels(nlcd_crop)[[1]]$Class_EN
cat_table <- levels(nlcd_crop)[[1]]
sum(dat_proj2$LabID==dat$LabID)
extract_nlcd <- terra::extract(nlcd_crop, dat_proj2, raw=TRUE) %>% as.data.frame()
head(extract_nlcd)
colnames(extract_nlcd)[2]<-"NLCD"
extract_nlcd$NLCD<-as.factor(extract_nlcd$NLCD)
#group NLCD classes into fewer groups
summary(extract_nlcd$NLCD)
cat_table
extract_nlcd <- extract_nlcd %>%
  mutate(
    NLCD = case_when(
      NLCD %in% c(17) ~ "developed",
      NLCD %in% c(1,2,3,4,5,6) ~ "forest",
      NLCD %in% c(15) ~ "agriculture",
      NLCD %in% c(14) ~ "wetland",
      TRUE ~ "other"
    ),
    NLCD = factor(NLCD)
  )
table(extract_nlcd$NLCD)

# aggregate (mode) NLCD over 1km and 5km to get landscape level effect
nlcd_1km<-aggregate(nlcd_crop, fact=33, fun='modal')
extract_nlcd_1km <- terra::extract(nlcd_1km, dat_proj2, raw=TRUE) %>% as.data.frame()
colnames(extract_nlcd_1km)[2]<-"NLCD_1km"
extract_nlcd_1km$NLCD_1km<-as.factor(extract_nlcd_1km$NLCD_1km)
summary(extract_nlcd_1km$NLCD_1km)
extract_nlcd_1km <- extract_nlcd_1km %>%
  mutate(
    NLCD_1km = case_when(
      NLCD_1km %in% c(17) ~ "developed",
      NLCD_1km %in% c(1,2,3,4,5,6) ~ "forest",
      NLCD_1km %in% c(15) ~ "agriculture",
      NLCD_1km %in% c(14) ~ "wetland",
      TRUE ~ "other"
    ), NLCD_1km = factor(NLCD_1km))
table(extract_nlcd_1km$NLCD_1km)
head(extract_nlcd_1km)
head(extract_nlcd)

nlcd_5km<-aggregate(nlcd_crop, fact=167, fun='modal')
nlcd_5km
extract_nlcd_5km <- terra::extract(nlcd_5km, dat_proj2, raw=TRUE) %>% as.data.frame()
colnames(extract_nlcd_5km)[2]<-"NLCD_5km"
extract_nlcd_5km$NLCD_5km<-as.factor(extract_nlcd_5km$NLCD_5km)
summary(extract_nlcd_5km$NLCD_5km)
extract_nlcd_5km <- extract_nlcd_5km %>%
  mutate(
    NLCD_5km = case_when(
      NLCD_5km %in% c(17) ~ "developed",
      NLCD_5km %in% c(1,2,3,4,5,6) ~ "forest",
      NLCD_5km %in% c(15) ~ "agriculture",
      NLCD_5km %in% c(14) ~ "wetland",
      TRUE ~ "other"), NLCD_5km = factor(NLCD_5km))
table(extract_nlcd_5km$NLCD_5km)
head(extract_nlcd_5km)
head(extract_nlcd)
head(extract_nlcd_1km)

# make custom dist to forest edge layer using NLCD land classifier
levels(nlcd_crop)
plot(nlcd_crop)
plot(nlcd_crop, xlim=c(1500000,1503000), ylim=c(-500000,-503000))

rcl <- matrix(c(
  1, 6, 1,   # classes 1–6 → forest
  7, 19, 0   # classes 4–5 → non-forest
), ncol=3, byrow=TRUE)
lc <- classify(nlcd_crop,rcl)

forest <- ifel(lc == 1, 1, NA)
nonforest <- ifel(lc != 1, 1, NA)

plot(nonforest)
plot(nonforest, xlim=c(1500000,1503000), ylim=c(-500000,-503000))
difeldist_to_edge <- distance(nonforest)

plot(difeldist_to_edge,xlim=c(1500000,1501000), ylim=c(-500000,-501000))
plot(nonforest,xlim=c(1500000,1501000), ylim=c(-500000,-501000))
plot(nlcd_crop, xlim=c(1500000,1501000), ylim=c(-500000,-501000))
plot(forest, xlim=c(1500000,1501000), ylim=c(-500000,-501000))
plot(lc,xlim=c(1500000,1501000), ylim=c(-500000,-501000))
plot(difeldist_to_edge,xlim=c(1500000,1501000), ylim=c(-500000,-501000))

dat_proj<-project(dat,difeldist_to_edge)
extract_disttoedge<-terra::extract(difeldist_to_edge, dat_proj)%>% as.data.frame()
extract_disttoedge$Class_EN<-as.numeric(extract_disttoedge$Class_EN)
colnames(extract_disttoedge)[2]<-"DistToEdge"
summary(extract_disttoedge$DistToEdge)
extract_disttoedge$ID<-NULL
head(extract_disttoedge)

# get dist to forest edge: downloaded from <https://www.sciencebase.gov/catalog/item/5540e3fce4b0a658d79395fe> and merged in arcGIS pro
foredge<-rast("../JUGCIN_git_externalFiles/foredge_merge/foredge_merge.tif")
dat_proj<-project(dat,foredge)
extract_foredge<-terra::extract(foredge, dat_proj)
colnames(extract_foredge)[2]<-"ForestEdge_30m"

#get tree canopy cover
tcc<-rast("../JUGCIN_git_externalFiles/nlcd_tcc_conus_wgs84_v2023-5_19850101_19851231.tif")
summary(tcc)
dat_proj2<-project(dat,tcc)
extract_tcc<-terra::extract(tcc, dat_proj2)%>% as.data.frame()
extract_tcc$category<-as.numeric(extract_tcc$category)
colnames(extract_tcc)[2]<-"TreeCanopyCover"
summary(extract_tcc$TreeCanopyCover)

# aggregate (average) TCC over 1 and 5km to get landscape level effect
tcc_1km<-aggregate(tcc, fact=33, fun='mean')
summary(tcc_1km)
extract_tcc_1km<-terra::extract(tcc_1km, dat_proj2)%>% as.data.frame()
extract_tcc_1km$category<-as.numeric(extract_tcc_1km$category)
colnames(extract_tcc_1km)[2]<-"TreeCanopyCover_1km"
summary(extract_tcc_1km$TreeCanopyCover_1km)

tcc_5km<-aggregate(tcc, fact=167, fun='mean')
extract_tcc_5km<-terra::extract(tcc_5km, dat_proj2)%>% as.data.frame()
extract_tcc_5km$category<-as.numeric(extract_tcc_5km$category)
colnames(extract_tcc_5km)[2]<-"TreeCanopyCover_5km"
summary(extract_tcc_5km$TreeCanopyCover_5km)

# compile all extracts
extract_nlcd
extract_nlcd_1km
extract_nlcd_5km
extract_disttoedge
extract_foredge
extract_tcc
extract_tcc_1km
extract_tcc_5km

allex<-cbind(extract_locs, extract_foredge$ForestEdge_30m, extract_disttoedge$DistToEdge, extract_tcc$TreeCanopyCover,extract_tcc_1km$TreeCanopyCover_1km,extract_tcc_5km$TreeCanopyCover_5km,extract_nlcd$NLCD, extract_nlcd_1km$NLCD_1km,extract_nlcd_5km$NLCD_5km)
head(allex)
colnames(allex)
colnames(allex)[25]<-"ForestEdge_30m"
colnames(allex)[26]<-"ForestEdge_NorAmer_custom"
colnames(allex)[27]<-"TCC"
colnames(allex)[28]<-"TCC_1km"
colnames(allex)[29]<-"TCC_5km"
colnames(allex)[30]<-"NLCD"
colnames(allex)[31]<-"NLCD_1km"
colnames(allex)[32]<-"NLCD_5km"
colnames(allex)

pops_health5<-left_join(pops_health4, allex, by="LabID")
dim(pops_health5)
colnames(pops_health5)

#############################################################################################
# SPATIALLY CLUSTERING POINTS
library(dbscan)
library(sf)
pts <- st_as_sf(coords, coords = c("x_new", "y_new"), crs = 4326)
# project to meters
pts_proj <- st_transform(pts, 3857)
coords_proj <- st_coordinates(pts_proj)
#points within 10k of each other are clustered
cl <- dbscan(coords_proj, eps = 10000, minPts = 2)
pops_health5$CLUSTER_ID<-cl$cluster

#############################################################################################
# get Q values from STRUCTURE and NewHybrids
genetics_summary<-read.csv("genetics_summary.csv")
colnames(genetics_summary)
colnames(pops_health5)
pops_health6<-left_join(pops_health5, genetics_summary %>% select("LabID","NewHyb_FinalAssignment"), by="LabID")
dim(pops_health6)
colnames(pops_health6)
pops_health6<-rename(pops_health6,"NewHyb_31cats"="NewHyb_FinalAssignment.y","NewHyb_msats"="NewHyb_FinalAssignment.x")
genetics_summary2<-read.csv("genetics_summary_37cat.csv")
colnames(genetics_summary2)
pops_health6<-left_join(pops_health6, genetics_summary2 %>% select('LabID', 'NewHyb_FinalAssignment'), by='LabID')
pops_health6<-rename(pops_health6,"NewHyb_37cats"="NewHyb_FinalAssignment")
pops_health6<-relocate(pops_health6, NewHyb_msats, .before=NewHyb_31cats)
summary(as.factor(pops_health6$NewHyb_msats))
pops_health6$NewHyb_msats[pops_health6$NewHyb_msats=='complexF']<-'complexBCF'
pops_health6$NewHyb_msats <- gsub("^NewHyb_", "", pops_health6$NewHyb_msats)
summary(as.factor(pops_health6$NewHyb_37cats))
pops_health6$NewHyb_37cats[pops_health6$NewHyb_37cats=='complexBC_JA']<-'complexBCJA'
pops_health6$NewHyb_37cats[pops_health6$NewHyb_37cats=='complexBC_JC']<-'complexBCJC'
pops_health6$NewHyb_31cats[pops_health6$NewHyb_31cats=='complexBC_JA']<-'complexBCJA'
pops_health6$NewHyb_31cats[pops_health6$NewHyb_31cats=='complexBC_JC']<-'complexBCJC'
pops_health6$NewHyb_31cats <- gsub("^NewHyb_", "", pops_health6$NewHyb_31cats)
summary(as.factor(pops_health6$NewHyb_37cats))
summary(as.factor(pops_health6$NewHyb_31cats))
summary(as.factor(pops_health6$NewHyb_37cats))

pops_health6$NewHyb_FinalAssignment <- coalesce(pops_health6$NewHyb_msats, pops_health6$NewHyb_37cats)

#############################################################################################
# AMMEND, CLEAN UP AND ADD ADMINISTRATIVE DESIGNATIONS
colnames(pops_health6)
pops_health6$y_old<-NULL
pops_health6$x_old<-NULL
pops_health6$Latitude<-NULL
pops_health6$Longitude<-NULL
colnames(pops_health6)

pts <- vect(pops_health6, geom = c("x_new", "y_new"), crs = "EPSG:4326")
## get admin boundaries

usa <- gadm(country = "USA", level = 2, path = "../JUGCIN_git_externalFiles")
can <- gadm(country = "CAN", level = 2, path = "../JUGCIN_git_externalFiles")
## combine USA and CAN
admin <- rbind(usa, can)
result <- terra::extract(admin, pts)
colnames(result)
pops_health6<-cbind(pops_health6, result %>% select(GID_0, COUNTRY, NAME_1, NAME_2, HASC_2))
pops_health6$Country<-NULL

## investigate the 7 points in Quebec that did not populate, they appear to be in the St. Lawrence River near Vaudreuil-Soulanges, update with that.
pops_health6[pops_health6$NAME_1=='',]
pops_health6[1691:1697, "HASC_2"]<-"CA.QC.VS"
pops_health6[1691:1697, "NAME_1"]<-"Québec"
pops_health6[1691:1697, "NAME_2"]<-"Vaudreuil-Soulanges"
pops_health6[1691:1697, "COUNTRY"]<-"Canada"
pops_health6[1691:1697, "GID_0"]<-"CAN"
## investigate point in NY missing data
pops_health6[pops_health6$HASC_2=='',]
pops_health6[1161, 'NAME_2']<-'Niagara'
pops_health6[1161, 'HASC_2']<-"US.NY.NI"

## fix discrepancy in NEWHYB categories
summary(as.factor(pops_health6$NewHyb_FinalAssignment))
#############################################################################################
# EXPORT
write.csv(pops_health6,"pops_health6.csv", na ="", row.names = FALSE)


