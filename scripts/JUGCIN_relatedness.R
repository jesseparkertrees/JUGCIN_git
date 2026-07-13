#export file for cervus
library(vcfR)
library(dartR)
library(SNPfiltR)
library(adegenet)
library(poppr)
library(pegas)
library(dplyr)
library(tidyr)
setwd("C:/Users/jpark107/BUTTERNUT")
list.files()
vc<-read.vcfR("./data/Butternut_Landscape_nJNn9.clean.sub3382.vcf.gz")
colnames(vc@gt)
names<-colnames(vc@gt)
newnames <- sub("\\R1_filtered$", "", basename(names))
length(names)
length(newnames)
head(newnames)
head(names)
colnames(vc@gt)<-newnames
head(colnames(vc@gt))

################################################################
########## extract data for each larger population #############
################################################################
# get ILM genotypes
health_merge<-read.csv("pops_health5_04272026.csv")
health_merge$CLUSTER_ID<-as.factor(health_merge$CLUSTER_ID)
summary(health_merge$CLUSTER_ID)
#get Isle La Motte trees
ILM_inds_all<-health_merge[health_merge$CLUSTER_ID=="3",]
head(ILM_inds_all)
ILM_inds_all$X<-NULL
plot(x=ILM_inds_all$x_new,y=ILM_inds_all$y_new)
ILM_inds_all$PlantHeight_ft<-as.numeric(ILM_inds_all$PlantHeight_ft)
ILM_inds<-ILM_inds_all[!is.na(ILM_inds_all$PlantHeight_ft),]
#ILM_inds<-ILM_inds[ILM_inds$NH_final_assignment=="NewHyb_Butternut",]
dim(ILM_inds)

#get CP trees
CP_inds_all<-health_merge[health_merge$CLUSTER_ID=="7",]
head(CP_inds_all)
CP_inds_all$X<-NULL
plot(x=CP_inds_all$x_new,y=CP_inds_all$y_new)
CP_inds_all<-CP_inds_all[CP_inds_all$LabID!="076795_727_S27_" & CP_inds_all$LabID!="101406_2308_S121_",]
plot(x=CP_inds_all$x_new,y=CP_inds_all$y_new)
CP_inds_all$PlantHeight_ft<-as.numeric(CP_inds_all$PlantHeight_ft)
CP_inds<-CP_inds_all[!is.na(CP_inds_all$PlantHeight_ft),]
#ILM_inds<-ILM_inds[ILM_inds$NH_final_assignment=="NewHyb_Butternut",]
dim(CP_inds)
dim(CP_inds_all)

#get CP trees
WCP_inds_all<-health_merge[health_merge$CLUSTER_ID=="2",]
head(WCP_inds_all)
WCP_inds_all$X<-NULL
plot(x=WCP_inds_all$x_new,y=WCP_inds_all$y_new)
WCP_inds_all$PlantHeight_ft<-as.numeric(WCP_inds_all$PlantHeight_ft)
WCP_inds<-WCP_inds_all[!is.na(WCP_inds_all$PlantHeight_ft),]
dim(WCP_inds)
dim(WCP_inds_all)


################################################################
############  create cervus files for ILM ######################
################################################################
### filter vcf file
vc_filt <- vc[, c("FORMAT", ILM_inds$LabID)]
vc_filt
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt2 <- vc_filt[loci_to_keep,]

# create CERVUS file from vcf
gt <- extract.gt(vc_filt2, element = "GT", as.numeric = FALSE)
# Split genotypes into two alleles
gt_long <- as.data.frame(gt) |>
  mutate(locus = rownames(gt)) |>
  pivot_longer(-locus, names_to = "ID", values_to = "gt") |>
  separate(gt, into = c("A", "B"), sep = "[/|]") |>
  mutate(
    A = ifelse(is.na(A), "0", A),
    B = ifelse(is.na(B), "0", B)
  )
gt_long$locus <- factor(gt_long$locus, levels = rev(unique(gt_long$locus)))

allele_map <- c("0" = "1", "1" = "3")

gt_long <- gt_long |>
  mutate(
    A = allele_map[A],
    B = allele_map[B]
  )

cervus_df <- gt_long |>
  pivot_wider(
    id_cols = ID,
    names_from = locus,
    values_from = c(A, B),
    names_glue = "{locus}{.value}",
    names_vary = "slowest"
  )
View(cervus_df)
dir.create("./cervus/ILM_withhybs_mod052626/")
write.csv(cervus_df,"./cervus/ILM_withhybs_mod052626/cervus_df.csv", row.names = FALSE)

#dir.create("./cervus/ILM_expandedOffspring/")
#write.csv(cervus_df,"./cervus/ILM_expandedOffspring/cervus_df.csv", row.names = FALSE)


################################################################
# generate list of candidate offspring and parents
lifehistory<-select(ILM_inds, c("LabID","PlantHeight_ft","DBH_cm", "Seedling_YN","ProducingSeed","Description","location"))
dim(lifehistory)
colnames(lifehistory)
summary(lifehistory$PlantHeight_ft)
lifehistory$PlantHeight_ft<-as.numeric(lifehistory$PlantHeight_ft)
offspring<-lifehistory[lifehistory$PlantHeight_ft<7, c(1,2)]
#offspring<-lifehistory[lifehistory$PlantHeight_ft<35, c(1,2)]
dim(offspring)
parents<-lifehistory[lifehistory$PlantHeight_ft>15, c(1,2)]
#parents<-lifehistory[lifehistory$PlantHeight_ft>10, c(1,2)]
dim(parents)

offspring<-cbind(offspring, t(parents$LabID))
colnames(offspring)[3:ncol(offspring)]<-parents$LabID
rownames(offspring)<-NULL
head(offspring, 1)
dim(offspring)
dim(parents)
offspring$PlantHeight_ft<-NULL

#write.csv(offspring, "./cervus/ILM_withhybs/offspring.csv", row.names=FALSE)
write.csv(offspring, "./cervus/ILM_withhybs_mod052626/offspring.csv", row.names=FALSE)

plot(x=lifehistory$PlantHeight_ft, y=lifehistory$DBH_cm)

# create files for colony
gl<-vcfR2genlight(vc_filt2)
gl<-gl.compliance.check(gl)

sum(gl@ind.names==datfile$LabID)==length(lifehistory$LabID)
lifehistory$offspring<-ifelse(as.numeric(lifehistory$PlantHeight_ft)<35,"yes","no")
lifehistory$mother<-ifelse(as.numeric(lifehistory$PlantHeight_ft)>10,"yes","no")
lifehistory$father<-ifelse(as.numeric(lifehistory$PlantHeight_ft)>10,"yes","no")


sum(gl@ind.names==datfile$LabID)==length(datfile$LabID)
sum(gl@ind.names==lifehistory$LabID)==length(datfile$LabID)

gl@other$ind.metrics<-cbind(gl@other$ind.metrics, lifehistory[,8:10])
library(dartR.captive)
gl2colony(gl,outfile="colony_input.dat", outpath=getwd(), di.mono.ecious = 1, windows.gui = 1)

################################################################
#########  create cervus files for Charlotte Park ##############
################################################################
### create .gen file
vc_filt <- vc[, c("FORMAT", CP_inds$LabID)]
vc_filt
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt2 <- vc_filt[loci_to_keep,]
gt <- extract.gt(vc_filt2, element = "GT", as.numeric = FALSE)
# Split genotypes into two alleles
gt_long <- as.data.frame(gt) |>
  mutate(locus = rownames(gt)) |>
  pivot_longer(-locus, names_to = "ID", values_to = "gt") |>
  separate(gt, into = c("A", "B"), sep = "[/|]") |>
  mutate(
    A = ifelse(is.na(A), "0", A),
    B = ifelse(is.na(B), "0", B)
  )
gt_long$locus <- factor(gt_long$locus, levels = rev(unique(gt_long$locus)))
allele_map <- c("0" = "1", "1" = "3")
gt_long <- gt_long |>
  mutate(
    A = allele_map[A],
    B = allele_map[B]
  )
cervus_df <- gt_long |>
  pivot_wider(
    id_cols = ID,
    names_from = locus,
    values_from = c(A, B),
    names_glue = "{locus}{.value}",
    names_vary = "slowest"
  )
View(cervus_df)
#dir.create("./cervus/CP_withhybs/")
#write.csv(cervus_df,"./cervus/CP_withhybs/cervus_df.csv", row.names = FALSE)
setwd("../")
getwd()
dir.create("./cervus/CP_expandedOffspring/")
write.csv(cervus_df,"./cervus/CP_expandedOffspring/cervus_df.csv", row.names = FALSE)

################################################################
# generate list of candidate offspring and parents
lifehistory<-select(CP_inds, c("LabID","PlantHeight_ft","DBH_cm", "Seedling_YN","ProducingSeed","Description","location"))
dim(lifehistory)
lifehistory$PlantHeight_ft<-as.numeric(lifehistory$PlantHeight_ft)
#sum(lifehistory$PlantHeight_ft<=15, na.rm=TRUE)
#48
#sum(lifehistory$PlantHeight_ft>15, na.rm=TRUE)
#48
#colnames(lifehistory)
#offspring<-lifehistory[lifehistory$PlantHeight_ft<=15, 1] %>% as.data.frame()
#parents<-lifehistory[lifehistory$PlantHeight_ft>15, 1] %>% as.data.frame()
#offspring<-cbind(offspring, t(parents))
#colnames(offspring)[1]<-"ID"
#colnames(offspring)[2:ncol(offspring)]<-parents$.

offspring<-lifehistory[lifehistory$PlantHeight_ft<35, c(1,2)]
dim(offspring)
parents<-lifehistory[lifehistory$PlantHeight_ft>10, c(1,2)]
dim(parents)
offspring<-cbind(offspring, t(parents$LabID))
colnames(offspring)[3:ncol(offspring)]<-parents$LabID
rownames(offspring)<-NULL
head(offspring, 1)
dim(offspring)
dim(parents)
offspring$PlantHeight_ft<-NULL

dim(offspring)
dim(parents)
write.csv(offspring, "./cervus/CP_expandedOffspring/offspring.csv", row.names=FALSE)




################################################################
########## create cervus files for WCP #########################
################################################################
### create .gen file
vc_filt <- vc[, c("FORMAT", WCP_inds$LabID)]
vc_filt
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt2 <- vc_filt[loci_to_keep,]
gt <- extract.gt(vc_filt2, element = "GT", as.numeric = FALSE)
# Split genotypes into two alleles
gt_long <- as.data.frame(gt) |>
  mutate(locus = rownames(gt)) |>
  pivot_longer(-locus, names_to = "ID", values_to = "gt") |>
  separate(gt, into = c("A", "B"), sep = "[/|]") |>
  mutate(
    A = ifelse(is.na(A), "0", A),
    B = ifelse(is.na(B), "0", B)
  )
gt_long$locus <- factor(gt_long$locus, levels = rev(unique(gt_long$locus)))
allele_map <- c("0" = "1", "1" = "3")
gt_long <- gt_long |>
  mutate(
    A = allele_map[A],
    B = allele_map[B]
  )
cervus_df <- gt_long |>
  pivot_wider(
    id_cols = ID,
    names_from = locus,
    values_from = c(A, B),
    names_glue = "{locus}{.value}",
    names_vary = "slowest"
  )
View(cervus_df)
dir.create("./cervus/WCP_expandedOffspring/")
write.csv(cervus_df,"./cervus/WCP_expandedOffspring/cervus_df.csv", row.names = FALSE)



################################################################
# generate list of candidate offspring and parents
lifehistory<-select(WCP_inds, c("LabID","PlantHeight_ft","DBH_cm", "Seedling_YN","ProducingSeed","Description","location"))
dim(lifehistory)
lifehistory$PlantHeight_ft<-as.numeric(lifehistory$PlantHeight_ft)
#hist(lifehistory$PlantHeight_ft, breaks=20)
#sum(lifehistory$PlantHeight_ft<=15, na.rm=TRUE)
#51
#sum(lifehistory$PlantHeight_ft>15, na.rm=TRUE)
#7
#colnames(lifehistory)
#offspring<-lifehistory[lifehistory$PlantHeight_ft<=15, 1] %>% as.data.frame()
#parents<-lifehistory[lifehistory$PlantHeight_ft>15, 1] %>% as.data.frame()
#offspring<-cbind(offspring, t(parents))
#colnames(offspring)[1]<-"ID"
#colnames(offspring)[2:ncol(offspring)]<-parents$.
#dim(offspring)
#dim(parents)
#write.csv(offspring, "./cervus/WCP_withhybs/offspring.csv", row.names=FALSE)


offspring<-lifehistory[lifehistory$PlantHeight_ft<35, c(1,2)]
dim(offspring)
parents<-lifehistory[lifehistory$PlantHeight_ft>10, c(1,2)]
dim(parents)
offspring<-cbind(offspring, t(parents$LabID))
colnames(offspring)[3:ncol(offspring)]<-parents$LabID
rownames(offspring)<-NULL
head(offspring, 1)
dim(offspring)
dim(parents)
offspring$PlantHeight_ft<-NULL

dim(offspring)
dim(parents)
write.csv(offspring, "./cervus/WCP_expandedOffspring/offspring.csv", row.names=FALSE)







################################################################
############## analyze cervus files ############################
################################################################
# first, open csv files and change column names so no spaces or duplicates, then:
library(dplyr)
library(igraph)
library(ggplot2)
# read files

library(readr)
cols1<-c("OffspringID","LociTyped1","FirstCandidateID","LociTyped2","PairLociCompared1","PairLociMismatching1","PairLODscore1","PairTopLOD1","PairConfidence1","SecondCandidateID","LociTyped3","PairLociCompared2","PairLociMismatching2","PairLODscore2","PairTopLOD2","PairConfidence2","TrioLociCompared","TrioLociMismatching","TrioLOD score","TrioTopLOD","TrioConfidence")
cols2<-c("OffspringID","LociTyped1","FatherID","LociTyped2","PairLociCompared1","PairLociMismatching1","PairLODscore1","CandidateMotherID","LociTyped3","PairLociCompared2","PairLociMismatching2","PairLODscore2","PairTopLOD","PairConfidence","TrioLociCompared","TrioLociMismatching","TrioLOD score","TrioTopLOD","TrioConfidence")

#setwd("C:/Users/jpark107/BUTTERNUT/cervus/ILM_withhybs/")
#setwd("C:/Users/jpark107/BUTTERNUT/cervus/ILM_expandedOffspring/")
setwd("C:/Users/jpark107/BUTTERNUT/cervus/ILM_withhybs_mod052626/")
datfile<-ILM_inds
datfile2<-ILM_inds_all

#setwd("C:/Users/jpark107/BUTTERNUT/cervus/CP_withhybs/")
setwd("C:/Users/jpark107/BUTTERNUT/cervus/CP_expandedOffspring/")
datfile<- CP_inds
datfile2<-CP_inds_all

#setwd("C:/Users/jpark107/BUTTERNUT/cervus/WCP_withhybs/")
setwd("C:/Users/jpark107/BUTTERNUT/cervus/WCP_expandedOffspring/")
datfile<- WCP_inds
datfile2<-WCP_inds_all


################################################################
# look at maternity assignment file
# significant maternity assignments
parpair <- read_csv("./parpairRun.csv", col_names=cols1) %>% slice(-1)
mattest <- read_csv("./matRun.csv", col_names=cols2) %>% slice(-1)

library(ggrepel)
mattest_sig <- mattest %>%
  filter(PairConfidence == "*")
length(mattest_sig$OffspringID)==length(unique(mattest_sig$OffspringID))
colnames(mattest_sig)
# add offspring coordinates
mattest_sig<-rename(mattest_sig,"LabID"="OffspringID")
moms<-left_join(mattest_sig %>% dplyr::select("CandidateMotherID") %>% rename("LabID"="CandidateMotherID"), datfile %>% dplyr::select("LabID","DBH_cm","PlantHeight_ft"), by="LabID") %>% rename("CandidateMotherID"="LabID","CandidateMother_DBH"="DBH_cm","CandidateMother_Height"="PlantHeight_ft")


mattest_sig$CandidateMotherID==moms$CandidateMotherID
mattest_sig<-cbind(mattest_sig, moms$CandidateMother_DBH, moms$CandidateMother_Height)
mattest_sig<-mattest_sig[,c(1:8,20:21,9:19)]
mattest_sig<-rename(mattest_sig,"CandidateMother_DBH" = "moms$CandidateMother_DBH", "CandidateMother_Height" = "moms$CandidateMother_Height")
colnames(mattest_sig)

mattest_sig1 <- mattest_sig %>%
  left_join(
    dplyr::select(datfile, LabID, PlantHeight_ft, DBH_cm, x_new, y_new),
    by = "LabID"
  )
mattest_sig1$DBH_cm_diff<-as.numeric(mattest_sig1$CandidateMother_DBH)-as.numeric(mattest_sig1$DBH_cm)
mattest_sig1$Height_diff<-as.numeric(mattest_sig1$CandidateMother_Height)-as.numeric(mattest_sig1$PlantHeight_ft)

#write.csv(select(mattest_sig1, LabID, CandidateMotherID), "OffspringParentPair.csv", row.names = FALSE)
#write.csv(mattest_sig1, "OffspringParentPair_full.csv", row.names = FALSE)
getwd()
dim(mattest_sig1)
# remove relationships 
mattest_sig2 <- mattest_sig1 %>%
  filter(
    (DBH_cm_diff >= 0 | is.na(DBH_cm_diff)) &
      (Height_diff >= 0 | is.na(Height_diff)))
dim(mattest_sig2)
dim(mattest_sig2)
length(unique(mattest_sig2$LabID))
length(unique(mattest_sig2$CandidateMotherID))

windows()
par(mar=c(9.5,2,1,1))
mother_counts <- sort(summary(as.factor(mattest_sig2$CandidateMotherID)),
                      decreasing = TRUE)
barplot(mother_counts, las = 2)

barplot(summary(as.factor(mattest_sig2$CandidateMotherID)), las=2)
length(unique(mattest$OffspringID))
hist(mattest_sig2$CandidateMother_Height, main= "Histogram of Candidates Parents with Offspring Counts",xlab="Height (ft)")


dim(mattest_sig1)

# EDGE TABLE
edges_cervus <- data.frame(
  from = as.character(mattest_sig2$CandidateMotherID),
  to   = as.character(mattest_sig2$LabID),
  stringsAsFactors = FALSE)
# VERTEX TABLE (ALL IDs NEED COORDS)
all_ids <- data.frame(
  LabID = unique(c(edges_cervus$from, edges_cervus$to)),
  stringsAsFactors = FALSE)

vertices_df <- all_ids %>%
  left_join(
    dplyr::select(datfile, LabID, x_new, y_new),
    by = "LabID"
  ) %>%
  rename(
    name = LabID,
    x = x_new,
    y = y_new
  ) %>%
  distinct(name, .keep_all = TRUE)

# remove vertices with missing coords
vertices_df <- vertices_df %>%
  filter(!is.na(x), !is.na(y))

# keep only edges where both endpoints exist
edges_cervus <- edges_cervus %>%
  filter(from %in% vertices_df$name,
         to   %in% vertices_df$name)
# GRAPH
g <- graph_from_data_frame(
  edges_cervus,
  vertices = vertices_df,
  directed = TRUE)
# EDGE COORDINATES (USE match, NOT row indexing)
edge_plot <- data.frame(
  x    = vertices_df$x[match(edges_cervus$from, vertices_df$name)],
  y    = vertices_df$y[match(edges_cervus$from, vertices_df$name)],
  xend = vertices_df$x[match(edges_cervus$to, vertices_df$name)],
  yend = vertices_df$y[match(edges_cervus$to, vertices_df$name)],
  col = vertices_df$name[match(edges_cervus$from, vertices_df$name)]
)
sum(summary(as.factor(edge_plot$col))>1)

node_plot <- vertices_df %>%
  dplyr::select(x, y, id = name)
allpts<-datfile %>% dplyr::select("x_new","y_new")
# PLOT
p <- ggplot() +
  geom_segment(
    data = edge_plot,
    aes(x = x, y = y, xend = xend, yend = yend, colour= col),
    linewidth = 0.4,
    alpha = 0.7,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed")) +
  geom_point(
    data = allpts,
    aes(x = x_new, y = y_new),
    size = 0.8, 
    shape=3) +
  geom_text_repel(
    data = node_plot,
    aes(x = x, y = y, label = id),
    size = 0.5,
    box.padding = 0.05,
    point.padding = 0.25,
    segment.color = "grey50",
    max.overlaps = Inf) +
  coord_fixed() +
  theme_void() +
  theme(legend.position="none")
p
#ggsave("./maternityPlot.png", plot=p, width = 8, height=6, bg='white', dpi=700)

library(geosphere)
edge_plot$dist_m <- distHaversine(
  p1 = edge_plot[, c("x", "y")],
  p2 = edge_plot[, c("xend", "yend")]
)
edge_plot<-rename(edge_plot,"Parent_ID"="col")
bplot<-ggplot(edge_plot, aes(x = Parent_ID, y = dist_m)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_bw() +
  geom_hline(
    yintercept = mean(edge_plot$dist_m, na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.8) +
  labs(y="Distance (m) between parent and offspring", x="Most Likely Parent") +
  theme(axis.text.x=element_blank())
#ggsave("./ParentOffspringDistPlot.png", plot=bplot, width = 8, height=6, bg='white')
histplot<-ggplot(edge_plot, aes(x = dist_m)) +
  geom_histogram(bins = 20, linewidth = 0.2) +
  geom_vline(
    xintercept = mean(edge_plot$dist_m, na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  theme_bw() +
  labs(
    x = "Distance (m) between parent and offspring",
    y = "Count")
histplot
#ggsave("./ParentOffspringDistHist.png", plot=histplot, width = 6, height=6, bg='white')

library(tidyr)
edge_plot_sum <- edge_plot %>% 
  select(Parent_ID, dist_m) %>%
  group_by(Parent_ID) %>%
  mutate(dist_n = row_number()) %>%   # 1st, 2nd, 3rd distance within each col
  ungroup() %>%
  pivot_wider(
    names_from = dist_n,
    values_from = dist_m,
    names_prefix = "dist_"
  )
colnames(edge_plot_sum)
edge_plot_sum<-rename(edge_plot_sum,"LabID"="Parent_ID")
par_sum<-as.data.frame(summary(as.factor(edge_plot$Parent_ID)))
par_sum$LabID<-rownames(par_sum)
par_sum<-rename(par_sum, "offspring_count"="summary(as.factor(edge_plot$Parent_ID))")
par_sum2<-left_join(par_sum, edge_plot_sum, by="LabID")
par_sum3<-left_join(par_sum2, ILM_inds %>% select(LabID,PlantHeight_ft, DBH_cm, InfectionSeverity, PercentLiveCanopy, AreaInfectedByCanker_Trunk_perc,NH_final_assignment, JC_JP_fastStrucK2), by="LabID")

write.csv(par_sum3,"parent_summary.csv", row.names=FALSE)

par_sum<-left_join(par_sum, datfile, by="LabID")
par_sum<-left_join(par_sum, edge_plot_sum, by="LabID")


par(mar=c(2,2,2,2))
plot(x=par_sum$DBH_cm,y=par_sum$offspring_count)
plot(x=par_sum$PlantHeight_ft,y=par_sum$offspring_count)
plot(x=par_sum$DBH_cm,y=par_sum$offspring_count)
plot(x=par_sum$InfectionSeverity,y=par_sum$offspring_count)
plot(x=par_sum$PercentLiveCanopy,y=par_sum$offspring_count)
par_sum

library(MASS)
plot(par_sum$offspring_count ~ par_sum$InfectionSeverity)
abline(lm(par_sum$offspring_count ~ par_sum$InfectionSeverity))
mod<-glm.nb(offspring_count ~ as.numeric(PlantHeight_ft) + as.numeric(DBH_cm)+ as.numeric(InfectionSeverity) + as.numeric(PercentLiveCanopy),data=par_sum)
summary(mod)

library(quantreg)
mod90 <- rq(par_sum$offspring_count ~ par_sum$InfectionSeverity, tau = 0.9)
summary(mod90)
taus <- c(0.5, 0.75, 0.9)

library(ggplot2)
library(quantreg)

# Remove rows with NA values
dat <- par_sum[complete.cases(par_sum[, c("InfectionSeverity",
                                          "offspring_count")]), ]


# Quantile regression plot
ggplot(dat,
       aes(x = InfectionSeverity,
           y = offspring_count)) +
  geom_point() +
  geom_quantile(quantiles = c(0.5, 0.75, 0.9, 0.99)) +
  theme_classic()
par_sum$AreaInfectedByCanker_Trunk_perc

health <- dat[, c(
  "DBH_cm",
  "PercentLiveCanopy",
  "InfectionSeverity",
  "AreaInfectedByCanker_Trunk_perc"
)]
health$DBH_cm<-as.numeric(health$DBH_cm)
health$PercentLiveCanopy<-health$PercentLiveCanopy %>% as.numeric()
health$InfectionSeverity<-health$InfectionSeverity %>% as.numeric()
health$AreaInfectedByCanker_Trunk_perc<-health$AreaInfectedByCanker_Trunk_perc %>% as.numeric()
pca <- prcomp(health,
              scale. = TRUE)

summary(pca)

dat$HealthPC1 <- pca$x[,1]
mod90 <- rq(offspring_count ~ HealthPC1,
            tau = 0.9,
            data = dat)

summary(mod90)
library(ggplot2)
ggplot(dat, aes(x = HealthPC1,
                y = offspring_count)) +
  geom_point() +
  geom_quantile(quantiles = 0.9) +
  theme_classic()

################################################################
# look at parent pair test
colnames(parpair)
parpair<-rename(parpair,"LabID"="OffspringID")
parpair_sig <- parpair %>%
  filter(TrioConfidence == "*,")
dim(parpair_sig)
length(unique(parpair_sig$LabID))==length(parpair_sig$LabID)

# add offspring coordinates
parpair_sig2 <- parpair_sig %>%
  left_join(datfile %>%
    dplyr::select(LabID, PlantHeight_ft, DBH_cm, x_new, y_new),
    by = "LabID"
  )
colnames(parpair_sig2)
length(unique(parpair_sig2$LabID))
length(unique(parpair$LabID))
length(unique(c(parpair_sig2$FirstCandidateID),parpair_sig2$SecondCandidateID))

moms<-parpair_sig2$FirstCandidateID %>% as.data.frame %>% rename("LabID" = ".")
moms2<-left_join(moms, datfile %>% dplyr::select(LabID,PlantHeight_ft, DBH_cm, NH_final_assignment, JC_JP_fastStrucK2, InfectionSeverity, PercentLiveCanopy, AreaInfectedByCanker_Trunk_perc), by="LabID")
moms2 <- moms2 %>%
  rename_with(
    ~ paste0("FirstCandidate.", .x),
    -LabID) %>% rename("FirstCandidateID"="LabID")
moms2<-moms2[!duplicated(moms2$FirstCandidateID),]

dads<-parpair_sig2$SecondCandidateID %>% as.data.frame %>% rename("LabID" = ".")
dads2<-left_join(dads, datfile %>% dplyr::select(LabID,PlantHeight_ft, DBH_cm, NH_final_assignment, JC_JP_fastStrucK2, InfectionSeverity, PercentLiveCanopy, AreaInfectedByCanker_Trunk_perc), by="LabID")
dads2 <- dads2 %>%
  rename_with(
    ~ paste0("SecondCandidate.", .x),
    -LabID) %>% rename("SecondCandidateID"="LabID")
dads2<-dads2[!duplicated(dads2$SecondCandidateID),]
parpair_sig3<-left_join(parpair_sig2, moms2, by="FirstCandidateID")
parpair_sig4<-left_join(parpair_sig3, dads2, by="SecondCandidateID")
write.csv(parpair_sig4, "parpair_sig4.csv")
getwd()

parents<-bind_rows(moms, dads) %>% distinct(LabID, .keep_all = TRUE)
parents2<-left_join(parents, datfile, by="LabID")
hist(as.numeric(parents2$PercentLiveCanopy), breaks=10)
hist(as.numeric(parents2$InfectionSeverity), breaks=10 )
hist(as.numeric(parents2$AreaInfectedByCanker_Trunk_perc), breaks=10)
par(mar=c(4,2,2,2))
hist(as.numeric(parents2$PlantHeight_ft), breaks=10, main = "Histogram of Height(ft) of Candidate Parents", xlab= "Height (ft)")
hist(as.numeric(parents2$DBH_cm), breaks=10)


# EDGE TABLE
edges_cervus_parpair <- data.frame(
  from = c(as.character(parpair_sig2$FirstCandidateID),as.character(parpair_sig2$SecondCandidateID)),
  to   = as.character(parpair_sig2$LabID),
  stringsAsFactors = FALSE
)
# VERTEX TABLE (ALL IDs NEED COORDS)
all_ids_parpair <- data.frame(
  LabID = unique(c(edges_cervus_parpair$from, edges_cervus_parpair$to)),
  stringsAsFactors = FALSE
)
vertices_df_parpair <- all_ids_parpair %>%
  left_join(
    dplyr::select(datfile, LabID, x_new, y_new),
    by = "LabID"
  ) %>%
  rename(
    name = LabID,
    x = x_new,
    y = y_new
  ) %>%
  distinct(name, .keep_all = TRUE)

# remove vertices with missing coords
vertices_df_parpair <- vertices_df_parpair %>%
  filter(!is.na(x), !is.na(y))

# keep only edges where both endpoints exist
edges_cervus_parpair <- edges_cervus_parpair %>%
  filter(from %in% vertices_df_parpair$name,
         to   %in% vertices_df_parpair$name)
# GRAPH
g <- graph_from_data_frame(
  edges_cervus_parpair,
  vertices = vertices_df_parpair,
  directed = TRUE
)
edge_plot_parpair <- data.frame(
  x    = vertices_df_parpair$x[match(edges_cervus_parpair$from, vertices_df_parpair$name)],
  y    = vertices_df_parpair$y[match(edges_cervus_parpair$from, vertices_df_parpair$name)],
  xend = vertices_df_parpair$x[match(edges_cervus_parpair$to, vertices_df_parpair$name)],
  yend = vertices_df_parpair$y[match(edges_cervus_parpair$to, vertices_df_parpair$name)],
  col = vertices_df_parpair$name[match(edges_cervus_parpair$to, vertices_df_parpair$name)]
)
sum(summary(as.factor(edge_plot_parpair$col))>1)

node_plot_parpair <- vertices_df_parpair %>%
  dplyr::select(x, y, id = name)
# PLOT
p2 <- ggplot() +
  geom_segment(
    data = edge_plot_parpair,
    aes(x = x, y = y, xend = xend, yend = yend, colour= col),
    linewidth = 0.2,
    alpha = 0.7,
    arrow = arrow(
      length = unit(0.05, "inches"))) +
  geom_point(
    data = allpts,
    aes(x = x_new, y = y_new),
    size = 0.8, 
    shape=3) +
  geom_text_repel(
    data = node_plot_parpair,
    aes(x = x, y = y, label = id),
    size = 0.5,
    box.padding = 0.05,
    point.padding = 0.25,
    segment.color = "grey50",
    max.overlaps = Inf) +
  coord_fixed() +
  theme_void() +
  theme(legend.position="none")
ggsave("./ParentPairPlot.png", plot=p2, width = 8, height=6, bg='white', dpi=700)



################################################################
# plot with connections between PARENTS instead
edges_cervus_mates <- data.frame(
  from = c(as.character(parpair_sig2$FirstCandidateID)),
  to   = as.character(parpair_sig2$SecondCandidateID),
  stringsAsFactors = FALSE
)
all_ids_mates <- data.frame(
  LabID = unique(c(edges_cervus_mates$from, edges_cervus_mates$to)),
  stringsAsFactors = FALSE
)

vertices_df_mates <- all_ids_mates %>%
  left_join(
    dplyr::select(datfile, LabID, x_new, y_new),
    by = "LabID"
  ) %>%
  rename(
    name = LabID,
    x = x_new,
    y = y_new
  ) %>%
  distinct(name, .keep_all = TRUE)
g <- graph_from_data_frame(
  edges_cervus_mates,
  vertices = vertices_df_mates,
  directed = TRUE
)
edge_plot_mates <- data.frame(
  x    = vertices_df_mates$x[match(edges_cervus_mates$from, vertices_df_mates$name)],
  y    = vertices_df_mates$y[match(edges_cervus_mates$from, vertices_df_mates$name)],
  par1 = vertices_df_mates$name[match(edges_cervus_mates$from, vertices_df_mates$name)],
  par2 = vertices_df_mates$name[match(edges_cervus_mates$to, vertices_df_mates$name)],
  xend = vertices_df_mates$x[match(edges_cervus_mates$to, vertices_df_mates$name)],
  yend = vertices_df_mates$y[match(edges_cervus_mates$to, vertices_df_mates$name)],
  col = sample(1:length(vertices_df_mates$y[match(edges_cervus_mates$to, vertices_df_mates$name)]))
)
node_plot_mates <- vertices_df_mates %>%
  dplyr::select(x, y, id = name)

p3 <- ggplot() +
  geom_segment(
    data = edge_plot_mates,
    aes(x = x, y = y, xend = xend, yend = yend, colour= col),
    linewidth = 0.2,
    alpha = 0.7) +
  geom_point(
    data = allpts,
    aes(x = x_new, y = y_new),
    size = 0.8, 
    shape=3) +
  geom_text_repel(
    data = node_plot_mates,
    aes(x = x, y = y, label = id),
    size = 1,
    box.padding = 0.05,
    point.padding = 0.25,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  coord_fixed() +
  theme_void()+
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red"))+
  theme(legend.position = "none")
ggsave("./MatesPlot.png", plot=p3, width = 8, height=6, bg='white', dpi=700)

################################################################
# calculate distance between mates
library(geosphere)
edge_plot_mates$dist_m <- distHaversine(
  p1 = edge_plot_mates[, c("x", "y")],
  p2 = edge_plot_mates[, c("xend", "yend")])
colnames(edge_plot_mates)
summary(edge_plot_mates$dist_m)
out <- data.frame(
  Statistic = names(summary(edge_plot_mates$dist_m)),
  Value = as.numeric(summary(edge_plot_mates$dist_m)))
write.csv(out, "./MateDistTable.csv", row.names = FALSE)


################################################################
########## relatedness  ########################################
################################################################
library(maps)
# Extract state (returns "state_name:region")
states <- map.where(database = "state", health_merge$x_new, health_merge$y_new)
states <- sapply(strsplit(states, ":"), "[", 1)
print(states) 
health_merge$state<-states
health_merge$state[is.na(health_merge$state)]<-"new york"

vc_filt <- vc[, c("FORMAT", health_merge$LabID)] #for all points
#vc_filt <- vc[, c("FORMAT", datfile2$LabID)] #for pop level tests
#vc_filt
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt3 <- vc_filt[loci_to_keep,]
gl<-vcfR2genlight(vc_filt3)
gl<-gl.compliance.check(gl)
gl@ind.names==health_merge$LabID
gl@pop<-health_merge$state %>% as.factor()


#sum(gl@ind.names==datfile2$LabID)==length(datfile2$LabID)
gl@strata<-health_merge$PlantHeight_ft %>% as.data.frame
#gl@strata<-datfile2$PlantHeight_ft %>% as.data.frame
colnames(gl@strata)<-"PlantHeight_ft"
rownames(gl@strata)<-indNames(gl)
#gl@pop<-health_merge$Seedling_YN %>% as.factor()
#gl@pop<-datfile2$Seedling_YN %>% as.factor()
#gl@pop<-NULL
png("ALLPTS_GRM.png", width = 20, height = 15, units="in" , res = 700)
grm_gl <- gl.grm(gl)
dev.off()
png("ALLPTS_GRMnetwork_point125.png", width = 20, height = 15, units="in" , res = 700)
grmnet<-gl.grm.network(grm_gl,gl, relatedness_factor = 0.125, method='mds')
dev.off()
summary(as.factor(health_merge$NH_final_assignment))
fil<-health_merge[health_merge$NH_final_assignment=="NewHyb_Butternut",]
vc_filt <- vc[, c("FORMAT", fil$LabID)] #for all points
#vc_filt <- vc[, c("FORMAT", datfile2$LabID)] #for pop level tests
#vc_filt
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt3 <- vc_filt[loci_to_keep,]
gl<-vcfR2genlight(vc_filt)
gl<-gl.compliance.check(gl)
gl@ind.names==fil$LabID
gl@pop<-fil$CLUSTER_ID %>% as.factor()
gl@pop<-fil$state %>% as.factor()
gl@pop==fil$state
gl@strata<-as.data.frame(fil$CLUSTER_ID)
#gl<-gl.keep.pop(gl,"vermont")
#gl@pop<-as.factor(gl@strata$`fil$CLUSTER_ID`)
png("purePTS_GRM.png", width = 20, height = 15, units="in" , res = 700)
grm_gl <- gl.grm(gl)
dev.off()
png("purePTS_GRMnetwork_point25.png", width = 20, height = 15, units="in" , res = 700)
grmnet<-gl.grm.network(grm_gl,gl, relatedness_factor = 0.25, method='mds', node.size = 1, node.label.size)
dev.off()
?gl.grm.network

################################################################
########## compile all pops  ###################################
################################################################
setwd("C:/Users/jpark107/BUTTERNUT/cervus/")

#ilm<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/ILM_withhybs/OffspringParentPair.csv")
#cp<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/CP_withhybs/OffspringParentPair.csv")
#wcp<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/WCP_withhybs/OffspringParentPair.csv")

ilm<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/ILM_expandedOffspring/OffspringParentPair.csv")
cp<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/CP_expandedOffspring/OffspringParentPair.csv")
wcp<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/WCP_expandedOffspring/OffspringParentPair.csv")


comb<-rbind(ilm,cp,wcp)
health_merge_paren<-left_join(health_merge, comb, by="LabID")
health_merge_paren$X<-NULL

ilm_parsum<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/ILM_withhybs/parent_summary.csv")
cp_parsum<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/CP_withhybs/parent_summary.csv")
wcp_parsum<-read.csv("C:/Users/jpark107/BUTTERNUT/cervus/WCP_withhybs/parent_summary.csv")
comb2 <- bind_rows(ilm_parsum, cp_parsum, wcp_parsum)
health_merge_parensum<-left_join(health_merge_paren, comb2, by="LabID")

#write.csv(health_merge_parensum,"pops_health_parentage_05012026.csv", row.names=FALSE)
parsumILM<-left_join(ilm_parsum, health_merge, by="LabID")
#write.csv(parsumILM, "ILM_parentsummary.csv")
plot(parsumILM$x_new, parsumILM$y_new)
cex_vals <- parsumILM$offspring_count / max(parsumILM$offspring_count) * 10  # scale sizes
plot(parsumILM$x_new, parsumILM$y_new, cex = cex_vals, pch = 1, xlab= "x", ylab="y")
legend("bottomright",
       legend = "Circle size scaled based on # of offspring",
       bty = "n")
cex_vals <- as.numeric(parsumILM$TCC) / as.numeric(max(parsumILM$TCC)) * 10  # scale sizes
plot(parsumILM$x_new, parsumILM$y_new, cex = cex_vals, pch = 1, xlab= "x", ylab="y")
legend("bottomright",
       legend = "Circle size scaled based on # of offspring",
       bty = "n")

# Scale point sizes by offspring count
cex_vals <- parsumILM$offspring_count /
  max(parsumILM$offspring_count, na.rm = TRUE) * 5

# Color scale based on TCC
tcc_vals <- parsumILM$TCC
cols <- colorRampPalette(c("red", "blue"))(100)

# Assign colors
col_index <- cut(tcc_vals,
                 breaks = 100,
                 labels = FALSE,
                 include.lowest = TRUE)

point_cols <- cols[col_index]

# Set margins to leave room for color legend
par(mar = c(5, 4, 3, 7))

# Main plot
plot(parsumILM$x_new,
     parsumILM$y_new,
     cex = cex_vals,
     pch = 16,
     col = point_cols,
     xlab = "x",
     ylab = "y")

# Size legend
legend("topleft",
       legend = c("Low offspring", "High offspring"),
       pt.cex = c(.75, 3),
       pch = 16,
       bty = "n",
       title = "# offspring")

# ----- Color legend -----
# Coordinates for legend
usr <- par("usr")

xleft  <- usr[2] - 0.02 * diff(usr[1:2])
xright <- usr[2] + 0.02 * diff(usr[1:2])

ybottom <- usr[3]
ytop    <- usr[4]

# Draw color bar
yseq <- seq(ybottom, ytop, length.out = length(cols) + 1)

for(i in 1:length(cols)) {
  rect(xleft, yseq[i],
       xright, yseq[i + 1],
       col = cols[i],
       border = NA,
       xpd = TRUE)
}

# Add axis labels for TCC scale
tcc_ticks <- pretty(range(tcc_vals, na.rm = TRUE))
tick_pos <- ybottom +
  (tcc_ticks - min(tcc_vals, na.rm = TRUE)) /
  diff(range(tcc_vals, na.rm = TRUE)) *
  (ytop - ybottom)

axis(4,
     at = tick_pos,
     labels = round(tcc_ticks, 2),
     las = 1)

mtext("% Tree\nCanopy\nCover", side = 4, line = 3)




rents <- health_merge_parensum %>%
  filter(!is.na(offspring_count), offspring_count > 0)
offs<- health_merge_parensum %>%
  filter(!is.na(CandidateMotherID))

plot(x= rents$InfectionSeverity,y=rents$offspring_count)
par(mar=c(4,4,2,2))
plot(x= rents$PlantHeight_ft,y=rents$offspring_count, main= "Plant Height vs Offspring Count", xlab= "Plant Height (ft)", ylab="Offspring Count")
plot(x= rents$DBH_cm,y=rents$offspring_count)
plot(x= rents$AreaInfectedByCanker_Trunk_perc,y=rents$offspring_count)
plot(x= rents$AreaInfectedByCanker_RootFlare_prec,y=rents$offspring_count)
plot(x= rents$PercentLiveCanopy,y=rents$offspring_count)



##############################################################

vc_filt <- vc[, c("FORMAT", ILM_inds$LabID)]
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.05
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt2 <- vc_filt[loci_to_keep,]
gt <- extract.gt(vc_filt2, element = "GT")  
ref <- vc_filt2@fix[, "REF"]
alt <- vc_filt2@fix[, "ALT"]
recode_gt <- function(gt_vec, ref, alt) {
  sapply(gt_vec, function(g) {
    if (is.na(g) || g == "./.") return(NA)
    
    alleles <- unlist(strsplit(g, "[/|]"))
    
    # map 0 → REF, 1 → ALT
    alleles <- ifelse(alleles == "0", ref,
                      ifelse(alleles == "1", alt, NA))
    
    paste(alleles, collapse = "/")
  })
}

geno_ac <- matrix(NA, nrow = nrow(gt), ncol = ncol(gt))

for (i in seq_len(nrow(gt))) {
  geno_ac[i, ] <- recode_gt(gt[i, ], ref[i], alt[i])
}
geno_ac_t <- t(geno_ac)

colnames(geno_ac_t) <- paste0(vc_filt2@fix[, "CHROM"], "_", vc_filt2@fix[, "POS"])
rownames(geno_ac_t) <- colnames(gt)
geno_df <- as.data.frame(geno_ac_t, stringsAsFactors = FALSE)
geno_df <- cbind(LabID = rownames(geno_df), geno_df)
rownames(geno_df) <- NULL

ordered_ILM <- ILM_inds[match(geno_df$LabID, ILM_inds$LabID), ]
sum(ordered_ILM$LabID==geno_df$LabID)
par(mar=c(2,2,2,2))
hist(ordered_ILM$PlantHeight_ft)
ordered_ILM$key<-ifelse(ordered_ILM$PlantHeight_ft<10,"Off",ifelse(ordered_ILM$DBH_cm>50,"Pa",ifelse(ordered_ILM$PlantHeight_ft>40,"Pa","All"))) %>% as.factor()
sum(is.na(ordered_ILM$key))
ordered_ILM$key[is.na(ordered_ILM$key)]<-"All"
summary(ordered_ILM$key)
geno_df$key<-ordered_ILM$key
tail(colnames(geno_df))
geno_df <- geno_df[, c(1, ncol(geno_df), 2:(ncol(geno_df) - 1))]
dim(geno_df)
geno_df$LabID.1<-NULL
set.seed(55)
subset<-sample(colnames(geno_df[,3:ncol(geno_df)]),size=300)
geno_df_sub<-cbind(geno_df[,c(1,2)], geno_df[,colnames(geno_df) %in% subset])
dim(geno_df_sub)
#run apparent
apparentOUT <- apparent(geno_df_sub, MaxIdent=0.05, alpha=0.01, nloci=100, self=FALSE, plot=TRUE, Dyad=FALSE)





##########################################################################
#######  try sequoia  ####################################################
##########################################################################
vc_filt <- vc[, c("FORMAT", ILM_inds$LabID)]
maf <- maf(vc_filt) %>% as.data.frame()
str(maf)
threshold <- 0.2
loci_to_keep <- maf$Frequency >= threshold
sum(loci_to_keep)
vc_filt2 <- vc_filt[loci_to_keep,]

gt <- extract.gt(vc_filt2, element = "GT")  

gt_to_numeric <- function(g) {
  if (is.na(g) || g %in% c("./.", ".|.")) return(NA)
  
  alleles <- unlist(strsplit(g, "[/|]"))
  
  # count number of "1" alleles (ALT copies)
  sum(alleles == "1")
}

geno_num <- apply(gt, c(1,2), gt_to_numeric)
dim(geno_num)
geno_sequoia <- t(geno_num)
dim(geno_sequoia)
rownames(geno_sequoia) <- colnames(gt)
geno_sequoia[is.na(geno_sequoia)] <- -9
dim(geno_sequoia)
View(geno_sequoia)
lifehistory<-select(ordered_ILM, c("LabID","PlantHeight_ft"))
lifehistory$Sex<-4
lifehistory$BirthYear <- cut(
  lifehistory$PlantHeight_ft,
  breaks = 4,
  labels = FALSE
)
lifehistory<-rename(lifehistory, "ID"="LabID")
lifehistory$PlantHeight_ft=NULL

library("sequoia")
ageprior<-MakeAgePrior(LifeHistData = data.frame(ID = lifehistory$ID, Sex=lifehistory$Sex, BirthYear=round(lifehistory$BirthYear)))

#Seq_HSg5 <- sequoia(
GenoM = geno_sequoia,
LifeHistData = lifehistory,
Module = "ped",
args.AP = list(Discrete = TRUE, MinAgeParent = 1),
CalcLLR = TRUE,
Plot = FALSE
)

Seq_HSg5 <- sequoia(
  GenoM = geno_sequoia,
  Module = "par",
  CalcLLR = TRUE,
  Plot = FALSE
)

SummarySeq(Seq_HSg5)
View(Seq_HSg5$Pedigree)
View(Seq_HSg5$PedigreePar)

summary(colMeans(geno_sequoia, na.rm = TRUE))



















