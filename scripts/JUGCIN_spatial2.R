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

## edit PlantHeight_ft numbers from Hoban microsatellite data that aren't exact numbers
pops_health6 <- pops_health6 %>%
  mutate(
    PlantHeight_ft = case_when(
      grepl("^\\d+-\\d+$", PlantHeight_ft) ~
        as.character(sapply(strsplit(PlantHeight_ft, "-"),
                            function(x) mean(as.numeric(x)))),
      grepl("^\\d+\\s*\\+$", PlantHeight_ft) ~
        sub("\\s*\\+$", "", PlantHeight_ft),
      TRUE ~ PlantHeight_ft
    ),
    PlantHeight_ft = as.numeric(PlantHeight_ft)
  )

##################################################################
# SUMMARIZE NEWHYBIDS CATEGORIES
## filter out trees that there were not NewHybrids assignments for/ make new final assigment with 31 cats
pops_health6$NewHyb_FinalAssignment<-NULL
pops_health6$NewHyb_31cats[pops_health6$NewHyb_31cats=='complexF']<-'complexBCF'
pops_health6$NewHyb_31cats[pops_health6$NewHyb_31cats=='BCxBC']<-'BCJCxBCJA'
pops_health6$NewHyb_msats[pops_health6$NewHyb_msats == ""] <- NA
pops_health6$NewHyb_msats[pops_health6$NewHyb_msats == "NA"] <- NA

pops_health6$NewHyb_FinalAssignment<-coalesce(pops_health6$NewHyb_msats, pops_health6$NewHyb_31cats)
summary(as.factor(pops_health6$NewHyb_FinalAssignment))


#CHOOSE TO ANALYZE ALL SAMPLES TOGETHER (MSATS and GBS, 2) or just GBS (1) 
#1
#pops_health7<-pops_health6[pops_health6$NewHyb_31cats!="",]
#dim(pops_health7)
#summary(as.factor(pops_health7$NewHyb_31cats))
#pops_health7$NewHyb_FinalAssignment<-NULL
#pops_health7<-rename(pops_health7, 'NewHyb_FinalAssignment'='NewHyb_31cats')

#2
pops_health7<-pops_health6[pops_health6$NewHyb_FinalAssignment!="",]
sum(pops_health6$NewHyb_FinalAssignment=="")
dim(pops_health6)
dim(pops_health7)

## summarize categories
nh_gens<-read.csv("./summaries/NewHybrids_categories_generations_31.csv")
order<-nh_gens$Hybrid_Category
nh_gens$Hybrid_Category <- factor(nh_gens$Hybrid_Category,levels = order)
pops_health7$NewHyb_FinalAssignmentFACTOR<-factor(pops_health7$NewHyb_FinalAssignment, levels=order)
newhyb_table <- pops_health7 %>% count(NewHyb_FinalAssignmentFACTOR, .drop=FALSE)
write.csv(newhyb_table,"./summaries/NewHybridsCategoriesSummaryTableALL.csv", row.names = FALSE)
#write.csv(newhyb_table,"./summaries/summaries_justGBS/NewHybridsCategoriesSummaryTableALL.csv", row.names = FALSE)
p_hist <- pops_health7 %>%
  filter(!NewHyb_FinalAssignment %in% c("JC", "JA")) %>%
  ggplot(aes(x = JC_struc)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(size = 5)
  ) +
  labs(
    x = "Estimated JC ancestry proportion",
    y = "Count",
    caption = "Excludes individuals assigned to the pure JC and JA classes."
  )
p_hist
#ggsave(filename="./summaries/summaries_justGBS/AncestryProportionHistogramNOHYBRIDS.png",plot=p_hist, width=5, height=5, units='in', bg='white' )
ggsave(filename="./summaries/AncestryProportionHistogramNOHYBRIDS.png",plot=p_hist, width=5, height=5, units='in', bg='white' )


newhyb_table<-rename(newhyb_table, Hybrid_Category=NewHyb_FinalAssignmentFACTOR)

nh_gens<-nh_gens %>% left_join(newhyb_table, by='Hybrid_Category')

nh_gens$Expected.Ancestry.Proportion <- ifelse(!is.na(as.numeric(nh_gens$Expected.Ancestry.Proportion)),round(as.numeric(nh_gens$Expected.Ancestry.Proportion), 3),nh_gens$Expected.Ancestry.Proportion)

p<-ggplot(nh_gens, aes(x = Hybrid_Category, y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption=element_text(size=5),
        axis.title.x = element_text(size=10)) +
  labs(x = "Hybrid Category (arranged in order of expected JC ancestry proportion)",y = "Count") +
  geom_text(aes(label = Expected.Ancestry.Proportion), size = 1,color = "hotpink",fontface='bold', vjust=-0.5) +
  labs(caption = "Labels indicate expected JC ancestry proportion\n*these categories are composed of lumped genotypes that failed to converge on a single category")
p
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogramALL.png",plot=p, width=5, height=5, units='in', bg='white' )
ggsave(filename="./summaries/HybridSummariesHistogramALL.png",plot=p, width=5, height=5, units='in', bg='white' )
pmod<-ggplot(nh_gens, aes(x = Hybrid_Category, y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption=element_text(size=5),
        axis.title.x = element_text(size=10)) +
  labs(x = "Hybrid Category (arranged in order of expected JC ancestry proportion)",y = "Count") +
  geom_text(aes(label = paste0(round((n/sum(n))*100,1), "%")), size = 1,color = "hotpink",fontface='bold', vjust=-0.5) +
  labs(caption = "Labels indicate percentage of total in each hybrid category\n*these categories are composed of lumped genotypes that failed to converge on a single category")
pmod
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogramALLperclabel.png",plot=pmod, width=5, height=5, units='in', bg='white' )
ggsave(filename="./summaries/HybridSummariesHistogramALLperclabel.png",plot=pmod, width=5, height=5, units='in', bg='white' )

### how many of total trees sampled are some sort of hybrid?
sum(newhyb_table$n[newhyb_table$Hybrid_Category!='JA' & newhyb_table$Hybrid_Category!='JC'])/sum(newhyb_table$n)
### how many are pure butternut?
sum(newhyb_table$n[newhyb_table$Hybrid_Category=='JC'])/sum(newhyb_table$n)
### how many are pure heartnut?
sum(newhyb_table$n[newhyb_table$Hybrid_Category=='JA'])/sum(newhyb_table$n)

hybrs<-nh_gens[nh_gens$Hybrid_Category!='JA' & nh_gens$Hybrid_Category!='JC',]

p2<-ggplot(hybrs, aes(x = Hybrid_Category, y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption=element_text(size=5)) +
  labs(x = "Hybrid Category (arranged by expected JC ancestry proportion)",y = "Count") +
  geom_text(aes(label = Expected.Ancestry.Proportion), size = 1,color = "hotpink",fontface='bold', vjust=-0.5) +
  labs(caption = "Labels indicate expected JC ancestry proportion\n*these categories are composed of lumped genotypes that failed to converge on a single category")
p2
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogram.png",plot=p2, width=5, height=5, units='in', bg='white' )
ggsave(filename="./summaries/HybridSummariesHistogram.png",plot=p2, width=5, height=5, units='in', bg='white' )

hybrs <- hybrs %>% arrange(Minimum_Number_Hybrid_Generations) %>% mutate(Hybrid_Category = factor(Hybrid_Category,levels = Hybrid_Category))
p3<-ggplot(hybrs, aes(x = Hybrid_Category, y = n)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption=element_text(size=5),
        axis.title.x=element_text(size=10)) +
  labs(x = "Hybrid Category (arranged by minimum number of generations required)",y = "Count") +
  geom_text(aes(label = Expected.Ancestry.Proportion), size = 1,color = "hotpink",fontface='bold', vjust=-0.5) +
  labs(caption = "Labels indicate expected JC ancestry proportion\n*these categories are composed of lumped genotypes that failed to converge on a single category")
p3
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogramArrByGens.png",plot=p3, width=5, height=5, units='in', bg='white' )
ggsave(filename="./summaries/HybridSummariesHistogramArrByGens.png",plot=p3, width=5, height=5, units='in', bg='white' )

library(ggrepel)
p4<-ggplot(hybrs, aes(x=Minimum_Number_Hybrid_Generations,y=n))+
  geom_point() +
  theme_bw() +
  labs(x="Number of hybrid generations required", y="Count")+
  geom_text_repel(aes(label=Hybrid_Category), max.overlaps=Inf, size=1, segment.colour = 'grey',segment.size = 0.4, min.segment.length = 0)
p4
ggsave(filename="./summaries/HybridSummariesHistogramArrByGens2.png",plot=p4, width=5, height=5, units='in', bg='white' )
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogramArrByGens2.png",plot=p4, width=5, height=5, units='in', bg='white' )
p5<-ggplot(hybrs, aes(x=AmongHybridMatings,y=n))+
  geom_point() +
  theme_bw() +
  labs(x="Among Hybrid Matings?", y="Count")+
  geom_text_repel(aes(label=Hybrid_Category), max.overlaps=Inf, size=1, segment.colour = 'grey',segment.size = 0.4, min.segment.length = 0)
p5
#ggsave(filename="./summaries/summaries_justGBS/amongHybridMatings.png",plot=p5, width=5, height=5, units='in', bg='white' )


## summarize categories by "cluster"
colnames(pops_health7)
summary_table_clus <- pops_health7 %>%
  count(CLUSTER_ID, NewHyb_FinalAssignment) %>%
  pivot_wider(
    names_from = NewHyb_FinalAssignment,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(
    pops_health7 %>%
      group_by(CLUSTER_ID) %>%
      summarize(
        mean_JC_struc = mean(JC_struc, na.rm = TRUE), 
        mean_lon = mean(x_new, na.rm = TRUE),
        mean_lat = mean(y_new, na.rm =TRUE),
        .groups="drop"),
    by = "CLUSTER_ID")
summary_table_clus <- summary_table_clus %>% mutate(Hybrids = rowSums(across(-c(CLUSTER_ID, JC, JA, mean_JC_struc, mean_lon, mean_lat))))
summary_table_clus <- summary_table_clus %>% mutate(Total = (rowSums(across(-c(CLUSTER_ID, mean_JC_struc, Hybrids, mean_lon, mean_lat)))))
summary_table_clus <- summary_table_clus %>% mutate(Hybrid_perc = (Hybrids/Total)*100)

summary(summary_table_clus$Hybrid_perc)
hist(summary_table_clus$Hybrid_perc, breaks=20)
#what percent of clusters have more than 0 hybrids?
sum(summary_table_clus$Hybrid_perc>0)/length(summary_table_clus$CLUSTER_ID)
ggplot(data=summary_table_clus, aes(x=Total, y=mean_JC_struc))+
  geom_point() +
  theme_bw() + 
  labs(x="Population Size", y="Average JC Ancestry")
plotclus<-ggplot(summary_table_clus, aes(log10(Total), mean_JC_struc)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x))+
  labs(x="Log10(Number of Individuals in Cluster)", y="Average JC Ancestry")
ggsave("./summaries/clusterSizeVSancesComp.png", plot=plotclus, width=5, height=5, units='in', bg='white')
write.csv(summary_table,"./summaries/NewHybridsCategoriesSummaryTableCLUSTERS.csv", row.names = FALSE)

## summarize categories by county unit
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
#what percent of counties have more than 0 hybrids?
sum(summary_table$Hybrid_perc>0)/length(summary_table$HASC_2)
write.csv(summary_table,"./summaries/summaries_justGBS/NewHybridsCategoriesSummaryTableCounties.csv", row.names = FALSE)
#write.csv(summary_table,"./summaries/NewHybridsCategoriesSummaryTableCounties.csv", row.names = FALSE)
ggplot(data=summary_table, aes(x=Total, y=mean_JC_struc))+
  geom_point() +
  theme_bw() + 
  labs(x="Population Size", y="Average JC Ancestry")
plotcounty<-ggplot(summary_table, aes(log10(Total), mean_JC_struc)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x))+
  labs(x="Log10(Number of Individuals in Cluster)", y="Average JC Ancestry")
# plot on map
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
  geom_sf_text(data = subset(label_pts, !is.na(Total)),aes(label = Total),size = 0.5,color = "hotpink",fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per county") +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5, size=2),
        legend.position=c(0.85,0.35))
plot
ggsave(filename="./summaries/summaries_justGBS/HybridSummariesByCounty_percentHybrid.png",plot=plot, width=5, height=4, units='in', bg='white',dpi=600)
#ggsave(filename="./summaries/HybridSummariesByCounty_percentHybrid.png",plot=plot, width=5, height=4, units='in', bg='white',dpi=600)
plot2<-ggplot(admin1_east) +
  geom_sf(aes(fill = mean_JC_struc), color = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "white", direction=-1) +
  theme_void() + 
  labs(fill="Average\nAncestry\nProportion") +
  geom_sf_text(data=subset(label_pts, !is.na(Total)), aes(label=Total), size=0.5, color='hotpink', fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per county") +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5, size=2),
        legend.position=c(0.85,0.35))
plot2
ggsave(filename="./summaries/summaries_justGBS/HybridSummariesByCounty_ancestryProportion.png",plot=plot2, width=5, height=4, units='in', bg='white', dpi=600)

#ggsave(filename="./summaries/HybridSummariesByCounty_ancestryProportion.png",plot=plot2, width=5, height=4, units='in', bg='white', dpi=600)

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
    by = "NAME_1")
summary_table2 <- summary_table2 %>% mutate(Hybrids = rowSums(across(-c(NAME_1, JC, JA, mean_JC_struc))))
summary_table2 <- summary_table2 %>% mutate(Total = (rowSums(across(-c(NAME_1, mean_JC_struc, Hybrids)))))
summary_table2 <- summary_table2 %>% mutate(Hybrid_perc = (Hybrids/Total)*100)
write.csv(summary_table2,"./summaries/summaries_justGBS/NewHybridsCategoriesSummaryTableStates.csv", row.names = FALSE)

#write.csv(summary_table2,"./summaries/NewHybridsCategoriesSummaryTableStates.csv", row.names = FALSE)
#what percent of states have more than 0 hybrids?
sum(summary_table2$Hybrid_perc>0)/length(summary_table2$NAME_1)
## plot on map
library(rnaturalearth)
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
  geom_sf_text(data = subset(label_pts, !is.na(Total)),aes(label = Hybrids),size = 2,color = "hotpink",fontface='bold') +
  labs(caption = "*Numbers indicate total number of hybrids found per state")+
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
  geom_sf_text(data=subset(admin_east, !is.na(Total)), aes(label=Total), size=2, color='hotpink', fontface='bold') +
  labs(caption = "*Numbers indicate total number of trees sampled per state") +
  theme(plot.caption = element_text(hjust = 1,vjust = 1.5),
        legend.position=c(0.89,0.35))
plot4
ggsave(filename="./summaries/HybridSummariesByState_ancestryProportion.png",plot=plot4, width=5, height=4, units='in', bg='white' )


## HOW DO HYBRID CLASSES CHANGE ACROSS TREE SIZE (using plant height and DBH when available)?
pops_health7
colnames(pops_health7)
pops_health7$NAME_1
summary(as.factor(pops_health7$NewHyb_FinalAssignment))
## create binary hybrid variable
pops_health7$HybridYN<-ifelse(pops_health7$NewHyb_FinalAssignment=="JC" | pops_health7$NewHyb_FinalAssignment=="JA", 0,1 )
summary(as.factor(pops_health7$HybridYN))
summary(as.factor(pops_health7$CLUSTER_ID))

pops_health7group<-pops_health7 %>%
  mutate(
    Hybrid_Grouping = case_when(
      NewHyb_FinalAssignment == "JC" ~ "JC",
      NewHyb_FinalAssignment == "JA" ~ "JA",
      NewHyb_FinalAssignment == "F1" ~ "F1",
      NewHyb_FinalAssignment == "F2" ~ "F2",
      NewHyb_FinalAssignment == "BCJC" ~ "BCJC",
      NewHyb_FinalAssignment == "BCJA" ~ "BCJA",
      NewHyb_FinalAssignment %in% c("BC2JC","BC3JC","BC4JC","BC5JC","BC6JC") ~ "Advanced JC backcross",
      NewHyb_FinalAssignment %in% c("BCJCxBC2JC","BC2JCxBC3JC","BCJCxBC3JC","BCJCxBCJC","complexBCJC") ~ "Advanced complex hybrid JC",
      NewHyb_FinalAssignment %in% c("complexBCF","BCJCxBCJA","BC2JAxBC3JC") ~ "Advanced complex hybrid F" ,
      NewHyb_FinalAssignment %in% c("BC2JA","BC3JA","BC4JA","BC5JA","BC6JA") ~ "Advanced JA backcross",
      NewHyb_FinalAssignment %in% c("complexBCJA","BCJAxBCJA") ~ "Advanced complex hybrid JA"
    )
  )

summary(as.factor(pops_health7group$Hybrid_Grouping))

pops_health7group %>%
  ggplot(aes(x = PlantHeight_ft, color = Hybrid_Grouping)) +
  geom_density(linewidth = 1) +
  geom_density(
    aes(x = PlantHeight_ft, color = "All individuals"),
    linewidth = 1,
    linetype = "dashed",
    color='black'
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Density",
    color = "Hybrid class"
  )

regroup<-pops_health7 %>%
  filter(NewHyb_FinalAssignment!='JA') %>%
  mutate(
    Hybrid_Grouping3 = case_when(
      NewHyb_FinalAssignment == "JC" ~ "JC",
      NewHyb_FinalAssignment == "F1" ~ "F1",
      NewHyb_FinalAssignment %in% c("BC2JC","BC3JC","BC4JC","BC5JC","BC6JC","BCJC") ~ "Mating with JC",
      NewHyb_FinalAssignment %in% c("complexBCF","BCJCxBCJA","BC2JAxBC3JC","BCJCxBC2JC","BC2JCxBC3JC","BCJCxBC3JC","BCJCxBCJC","complexBCJC","complexBCJA","BCJAxBCJA","F2") ~ "Mating among hybrids" ,
      NewHyb_FinalAssignment %in% c("BCJA","BC2JA","BC3JA","BC4JA","BC5JA") ~ "Mating with JA"))
summary(as.factor(regroup$Hybrid_Grouping3))
hybrid_countsregroup <- regroup %>%
  filter(!is.na(PlantHeight_ft)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )
plotdataregroup <- regroup %>%
  left_join(hybrid_countsregroup, by = "Hybrid_Grouping3")
plottimeregroup <- ggplot(plotdataregroup, aes(x = PlantHeight_ft, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = plotdataregroup,
    aes(x = PlantHeight_ft, color = "Average of all trees with height data ()"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Density",
    color = "Hybrid class"
  )

plottimeregroup
















hybridsonly<-pops_health7[pops_health7$HybridYN==1,]
dim(hybridsonly)

hybridsonlygroup<-hybridsonly %>%
  mutate(
    Hybrid_Grouping3 = case_when(
      NewHyb_FinalAssignment == "F1" ~ "F1",
      NewHyb_FinalAssignment %in% c("BC2JC","BC3JC","BC4JC","BC5JC","BC6JC","BCJC") ~ "Mating with JC",
      NewHyb_FinalAssignment %in% c("complexBCF","BCJCxBCJA","BC2JAxBC3JC","BCJCxBC2JC","BC2JCxBC3JC","BCJCxBC3JC","BCJCxBCJC","complexBCJC","complexBCJA","BCJAxBCJA","F2") ~ "Mating among hybrids" ,
      NewHyb_FinalAssignment %in% c("BCJA","BC2JA","BC3JA","BC4JA","BC5JA") ~ "Mating with JA"))
summary(as.factor(hybridsonlygroup$Hybrid_Grouping3))


hybrid_counts <- hybridsonlygroup %>%
  filter(!is.na(PlantHeight_ft)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )

plotdata <- hybridsonlygroup %>%
  left_join(hybrid_counts, by = "Hybrid_Grouping3")

plottime <- ggplot(plotdata, aes(x = PlantHeight_ft, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonlygroup,
    aes(x = PlantHeight_ft, color = "All hybrids (157 hybrids with height data)"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Density",
    color = "Hybrid class"
  )

plottime
ggsave("./summaries/PlantHeightThroughTime_hybridMatings.png", plottime, width=8, height=4, units='in', bg='white' )

hybrid_counts2 <- hybridsonlygroup %>%
  filter(!is.na(DBH_cm)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )

plotdata2 <- hybridsonlygroup %>%
  left_join(hybrid_counts2, by = "Hybrid_Grouping3")
plotdata2<-plotdata2[plotdata2$Hybrid_Grouping3!='Mating with JA',]
plottime2 <- ggplot(plotdata2, aes(x = DBH_cm, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonlygroup,
    aes(x = PlantHeight_ft, color = "All hybrids (157 hybrids with height data)"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Density",
    color = "Hybrid class"
  )

plottime2
ggsave("./summaries/PlantHeightThroughTime_hybridMatings.png", plottime2, width=8, height=4, units='in', bg='white' )



hybridsonlygroup2<-hybridsonly %>%
  mutate(
    Hybrid_Grouping3 = case_when(
      JC_struc < 0.4 ~ "JA majority (<0.4 JC ances.",
      JC_struc < 0.6 ~ "approx. equal (0.4-0.6 JC ances.)",
      JC_struc > 0.6 ~ "JC majority (>0.6 JC ances.)"
      ))
hybrid_counts2 <- hybridsonlygroup2 %>%
  filter(!is.na(PlantHeight_ft)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )

plotdata3 <- hybridsonlygroup2 %>%
  left_join(hybrid_counts2, by = "Hybrid_Grouping3")

plottime3 <- ggplot(plotdata3, aes(x = PlantHeight_ft, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonlygroup2,
    aes(x = PlantHeight_ft, color = "All hybrids (157 hybrids with height data)"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Normalized Density",
    color = "Hybrid class"
  )

plottime3
ggsave("./summaries/PlantHeightThroughTime_structureDensity.png", plottime3, width=8, height=4, units='in', bg='white' )



hybrid_counts3 <- hybridsonlygroup2 %>%
  filter(!is.na(DBH_cm)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )

plotdata4 <- hybridsonlygroup2 %>%
  left_join(hybrid_counts3, by = "Hybrid_Grouping3")

plottime4 <- ggplot(plotdata4, aes(x = DBH_cm, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonlygroup2,
    aes(x = DBH_cm, color = "All hybrids (157 hybrids with height data)"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "DBH (cm)",
    y = "Normalized Density",
    color = "Hybrid class"
  )

plottime4

##################################################################
##################################################################
# SPATIAL MODELING


#check correlation and select relevant variables
library(tidyr)
colnames(pops_health7)
pops_health7_var<-pops_health7 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_1,wc2.1_30s_bio_11,wc2.1_30s_bio_12,wc2.1_30s_bio_3) %>% drop_na() %>% as.data.frame()
cor(pops_health7_var)
library(usdm)
usdm::vif(pops_health7_var)
pops_health7_var<-pops_health7 %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_12, wc2.1_30s_bio_15) %>% drop_na() %>% as.data.frame()
usdm::vif(pops_health7_var)

#rescale data
pops_health7sc <- pops_health7 |>
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
## remove pure JA for modeling
pops_health8<-pops_health7sc[pops_health7sc$NewHyb_FinalAssignment!="JA",]
dim(pops_health7sc)
dim(pops_health8)

# SPATIAL MODELING WITH BINARY HYBRID RESPONSE
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

mod_beta9 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta9)

mod_beta10 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+ph_z+TCC_5km_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta10)

mod_beta11 <- glmmTMB(HybridYN ~  bio6_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta6)

mod_beta12 <- glmmTMB(HybridYN ~  bio15_z*ForestEdge_NorAmer_custom_z*nit_z*TCC_5km_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta12) #no 'significant' interactions

mod_beta13 <- glmmTMB(HybridYN ~  ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta13)

mod_beta14 <- glmmTMB(HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+TCC_5km_z+NLCD+(1|CLUSTER_ID),family = binomial(),data = pops_health8)
summary(mod_beta14)

r2(mod_beta6)
r2(mod_beta7)
r2(mod_beta8)
r2(mod_beta9)
r2(mod_beta10)
r2(mod_beta11)
r2(mod_beta12)
r2(mod_beta13)
r2(mod_beta14)

AIC(mod_beta, mod_beta2, mod_beta3, mod_beta4, mod_beta5, mod_beta6, mod_beta7, mod_beta8, mod_beta9, mod_beta10, mod_beta11, mod_beta12, mod_beta13, mod_beta14)

summary(mod_beta6)
plot(ggpredict(mod_beta6, terms = "TCC_5km_z"))
plot(ggpredict(mod_beta6, terms = c("ForestEdge_NorAmer_custom_z","bio15_z")))
plot(ggpredict(mod_beta6, terms = "bio15_z"))
plot(ggpredict(mod_beta6, terms = "TCC_5km_z"))
plot(ggpredict(mod_beta6, terms = "nit_z"))

pred1 <- ggpredict(mod_beta6, terms = "TCC_5km_z") %>%
  mutate(variable = "TCC")
pred2 <- ggpredict(mod_beta6, terms = "bio15_z") %>%
  mutate(variable = "BIO15")
pred3 <- ggpredict(mod_beta6, terms = "nit_z") %>%
  mutate(variable = "Nitrogen")
pred4 <- ggpredict(mod_beta6, terms = "ForestEdge_NorAmer_custom_z") %>%
  mutate(variable = "Distance to Forest Edge")
preds <- bind_rows(pred1, pred2, pred3, pred4)
preds <- preds |> 
  dplyr::filter(!is.na(variable))
ggplot(preds, aes(x, predicted)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_x") +
  theme_bw()


par(mar=c(8,8,8,8))
boxplot(pops_health8$human_footprint_z~pops_health8$NewHyb_FinalAssignment, las=2)
boxplot(pops_health8$human_footprint_z~pops_health8$HybridYN, las=2)
boxplot(pops_health8$TCC_5km~pops_health8$HybridYN, las=2)
boxplot(pops_health8$TCC_1km~pops_health8$NewHyb_FinalAssignment, las=2)

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




#####################################ARCHIVE###############################################
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

