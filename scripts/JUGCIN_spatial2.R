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
library(ggeffects)
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
pops_health6<-pops_health6 %>% mutate(seq_type = ifelse(substr(LabID, 1,5)=='Hoban',1,0))
colnames(pops_health6)
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
## add variable with size of population and change 1 ind clusters to unique values
max_id <- max(pops_health6$CLUSTER_ID, na.rm = TRUE)
pops_health6 <- pops_health6 %>%
  mutate(
    CLUSTER_ID = if_else(
      CLUSTER_ID == 0,
      max_id + cumsum(CLUSTER_ID == 0),
      CLUSTER_ID) %>% as.factor()) %>%
  add_count(CLUSTER_ID, name = "CLUSTER_SIZE") %>%
  mutate(CLUSTER_SIZE_log = log10(CLUSTER_SIZE))
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

pops_health6 %>% filter(NewHyb_FinalAssignment!='JA' & NewHyb_FinalAssignment!='JC') %>% count(NewHyb_FinalAssignment) %>% mutate(perc = 100*n/sum(n))


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
newhyb_table <- pops_health7 %>% count(NewHyb_FinalAssignmentFACTOR, .drop=FALSE) %>%
  mutate(n_perc = 100*n/sum(n))
#write.csv(newhyb_table,"./summaries/NewHybridsCategoriesSummaryTableALL.csv", row.names = FALSE)
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
#ggsave(filename="./summaries/AncestryProportionHistogramNOHYBRIDS.png",plot=p_hist, width=5, height=5, units='in', bg='white' )


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
#ggsave(filename="./summaries/HybridSummariesHistogramALL.png",plot=p, width=5, height=5, units='in', bg='white' )
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
#ggsave(filename="./summaries/HybridSummariesHistogramALLperclabel.png",plot=pmod, width=5, height=5, units='in', bg='white' )

### how many of total trees sampled are some sort of hybrid?
sum(newhyb_table$n[newhyb_table$Hybrid_Category!='JA' & newhyb_table$Hybrid_Category!='JC'])/sum(newhyb_table$n)
### how many are pure butternut?
sum(newhyb_table$n[newhyb_table$Hybrid_Category=='JC'])/sum(newhyb_table$n)
### how many are pure heartnut?
sum(newhyb_table$n[newhyb_table$Hybrid_Category=='JA'])/sum(newhyb_table$n)

hybrs<-nh_gens[nh_gens$Hybrid_Category!='JA' & nh_gens$Hybrid_Category!='JC',]
hybrs$n_perc<-100*hybrs$n/sum(hybrs$n)
#write.csv(hybrs,"./summaries/NewHybridsCategoriesSummaryTableHYBRIDSonly.csv", row.names = FALSE)
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
#ggsave(filename="./summaries/HybridSummariesHistogram.png",plot=p2, width=5, height=5, units='in', bg='white' )

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
#ggsave(filename="./summaries/HybridSummariesHistogramArrByGens.png",plot=p3, width=5, height=5, units='in', bg='white' )

library(ggrepel)
p4<-ggplot(hybrs, aes(x=Minimum_Number_Hybrid_Generations,y=n))+
  geom_point() +
  theme_bw() +
  labs(x="Number of hybrid generations required", y="Count")+
  geom_text_repel(aes(label=Hybrid_Category), max.overlaps=Inf, size=1, segment.colour = 'grey',segment.size = 0.4, min.segment.length = 0)
p4
#ggsave(filename="./summaries/HybridSummariesHistogramArrByGens2.png",plot=p4, width=5, height=5, units='in', bg='white' )
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesHistogramArrByGens2.png",plot=p4, width=5, height=5, units='in', bg='white' )
p5<-ggplot(hybrs, aes(x=AmongHybridMatings,y=n))+
  geom_point() +
  theme_bw() +
  labs(x="Among Hybrid Matings?", y="Count")+
  geom_text_repel(aes(label=Hybrid_Category), max.overlaps=Inf, size=1, segment.colour = 'grey',segment.size = 0.4, min.segment.length = 0)
p5
#ggsave(filename="./summaries/summaries_justGBS/amongHybridMatings.png",plot=p5, width=5, height=5, units='in', bg='white' )
hybrs<-hybrs %>% mutate(
  perc = 100*n/sum(n)
)

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
#ggsave("./summaries/clusterSizeVSancesComp.png", plot=plotclus, width=5, height=5, units='in', bg='white')
#write.csv(summary_table,"./summaries/NewHybridsCategoriesSummaryTableCLUSTERS.csv", row.names = FALSE)

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
#write.csv(summary_table,"./summaries/summaries_justGBS/NewHybridsCategoriesSummaryTableCounties.csv", row.names = FALSE)
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
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesByCounty_percentHybrid.png",plot=plot, width=5, height=4, units='in', bg='white',dpi=600)
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
#ggsave(filename="./summaries/summaries_justGBS/HybridSummariesByCounty_ancestryProportion.png",plot=plot2, width=5, height=4, units='in', bg='white', dpi=600)

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
#write.csv(summary_table2,"./summaries/summaries_justGBS/NewHybridsCategoriesSummaryTableStates.csv", row.names = FALSE)

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
#ggsave(filename="./summaries/HybridSummariesByState_percentHybrid.png",plot=plot3, width=5, height=4, units='in', bg='white' )
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
#ggsave(filename="./summaries/HybridSummariesByState_ancestryProportion.png",plot=plot4, width=5, height=4, units='in', bg='white' )


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

pops_health7group<-pops_health7group %>%
  mutate(
    Hybrid_Grouping3 = case_when(
      NewHyb_FinalAssignment == "JA" ~ "JA",
      NewHyb_FinalAssignment == "JC" ~ "JC",
      NewHyb_FinalAssignment == "F1" ~ "F1",
      NewHyb_FinalAssignment %in% c("BC2JC","BC3JC","BC4JC","BC5JC","BC6JC","BCJC") ~ "Mating with JC",
      NewHyb_FinalAssignment %in% c("complexBCF","BCJCxBCJA","BC2JAxBC3JC","BCJCxBC2JC","BC2JCxBC3JC","BCJCxBC3JC","BCJCxBCJC","complexBCJC","complexBCJA","BCJAxBCJA","F2") ~ "Mating among hybrids" ,
      NewHyb_FinalAssignment %in% c("BCJA","BC2JA","BC3JA","BC4JA","BC5JA") ~ "Mating with JA"),
    Hybrid_Grouping4 = case_when(
      JC_struc < 0.4 ~ "JA majority (<0.4 JC ances.",
      JC_struc < 0.6 ~ "approx. equal (0.4-0.6 JC ances.)",
      JC_struc > 0.6 ~ "JC majority (>0.6 JC ances.)")
  )
regroup<-filter(pops_health7group, Hybrid_Grouping3!='JA')
summary(as.factor(regroup$Hybrid_Grouping4))
variable<-"Hybrid_Grouping4"
#variable<-"Hybrid_Grouping3"
hybrid_countsregroup <- regroup %>%
  filter(!is.na(PlantHeight_ft)) %>%
  count(.data[[variable]]) %>%
  mutate(
    Hybrid_Grouping_label = paste0(.data[[variable]], " (n=", n, ")")
  )
plotdataregroup <- regroup %>%
  left_join(hybrid_countsregroup, by = variable)
plottimeregroup <- ggplot(plotdataregroup, aes(x = PlantHeight_ft, color = Hybrid_Grouping_label)) +
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

hybridsonly<-pops_health7group[pops_health7group$HybridYN==1,]
dim(hybridsonly)
summary(as.factor(hybridsonly$Hybrid_Grouping3))
hybrid_counts <- hybridsonly %>%
  filter(!is.na(PlantHeight_ft)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )
plotdata <- hybridsonly %>%
  left_join(hybrid_counts, by = "Hybrid_Grouping3")
plottime <- ggplot(plotdata, aes(x = PlantHeight_ft, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonly,
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
#ggsave("./summaries/PlantHeightThroughTime_hybridMatings.png", plottime, width=8, height=4, units='in', bg='white' )
hybrid_counts2 <- hybridsonly %>%
  filter(!is.na(DBH_cm)) %>%
  count(Hybrid_Grouping3) %>%
  mutate(
    Hybrid_Grouping3_label = paste0(Hybrid_Grouping3, " (n=", n, ")")
  )
plotdata2 <- hybridsonly %>%
  left_join(hybrid_counts2, by = "Hybrid_Grouping3")
plotdata2<-plotdata2[plotdata2$Hybrid_Grouping3!='Mating with JA',]
plottime2 <- ggplot(plotdata2, aes(x = DBH_cm, color = Hybrid_Grouping3_label)) +
  geom_density(linewidth = 1) +
  geom_density(
    data = hybridsonly,
    aes(x = PlantHeight_ft, color = "All hybrids (157 hybrids with height data)"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Plant height (ft)",
    y = "Density",
    color = "Hybrid class")
plottime2
#ggsave("./summaries/PlantHeightThroughTime_hybridMatings.png", plottime2, width=8, height=4, units='in', bg='white' )
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
#ggsave("./summaries/PlantHeightThroughTime_structureDensity.png", plottime3, width=8, height=4, units='in', bg='white' )


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
##################################################################
# SPATIAL MODELING
#check correlation and select relevant variables
library(tidyr)

colnames(pops_health7group)
pops_health7_var<-pops_health7group %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_1,wc2.1_30s_bio_11,wc2.1_30s_bio_12,wc2.1_30s_bio_3) %>% drop_na() %>% as.data.frame()
cor(pops_health7_var)
library(usdm)
usdm::vif(pops_health7_var)
pops_health7_var<-pops_health7group %>% select(wildareas.v3.2009.human.footprint, nitrogen_0.5cm, ocd_0.5cm, phh2o_0.5cm, ForestEdge_30m, ForestEdge_NorAmer_custom, TCC, wc2.1_30s_bio_6,wc2.1_30s_bio_12, wc2.1_30s_bio_15) %>% drop_na() %>% as.data.frame()
usdm::vif(pops_health7_var)

# remove the three sites that only have JA...not relevant
pops_health7FILT<-pops_health7group %>% filter(!(NewHyb_FinalAssignment=='JA' & CLUSTER_SIZE==1))
dim(pops_health7FILT)
dim(pops_health7group)

#remove city level coordinate and recalculate cluster size
pops_health7FILT<-pops_health7FILT %>% filter(Situation!='City level')
summary(as.factor(pops_health7FILT$Situation))
## update cluster size
pops_health7FILT2 <- pops_health7FILT %>%
  add_count(CLUSTER_ID, name = "CLUSTER_SIZE_new") %>%
  mutate(CLUSTER_SIZE_new_log = log10(CLUSTER_SIZE_new))
dim(pops_health7FILT2)

# is population size skewed?
hist(pops_health7FILT2$CLUSTER_SIZE) #yes
hist(log10(pops_health7FILT2$CLUSTER_SIZE))
summary(pops_health7FILT2$CLUSTER_SIZE)
hist(scale(log10(pops_health7FILT2$CLUSTER_SIZE)))
# rescale data
colnames(pops_health7FILT2)
pops_health7sc <- pops_health7FILT2 %>%
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
    seq_type = as.factor(seq_type),
    CLUSTER_SIZE_log_z = scale(CLUSTER_SIZE_log),
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

## does population size have an effect on hybrid status?
summary(as.factor(pops_health7sc$Hybrid_Grouping3))
pops_health7sc<-pops_health7sc %>% mutate(
  Hybrid_Grouping3B = case_when(Hybrid_Grouping3 == 'Mating among hybrids' ~ 0,
                                Hybrid_Grouping3 == 'Mating with JC' ~ 1,
                                Hybrid_Grouping3 == 'Mating with JA' ~ -1,
                                Hybrid_Grouping3 == 'JC'~NA))
pop_ave <- pops_health7sc %>%
  group_by(CLUSTER_ID) %>%
  summarise(
    avg = mean(JC_struc, na.rm = TRUE),
    n = n(),
    CLUSTER_SIZE = first(CLUSTER_SIZE),
    CLUSTER_SIZE_log = first(CLUSTER_SIZE_log),
    .groups = "drop")
pop_ave
library(glmmTMB)
mod_beta<-glmmTMB(avg ~ CLUSTER_SIZE_log,family=beta_family(),data = pop_ave)
summary(mod_beta)
pred <- ggpredict(mod_beta, terms = "CLUSTER_SIZE_log")
# plot transformed data
rego<-ggplot(pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  theme_bw() +
  labs(x = "log10(Population size)", y = "JC ancestry proportion") +
  geom_point(
    data = pop_ave,
    aes(x = CLUSTER_SIZE_log, y = avg),
    inherit.aes = FALSE,
    alpha = 0.5)
#ggsave("clusterSizeVSancesComp.png",plot=rego, width=5, height=4, units='in', bg='white')
## when including all individuals, with increasing population size comes increased chance of admixture (lower JC ances. propo.) (because there are many small populations (1-5) that are isolated from JA and pure...larger populations are more likely to have SOME admixture), BUT if only looking at populations that contain at least one hybrid:
hist(pop_ave$avg[pop_ave$CLUSTER_SIZE==1]) #weakly bimodal...amongst singleton pops there are many pure individuals (no opportunity for interspecific mating due to isolation(?), but there are also many relatively highly admixed pops (smaller/ isolated pops more susceptible to introgression due to swamping)
#only include populations that have at least one hybrid

pop_ave2 <- pops_health7sc %>%
  group_by(CLUSTER_ID) %>%
  filter(any(!NewHyb_FinalAssignment %in% c("JC", "JA"))) %>%
  summarise(
    avg = mean(JC_struc, na.rm = TRUE),
    n = n(),
    n_HYB = sum(!NewHyb_FinalAssignment %in% c("JC", "JA"), na.rm=TRUE),
    CLUSTER_SIZE = first(CLUSTER_SIZE),
    CLUSTER_SIZE_log = first(CLUSTER_SIZE_log),
    .groups = "drop"
  )
sum(pop_ave2$n_HYB)
hist(pop_ave2$avg[pop_ave2$CLUSTER_SIZE==1])
mod_betaB<-glmmTMB(avg ~ CLUSTER_SIZE_log,family=beta_family(),data = pop_ave2)
summary(mod_betaB)
pval <- summary(mod_betaB)$coefficients$cond["CLUSTER_SIZE_log", "Pr(>|z|)"]
library(performance)
r2_vals <- r2(mod_betaB)
label <- sprintf("p = %.3g\nFerrari R² = %.3f",pval,r2_vals$R2)
pred <- ggpredict(mod_betaB, terms = "CLUSTER_SIZE_log")
reg1<-ggplot(pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(data = pop_ave2,aes(x = CLUSTER_SIZE_log, y = avg),inherit.aes = FALSE,alpha = 0.5) +
  annotate("text",x = Inf,y = 0.25,label = label,hjust = 1.1,vjust = 1.1,size = 4) +
  theme_bw() +
  labs(x = "log10(Population size) of populations containing > 0 hybrids",y = "Population average JC ancestry proportion")
reg1
ggsave("./summaries/populationSizeVSJCances.png",plot=reg1, width=5, height=4, units='in', bg='white')
##so, amongst populations containing hybrids, population size is strongly correlated with average JC ancestry...greater admixture among hybridized genotypes (F2, F3, etc.) may depress average JC ancestry. Larger populations have more JC with which to backcross and elevate the JC proportion over generations... 
dim(pop_ave2)
dim(pop_ave)
dim(pop_ave2[pop_ave2$CLUSTER_SIZE>1,])

pop_aveJA <- pops_health7sc %>%
  group_by(CLUSTER_ID) %>%
  filter(any(NewHyb_FinalAssignment=="JA")) %>%
  summarise(
    avg = mean(JC_struc, na.rm = TRUE),
    n = n(),
    n_JA = sum(NewHyb_FinalAssignment == 'JA', na.rm=TRUE),
    CLUSTER_SIZE = first(CLUSTER_SIZE),
    CLUSTER_SIZE_log = first(CLUSTER_SIZE_log),
    .groups = "drop"
  )

hyb_no_JA <- setdiff(pop_ave2$CLUSTER_ID, pop_aveJA$CLUSTER_ID)

site_summary <- data.frame(
  Category = c(
    "All clusters",
    "Sites with at least one hybrid",
    "Sites with JA",
    "Sites with at least one hybrid but no JA"),
  N = c(
    length(pop_ave$CLUSTER_ID),
    length(pop_ave2$CLUSTER_ID),
    length(pop_aveJA$CLUSTER_ID),
    length(hyb_no_JA))) %>% 
  mutate(Percent = round(100 * N / length(pop_ave$CLUSTER_ID), 1))
write.csv(site_summary,"./summaries/cluster_summary.csv", row.names = FALSE)

#############################################################
# SPATIAL MODELING WITH BINARY HYBRID RESPONSE
## remove pure JA for modeling
summary(as.factor(pops_health7sc$Natural.Forest.Plantation))
summary(as.factor(pops_health7sc$Situation))
# filter out pure JA and trees without exact coordinates
pops_health8<-pops_health7sc %>% filter(NewHyb_FinalAssignment!='JA' & Situation!='City level')
colnames(pops_health8)
library("glmmTMB")
library(performance)
library("ggeffects")
#make forest the reference variable for NLCD land cover (makes this the dummy variable)
pops_health8$NLCD_5km <- relevel(factor(pops_health8$NLCD_5km),ref = "forest")
models<-list(
  mod_beta1 = HybridYN ~ ForestEdge_30m_z*TCC_5km_z*bio6_z*nit_z+(1|CLUSTER_ID)+seq_type,
  mod_beta2 = HybridYN ~ ForestEdge_NorAmer_custom_z*TCC_5km_z*bio12_z*nit_z+(1|CLUSTER_ID)+seq_type,
  mod_beta3 = HybridYN ~ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type,
  mod_beta5 = HybridYN ~ bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type,
  mod_beta6 = HybridYN~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type,
  mod_beta7 = HybridYN ~ bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|CLUSTER_ID)+seq_type,
  mod_beta8 = HybridYN ~ bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|NAME_1)+seq_type,
  mod_beta9 = HybridYN ~  bio15_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+(1|CLUSTER_ID)+seq_type,
  mod_beta10 =HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+nit_z+ph_z+TCC_5km_z+(1|CLUSTER_ID)+seq_type,
  mod_beta11 =HybridYN ~  bio6_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type,
  mod_beta12 =HybridYN ~  bio15_z*ForestEdge_NorAmer_custom_z*nit_z*TCC_5km_z+(1|CLUSTER_ID)+seq_type,
  mod_beta13 =HybridYN ~  ForestEdge_NorAmer_custom_z+nit_z+TCC_5km_z+(1|CLUSTER_ID)+seq_type,
  mod_beta14 =HybridYN ~  bio15_z+ ForestEdge_NorAmer_custom_z+TCC_5km_z+NLCD+(1|CLUSTER_ID)+seq_type,
  mod_beta15 =HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type,
  mod_beta16 =HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_NorAmer_custom_z+nit_z+TCC_z+NLCD+(1|CLUSTER_ID)+seq_type,
  mod_beta17 =HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_z+NLCD_5km+(1|CLUSTER_ID)+seq_type
)
mods <- lapply(models, glmmTMB, family = binomial(), data = pops_health8)
lapply(mods, summary)
summary(mods$mod_beta12) #no 'significant' interactions
lapply(mods,r2, verbose=FALSE)
lapply(mods,AIC)

modtest1<-glmmTMB(HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID)+seq_type, family = binomial(), data = pops_health8)
modtest2<-glmmTMB(HybridYN ~  bio12_z+bio6_z+bio15_z+human_footprint_z+ ForestEdge_30m_z+nit_z+TCC_5km_z+NLCD_5km+(1|CLUSTER_ID), family = binomial(), data = pops_health8)
AIC(modtest1, modtest2)

summary(mods$mod_beta17)

summary(mods$mod_beta15)
selectedMod<-mods$mod_beta15
pred1 <- ggpredict(selectedMod, terms = "TCC_5km_z") %>%
  mutate(variable = "TCC")
pred2 <- ggpredict(selectedMod, terms = "bio15_z") %>%
  mutate(variable = "BIO15")
pred3 <- ggpredict(selectedMod, terms = "nit_z") %>%
  mutate(variable = "Nitrogen")
pred4 <- ggpredict(selectedMod, terms = "ForestEdge_30m_z") %>%
  mutate(variable = "Distance to Forest Edge")
preds <- bind_rows(pred1, pred2, pred3, pred4)
preds <- preds |> 
  dplyr::filter(!is.na(variable))
pvals<-summary(selectedMod)$coefficients$cond[,'Pr(>|z|)']
pval_df <- data.frame(
  variable = c("TCC","BIO15","Nitrogen","Distance to Forest Edge"),
  label = paste0("p = ",format.pval(pvals[c("TCC_5km_z","bio15_z","nit_z","ForestEdge_30m_z")],digits = 3)))
pval_df <- pval_df %>% mutate(x = -Inf,y = Inf)
r2_value <- as.numeric(r2(selectedMod)[2])
r2_df <- data.frame(
  variable = unique(preds$variable),
  label = paste0("Full Model R² = ", round(r2_value, 3)),
  x = -Inf,
  y = Inf
)
regplot1<-ggplot(preds, aes(x, predicted)) +
  geom_line(linewidth = .5) +
  facet_wrap(~variable, scales = "free_x") +
  geom_text(
    data = pval_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.73,
    vjust = 1.6,
    size = 2,
    fontface = 'italic'
  ) +
  geom_text(
    data = r2_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.45,
    vjust = 3.0,
    size = 2,
    fontface = "italic"
  ) +
  labs(x = "Scaled variables",y = "Probability of being a hybrid") + 
  theme_bw()
regplot1

hist(pops_health8$ForestEdge_30m)
hist(pops_health8$ForestEdge_30m[pops_health8$HybridYN==1])
hist(pops_health8$ForestEdge_30m[pops_health8$HybridYN==0])

#plot untransformed predictions
unscale <- function(z, original) {z * sd(original, na.rm = TRUE) + mean(original, na.rm = TRUE)}
pred1 <- ggpredict(selectedMod, terms = "TCC_5km_z") %>%
  mutate(
    x = unscale(x, pops_health8$TCC_5km),
    variable = "TCC")
pred2 <- ggpredict(selectedMod, terms = "bio15_z") %>%
  mutate(
    x = unscale(x, pops_health8$bio15),
    variable = "BIO15")
pred3 <- ggpredict(selectedMod, terms = "nit_z") %>%
  mutate(
    x = unscale(x, pops_health8$nitrogen_0.5cm),
    variable = "Nitrogen")
pred4 <- ggpredict(selectedMod, terms = "ForestEdge_30m_z") %>%
  mutate(
    x = unscale(x, pops_health8$ForestEdge_30m),
    variable = "Distance to Forest Edge")
preds <- bind_rows(pred1, pred2, pred3, pred4)
preds <- preds |> 
  dplyr::filter(!is.na(variable))
pvals<-summary(selectedMod)$coefficients$cond[,'Pr(>|z|)']
pval_df <- data.frame(
  variable = c("TCC","BIO15","Nitrogen","Distance to Forest Edge"),
  label = paste0("p = ",format.pval(pvals[c("TCC_5km_z","bio15_z","nit_z","ForestEdge_30m_z")],digits = 3)))
pval_df <- pval_df %>% mutate(x = -Inf,y = Inf)
r2_value <- as.numeric(r2(selectedMod)[2])
r2_df <- data.frame(
  variable = unique(preds$variable),
  label = paste0("Full Model R² = ", round(r2_value, 3)),
  x = -Inf,
  y = Inf
)
regplot1B<-ggplot(preds, aes(x, predicted)) +
  geom_line(linewidth = .5) +
  facet_wrap(~variable, scales = "free_x") +
  geom_text(
    data = pval_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.73,
    vjust = 1.6,
    size = 2,
    fontface = 'italic'
  ) +
  geom_text(
    data = r2_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.45,
    vjust = 3.0,
    size = 2,
    fontface = "italic"
  ) +
  labs(x = "Unscaled variables",y = "Probability of being a hybrid") + 
  theme_bw()
regplot1B

ggsave("./summaries/HybridLandscapeModel.png",plot=regplot1B, width=5, height=4, units='in', bg='white')
pred_NLCD <- ggpredict(selectedMod, terms = "NLCD_5km")
levels(as.factor(pred_NLCD$x))
pvals2<-summary(selectedMod)$coefficients$cond[,'Pr(>|z|)']

pval_df2 <- data.frame(
  x = c("agriculture","developed", "other", "wetland"),
  label = paste0(
    "p = ",
    format.pval(
      pvals2[c(
        "NLCD_5kmagriculture",
        "NLCD_5kmdeveloped",
        "NLCD_5kmother",
        "NLCD_5kmwetland"
      )],
      digits = 3
    )
  ),
  y = pred_NLCD$conf.high[match(
    c("agriculture","developed", "other", "wetland"),
    pred_NLCD$x
  )] + 0.04
)
regplot2<-ggplot(pred_NLCD, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.1
  ) +
  geom_text(
    data = pval_df2,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2,
    fontface = "italic"
  ) +
  theme_bw() +
  labs(
    x = "Land cover (NLCD)",
    y = "Probability of being a hybrid"
  )
library(patchwork)
combined_plot <- regplot1B / regplot2
ggsave("./summaries/HybridLandscapeModelplots.png",plot=combined_plot, width=4, height=6, units='in', bg='white')

summary(mods$mod_beta15)
# odds of being a hybrid increases with increasing nitrogen (p = 0.00744**), decreasing precipitation seasonality (p = 4.14e-06***), decreasing distance from a forest edge (p = 0.01800*), decreasing average tree canopy cover (5km focal area) (p = 0.01778*).  Conditional R2: 0.615, Marginal R2: 0.089 (very low, but maybe acceptable for landscape scale studies)


par(mar=c(8,8,8,8))
boxplot(pops_health8$human_footprint_z~pops_health8$NewHyb_FinalAssignment, las=2)
boxplot(pops_health8$human_footprint_z~pops_health8$HybridYN, las=2)
boxplot(pops_health8$TCC_5km~pops_health8$HybridYN, las=2)
boxplot(pops_health8$TCC_1km~pops_health8$NewHyb_FinalAssignment, las=2)

library("DHARMa")
sim <- simulateResiduals(mods$mod_beta6)
plot(sim)
testDispersion(sim) # no dispersion problem


##################################################################
##################################################################
# health analysis
health_vars <- pops_health8 %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))
colnames(pops_health8)
#ancestry <- select(pops_health8, "HybridYN")
#ancestry$HybridYN<-as.factor(ancestry$HybridYN)

ancestry <- select(pops_health8, "Hybrid_Grouping4")
ancestry$HybridYN<-as.factor(ancestry$Hybrid_Grouping4)

envdata<-select(pops_health8,"ForestEdge_NorAmer_custom_z","ph_z","nit_z","PlantHeight_ft", "bio15_z", "TCC_z")
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-pops_health8[complete_rows,]
dim(ancestry_clean)
health_vars_scaled <- scale(health_vars_clean)
envdata_scaled<-scale(envdata_clean) %>% as.data.frame
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
ind_scores_df <- as.data.frame(ind_scores)
ind_scores_df$GroupFactor <- ids$JC_struc

# how many hybrids are there in this dataset?
h <- hist(ind_scores_df$GroupFactor, plot = FALSE)
data.frame(Lower = head(h$breaks, -1),Upper = tail(h$breaks, -1),Count = h$counts)
## not very many...

ind_scores_df$GroupFactor <- as.factor(ids$HybridYN)
env_arrows <- scores(hyb_rda, display = "bp", scaling = 3)[, 1:2]
env_arrows_df <- as.data.frame(env_arrows)
env_arrows_df$Variable <- rownames(env_arrows_df)

# rescale arrows to fit the RDA plot scale
arrow_multiplier <- 2
env_arrows_df$RDA1 <- env_arrows_df$RDA1 * arrow_multiplier
env_arrows_df$RDA2 <- env_arrows_df$RDA2 * arrow_multiplier
summ<-summary(hyb_rda)
percexpl_rda1<-round((summ$concont$importance[2,1])*100,2)
percexpl_rda2<-round((summ$concont$importance[2,2])*100,2)
# 3. Plot with ggplot2
library(ggplot2)
plot<-ggplot() +
  geom_point(data = ind_scores_df, aes(x = RDA1, y = RDA2, color = GroupFactor), size = 3) +
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
colnames(pops_health8)
pops_health8_adults<-pops_health8 %>% filter(PlantHeight_ft>10)
dim(pops_health8_adults)
health_vars <- pops_health8_adults %>%
  select(InfectionSeverity,
         AreaInfectedByCanker_Trunk_perc,
         AreaInfectedByCanker_RootFlare_prec) %>%
  mutate(across(everything(), as.numeric))
ancestry <- select(pops_health8_adults, "HybridYN")
ancestry$HybridYN<-as.factor(ancestry$HybridYN)
envdata<-select(pops_health8_adults,"ForestEdge_NorAmer_custom_z","ph_z","nit_z","PlantHeight_ft", "bio15_z", "TCC_z","bio6_z")
complete_rows <- complete.cases(health_vars, ancestry, envdata)
sum(complete_rows)

health_vars_clean <- health_vars[complete_rows, ]
ancestry_clean <- ancestry[complete_rows,]
envdata_clean<-envdata[complete_rows,]
ids<-pops_health8_adults[complete_rows,]
dim(ancestry_clean)
health_vars_scaled <- scale(health_vars_clean)
envdata_scaled<-scale(envdata_clean) %>% as.data.frame
colnames(envdata_scaled)
predictors<-cbind(envdata_scaled,ancestry_clean)
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

