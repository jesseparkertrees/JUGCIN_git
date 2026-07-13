#############################################################################################
#####  file generation and analysis (STRUCTURE and NewHybrids) for Hoban microsat data  #####
#############################################################################################

# LOAD PACKAGES
library(dplyr)
library(tidyr)
library(stringr)

# LOAD AND PROCESS MICROSATELLITE METADATA
msat_ids_master<-read.csv("./data/JUGCIN_msatMasterSampleSheet.csv")
colnames(msat_ids_master)
msat_ids_master$SAMPLE.ID
msat_ids_master <- msat_ids_master %>%
  mutate(
    id_num = if_else(
      str_starts(SAMPLE.ID, "JC"),
      str_remove(SAMPLE.ID, "^JC"),
      SAMPLE.ID
    )
  ) %>%
  relocate(id_num, .after = 1)

length(unique(msat_ids_master$id_num))
length(msat_ids_master$id_num)

## check for duplicates
dups<-msat_ids_master[duplicated(msat_ids_master$id_num) |
                  duplicated(msat_ids_master$id_num, fromLast = TRUE), ]
dups

msat_ids_geno<-read.csv("./data/Nuclear_Data_F_5_1_10good.csv")
msat_ids_geno<-msat_ids_geno[,1:2]
head(msat_ids_geno)

popfile<-read.csv("./data/msat_popfile.csv", header=FALSE) %>% rename(pop_name="V1", pop="V2")
head(popfile)

msat_ids_geno2<-left_join(msat_ids_geno,popfile,by="pop")
head(msat_ids_geno2)

sum(msat_ids_master$id_num %in% msat_ids_geno2$ID)
length(msat_ids_geno2$ID)
msat_ids_master$id_num[!msat_ids_master$id_num %in% msat_ids_geno2$ID]
msat_ids_geno2$ID[!msat_ids_geno2$ID %in% msat_ids_master$id_num]

colnames(msat_ids_master)
msat_ids_master<-rename(msat_ids_master,"ID"="id_num")

msat_all<-left_join(msat_ids_geno2, msat_ids_master, by="ID")
length(unique(msat_all$ID))


# GET RESULTS FROM STRUCTURE RUN (performed in GUI and merged with metadata) 
k3qfull<-read.csv("./outputs/structure/msats/50100k_run_7_f_Qsummary_fulldata.csv")

## compare structure results with Hoban's reported structure results from 2010
colnames(k3qfull)
k3qfull$JAdiff_from_hoban
hist(as.numeric(k3qfull$JAdiff_from_hoban), breaks=10)
summary(as.numeric(k3qfull$JCdiff_from_hoban))
hist(as.numeric(k3qfull$JAdiff_from_hoban))
summary(as.numeric(k3qfull$JAdiff_from_hoban))

#############################################################################################
# GENERATE NEWHYBRIDS FILE WITH PURE SPECIES DESIGNATIONS
## (STRUCTURE file (.str) converted to newhybrids format using PGDspider)
## load newhybrids formatted data 
nh <- read.table("./data/msats_newhybs.txt",skip = 5,header = FALSE)
dim(nh)
## load STRUCTURE file 
strucfile<-read.delim("./data/project_data.str", header=FALSE)
dim(strucfile)
## assign IDs to newhybrids file
ids<-strucfile$V1[seq(1, nrow(strucfile), by = 2)]
length(ids)
sum(msat_all$ID %in% ids)
nh$ID<-ids
## remove duplicate and add additional metadata to table (including STRUCTURE results)
nh<-nh[!nh$ID=='JA1147',]
dim(nh)
k3qfull1<-k3qfull[!k3qfull$ID=='JA1147',]
dim(k3qfull1)
sum(duplicated(k3qfull1$ID))
sum(duplicated(nh$ID))
nh2<-left_join(nh, k3qfull1, by="ID")
## add z column indicating "pure" parental species (based on STRUCTURE results) 
nh2$z<- ifelse(nh2$SSR.p.JC..1 > 0.995 & nh2$k1andk2 > 0.995,"z0", ifelse(nh2$SSR.p.JA..1 >= 0.99 & nh2$k3_3 >= 0.99,"z1",NA))
nh2 <- relocate(nh2, "z", .after = "V1")
nh2$V1<-rownames(nh2)
summary(as.factor(nh2$z))
## save newhybrid file IDs and row index for future reference
newhybsIDS<-nh2[,c("ID","z")] 
newhybsIDS$index<-rownames(newhybsIDS)
write.csv(newhybsIDS,"./misc/NewhybridsID_key.csv")
## format table for export to newhybrids
colnames(nh2)
nh3<-nh2[,1:26]
header <- c(
  "NumIndivs 1741",
  "NumLoci 12",
  "Digits 3",
  "Format NonLumped",
  "LocusNames Locus_1 Locus_2 Locus_3 Locus_4 Locus_5 Locus_6 Locus_7 Locus_8 Locus_9 Locus_10 Locus_11 Locus_12"
  )
outfile <- "./data/newhybrids_input_withZ.txt"
## write the header and append genotype table
writeLines(header, outfile)
write.table(nh3,file = outfile,append = TRUE,quote = FALSE,na = "",  sep = " ",row.names = FALSE,col.names = FALSE)

# RUN NEWHYBRIDS IN COMMAND LINE (see: ./misc/msatNEWHYBRIDSparams.txt)

#############################################################################################
# GET NEWHYBRIDS RESULTS
nhq <- read.delim("./outputs/msatNEWHYBRIDS/aa-PofZ.txt", header=TRUE)
nhq$IndivName<-newhybsIDS$ID

nhq<-rename(nhq, c("NewHyb_Butternut"=X1.000.0.000.0.000.0.000,"NewHyb_Heartnut"=X0.000.0.000.0.000.1.000,"NewHyb_F1"=X0.000.0.500.0.500.0.000,"NewHyb_F2"=X0.250.0.250.0.250.0.250,"NewHyb_BC_Butternut"=X0.500.0.250.0.250.0.000,"NewHyb_BC_Heartnut"=X0.000.0.250.0.250.0.500))

## summary of 6 category newhybrids run
colnames(nhq)
sum(nhq$NewHyb_Heartnut>0.75)
sum(nhq$NewHyb_Butternut>0.75)
sum(nhq$NewHyb_F1>0.75)
sum(nhq$NewHyb_F2>0.75)
sum(nhq$NewHyb_BC_Heartnut>0.75)
sum(nhq$NewHyb_BC_Butternut>0.75)

colnames(nhq)
sum(apply(nhq[,3:8] > 0.75, 1, any))
sum(!apply(nhq[,3:8] > 0.75, 1, any))
unassigned <- nhq[!apply(nhq[,3:8] > 0.75, 1, any), ]

length(msat_all$ID)
sum(k3qfull1$ID %in% nhq$IndivName)
nhq<-rename(nhq, "ID"="IndivName")
nhq_all<-left_join(k3qfull1, nhq, by="ID")
getwd()
#write.csv(nhq_all, "nhq_all.csv")
#manually edited nhq_all to include complexF, complexBCJA, and complexBCJC categories

## check out run with 31 categories
nhq31 <- read.delim("./outputs/msatNEWHYBRIDS/31cats/aa-PofZ.txt", header=TRUE)
nhq31$IndivName<-newhybsIDS$ID

## summary of newhybrids runs with 31 categories
nhq31<-rename(nhq31, c(
  "NewHyb_JC"=X1.000.0.000.0.000.0.000,
  "NewHyb_JA"=X0.000.0.000.0.000.1.000,
  "NewHyb_F1"=X0.000.0.500.0.500.0.000,
  "NewHyb_F2"=X0.250.0.250.0.250.0.250,
  "NewHyb_BCJC"=X0.500.0.250.0.250.0.000,
  "NewHyb_BCJA"=X0.000.0.250.0.250.0.500, 
  "NewHyb_BC2JC"=X0.750.0.125.0.125.0.000, 
  "NewHyb_BC2JA"=X0.000.0.125.0.125.0.750, 
  "NewHyb_BC3JC"=X0.875.0.062.0.062.0.000, 
  "NewHyb_BC3JA"=X0.000.0.062.0.062.0.875, 
  "NewHyb_BC4JC"=X0.938.0.031.0.031.0.000, 
  "NewHyb_BC4JA"=X0.000.0.031.0.031.0.938, 
  "NewHyb_BC5JC"=X0.969.0.016.0.016.0.000, 
  "NewHyb_BC5JA"=X0.000.0.016.0.016.0.969, 
  "NewHyb_BC6JC"=X0.984.0.008.0.008.0.000, 
  "NewHyb_BC6JA"=X0.000.0.008.0.008.0.984, 
  "NewHyb_FxJC"=X0.879.0.059.0.059.0.004, 
  "NewHyb_BC2JCxBC3JC"=X0.820.0.086.0.086.0.008,
  "NewHyb_BC2JCxBC2JC"=X0.766.0.109.0.109.0.016, 
  "NewHyb_BCJCxBC3JC"=X0.703.0.141.0.141.0.016, 
  "NewHyb_BCJCxBC2JC"=X0.656.0.156.0.156.0.031, 
  "NewHyb_BCJCxBCJC"=X0.562.0.188.0.188.0.062, 
  "NewHyb_FxJA"=X0.004.0.059.0.059.0.879, 
  "NewHyb_BC2JAxBC3JA"=X0.008.0.086.0.086.0.820, 
  "NewHyb_BC2JAxBC2JA"=X0.016.0.109.0.109.0.766, 
  "NewHyb_BCJAxBC3JA"=X0.016.0.141.0.141.0.703, 
  "NewHyb_BCJAxBC2JA"=X0.031.0.156.0.156.0.656, 
  "NewHyb_BCJAxBCJA"=X0.062.0.188.0.188.0.562, 
  "NewHyb_BCxBC"=X0.059.0.441.0.441.0.059, 
  "NewHyb_BC2JAxBC3JC"=X0.117.0.414.0.414.0.055, 
  "NewHyb_BC3JAxBC2JC"=X0.055.0.414.0.414.0.117))


colnames(nhq31)

summary_table <- nhq31 %>%
  summarise(across(starts_with("NewHyb_"),
                   ~sum(.x > 0.75, na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to = "Category",
               values_to = "N_gt_0.75") %>%
  mutate(Category = sub("^NewHyb_", "", Category))

summary_table

colnames(nhq31)
sum(apply(nhq31[,3:33] > 0.75, 1, any))
sum(!apply(nhq31[,3:33] > 0.75, 1, any))
unassigned <- nhq31[!apply(nhq31[,3:33] > 0.75, 1, any), ]

length(msat_all$ID)
sum(k3qfull1$ID %in% nhq31$IndivName)
nhq31<-rename(nhq, "ID"="IndivName")
nhq31_all<-left_join(k3qfull1, nhq31, by="ID")
getwd()
#write.csv(nhq31_all, "./outputs/msatNEWHYBRIDS/nhq_all.csv")

#############################################################################################
