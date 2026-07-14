####################################################
########  SET UP DIRECTORIES  ######################
####################################################
getwd() #should be "JUGCIN_git/"
dir_path <- "../JUGCIN_git_externalFiles/"
# if it doesn't exist, create it
if (!dir.exists(dir_path)) {
  dir.create(dir_path)
}
# this is where you will put the file "Butternut_Landscape_nJNn9.clean.sub3382.vcf.gz" (acquired separately)
# you will also download environmental datasets into this folder later

####################################################
########  BUTTERNUT GENETIC ANALYSIS  ##############
####################################################
library(vcfR)
library(dartR)
library(SNPfiltR)
library(adegenet)
library(poppr)
library(pegas)

## load vcf (external file, change path accordingly)
vc<-read.vcfR("../JUGCIN_git_externalFiles/Butternut_Landscape_nJNn9.clean.sub3382.vcf.gz")

## convert to genlight
g<-vcfR2genlight(vc)
g<-gl.compliance.check(g)
g@ind.names

## export .str file for running fastStructure in the command line
gl2faststructure(g, outfile = "../JUGCIN_git_externalFiles/Butternut_Landscape_vcf2str.str", outpath=getwd())
# go to command line and run fastSTRUCTURE. See Part 1 of: "./misc/fastSTRUCTUREparams.txt" for instructions/code for running fastSTRUCTURE

######################################################################
########  IMPORT FASTSTRUCTURE RESULTS AND SELECT BEST K  ############
######################################################################
library(dplyr)
library(stringr)
library(ggplot2)
# folder containing logs
log_dir <- "../JUGCIN_git_externalFiles/structure_output/"  
files <- list.files(log_dir, pattern="AZIZrun.*\\.log$", full.names=TRUE)
# extract K and marginal likelihood
dat <- lapply(files, function(f){
  lines <- readLines(f)
  ml <- lines[grep("Marginal Likelihood", lines)]
  ml <- as.numeric(str_extract(ml, "-?[0-9]+\\.?[0-9]*"))
  k <- str_extract(basename(f), "K[0-9]+")
  k <- as.numeric(str_remove(k, "K"))
  rep <- str_extract(basename(f), "rep[0-9]+")
  rep <- as.numeric(str_remove(rep, "rep"))
  data.frame(file=f, K=k, replicate=rep, marginal_likelihood=ml)
})
dat <- bind_rows(dat)
summary_K <- dat %>%
  group_by(K) %>%
  summarise(
    mean_ML = mean(marginal_likelihood),
    sd_ML = sd(marginal_likelihood),
    n = n()
  )
summary_K
ggplot(dat, aes(K, marginal_likelihood)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun=mean, geom="point", size=4, color="red") +
  theme_classic() +
  labs(title="fastStructure K selection PURE butternut")

best_K <-summary_K$K[which.min(summary_K$sd_ML)]
best_rep <- dat %>%
  filter(K == 2) %>%
  arrange(desc(marginal_likelihood)) %>%
  slice(1)
best_rep
best_file <- best_rep$file
meanQ_file <- gsub("\\.log$", ".meanQ", best_file)
Q <- read.table(meanQ_file)
length(g@ind.names)==length(Q$V1)
names<-indNames(g)
newnames <- sub("\\R1_filtered$", "", basename(names))
length(names)
length(newnames)
head(newnames)
head(names)
g@ind.names<-newnames
length(g@ind.names)
head(g@ind.names)
rownames(Q)<-g@ind.names
Q<-rename(Q,"JA_struc"="V1","JC_struc"="V2")
Q$LabID<-rownames(Q)
## export and merge with coordinate data


# MAKE STRUCTURE FILE FOR RUNNING FASTSTRUCTURE ON PURE INDIVIDUALS
sum(rownames(Q)==indNames(g))
pure_ind<-rownames(Q[Q$JC_struc>0.96,])
length(pure_ind)
g_filt_maf <- gl.filter.maf(g, threshold = 0.05)
g_pure<-gl.keep.ind(g_filt_maf, ind.list = pure_ind)
g_pure<-gl.compliance.check(g_pure)

# EXPORT fastSTRUCTURE file
gl2faststructure(g_pure, outfile = "../JUGCIN_git_externalFiles/pure_2524inds_faststr.str", outpath=getwd())
## run fastSTRUCTURE in the command line (see part 2 of "./misc/fastSTRUCTUREparams.txt")

#######################################################################
######   analyze results of fastStructure for pure individuals ########
#######################################################################
library(dplyr)
library(stringr)
library(ggplot2)
## folder containing logs
log_dir <- "../JUGCIN_git_externalFiles/structure_output/pure/" 
files <- list.files(log_dir, pattern="run1.*\\.log$", full.names=TRUE)
## extract K and marginal likelihood
dat <- lapply(files, function(f){
  lines <- readLines(f)
  ml <- lines[grep("Marginal Likelihood", lines)]
  ml <- as.numeric(str_extract(ml, "-?[0-9]+\\.?[0-9]*"))
  k <- str_extract(basename(f), "K[0-9]+")
  k <- as.numeric(str_remove(k, "K"))
  rep <- str_extract(basename(f), "rep[0-9]+")
  rep <- as.numeric(str_remove(rep, "rep"))
  data.frame(file=f, K=k, replicate=rep, marginal_likelihood=ml)})
dat <- bind_rows(dat)
summary_K <- dat %>%
  group_by(K) %>%
  summarise(
    mean_ML = mean(marginal_likelihood),
    sd_ML = sd(marginal_likelihood),
    n = n())
summary_K
best_K <-summary_K$K[which.max(summary_K$mean_ML)] #k=2, least spread between reps
ggplot(dat, aes(K, marginal_likelihood)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun=mean, geom="point", size=4, color="red") +
  theme_classic() +
  labs(title="fastStructure K selection PURE butternut")
best_rep <- dat %>%
  filter(K == best_K) %>%
  arrange(desc(marginal_likelihood)) %>%
  slice(1)
best_rep
best_file <- best_rep$file
meanQ_file <- gsub("\\.log$", ".meanQ", best_file)
Q_pure <- read.table(meanQ_file)
length(g_pure@ind.names)==length(Q_pure$V1)
head(g_pure@ind.names)
Q_pure$LabID<-g_pure@ind.names
Q_pure <- Q_pure[, c(3,1,2)]
Q_pure<-rename(Q_pure, c("pure_K2_v1_JPfastStruc"="V1","pure_K2_v2_JPfastStruc"="V2"))
#write.csv(Q_pure,"Q_pure_K2_fastStruc_JP.csv")

## add for K=3
best_rep <- dat %>%
  filter(K == 3) %>%
  arrange(desc(marginal_likelihood)) %>%
  slice(1)
best_rep
best_file <- best_rep$file
meanQ_file <- gsub("\\.log$", ".meanQ", best_file)
Q_pure_k3 <- read.table(meanQ_file)
length(g_pure@ind.names)==length(Q_pure_k3$V1)
rownames(Q_pure_k3)<-g_pure@ind.names
Q_pure_k3<-rename(Q_pure_k3, c("pure_K3_v1_JPfastStruc"="V1","pure_K3_v2_JPfastStruc"="V2","pure_K3_v3_JPfastStruc"="V3"))
Q_pure<-cbind(Q_pure,Q_pure_k3)
dim(Q_pure)
## add for K=4
best_rep <- dat %>%
 filter(K == 4) %>%
 arrange(desc(marginal_likelihood)) %>%
 slice(1)
best_rep
best_file <- best_rep$file
meanQ_file <- gsub("\\.log$", ".meanQ", best_file)
Q_pure_k4 <- read.table(meanQ_file)
length(g_pure@ind.names)==length(Q_pure_k4$V1)
rownames(Q_pure_k4)<-g_pure@ind.names
Q_pure_k4<-rename(Q_pure_k4, c("pure_K4_v1_JPfastStruc"="V1","pure_K4_v2_JPfastStruc"="V2","pure_K4_v3_JPfastStruc"="V3","pure_K4_v4_JPfastStruc"="V4"))
Q_pure<-cbind(Q_pure,Q_pure_k4)
dim(Q_pure)
#add for K=5
best_rep <- dat %>%
  filter(K == 5) %>%
  arrange(desc(marginal_likelihood)) %>%
  slice(1)
best_rep
best_file <- best_rep$file
meanQ_file <- gsub("\\.log$", ".meanQ", best_file)
Q_pure_k5 <- read.table(meanQ_file)
length(g_pure@ind.names)==length(Q_pure_k5$V1)
rownames(Q_pure_k5)<-g_pure@ind.names
Q_pure_k5<-rename(Q_pure_k5, c("pure_K5_v1_JPfastStruc"="V1","pure_K5_v2_JPfastStruc"="V2","pure_K5_v3_JPfastStruc"="V3","pure_K5_v4_JPfastStruc"="V4", "pure_K5_v5_JPfastStruc"="V5"))
Q_pure<-cbind(Q_pure,Q_pure_k5)
dim(Q_pure)
Q_all<-left_join(Q,Q_pure, by="LabID")
dim(Q_all)
########################################################################
######  Generate NewHybrids file for all 3382 samples    ###############
########################################################################
## designate pure parents using fastSTRUCTURE results
sum(rownames(Q)==indNames(g))
Q$pure<-ifelse(Q$JC_struc>0.9999, "JC", ifelse(Q$JA_struc>0.9999,"JA",NA)) %>% as.factor()
g@strata<-Q
g@pop<-Q$pure %>% as.factor()
summary(g@pop)

g_filt_bal <- gl.filter.maf(g, threshold = 0.05)

## get highly informative SNPs for NewHybrids
af <- gl.percent.freq(g_filt_bal)
colnames(af)
af$freq_prop <- af$frequency / 100
summary(af$popn)
library(tidyr)
af_wide <- af %>%
  dplyr::select(popn, locus, freq_prop) %>%
  pivot_wider(names_from = popn, values_from = freq_prop)

af_wide$delta_af <- abs(af_wide$JA - af_wide$JC)
diag_loci <- af_wide$locus[af_wide$delta_af > 0.99999]
length(diag_loci)

g_fixed_bal <- gl.keep.loc(g_filt_bal, loc.list = diag_loci)
summary(g_fixed_bal@pop)

## generate 4 subsets of 200 randomly sampled SNPs for running NewHybrids (fails with more SNPs)
set.seed(17)
keep1 <- sample(locNames(g_fixed_bal), 200)
g_fixed_bal1 <- gl.keep.loc(g_fixed_bal, keep1)

set.seed(47)
keep2 <- sample(locNames(g_fixed_bal), 200)
g_fixed_bal2 <- gl.keep.loc(g_fixed_bal, keep2)

set.seed(56)
keep3 <- sample(locNames(g_fixed_bal), 200)
g_fixed_bal3 <- gl.keep.loc(g_fixed_bal, keep3)

set.seed(137)
keep4 <- sample(locNames(g_fixed_bal), 200)
g_fixed_bal4 <- gl.keep.loc(g_fixed_bal, keep4)

set.seed(222)
keep5 <- sample(locNames(g_fixed_bal), 200)
g_fixed_bal5 <- gl.keep.loc(g_fixed_bal, keep5)

## create function to convert genotypes for newhybrids:
gl_to_nh <- function(gl_mat) {
  n_ind <- nrow(gl_mat)
  n_loci <- ncol(gl_mat)
  out <- matrix(0, nrow = n_ind, ncol = n_loci * 2)
  for (l in 1:n_loci) {
    
    g <- gl_mat[, l]
    
    a1 <- ifelse(is.na(g), 0,
                 ifelse(g == 0, 1,
                        ifelse(g == 1, 1,
                               ifelse(g == 2, 2, 0))))
    
    a2 <- ifelse(is.na(g), 0,
                 ifelse(g == 0, 1,
                        ifelse(g == 1, 2,
                               ifelse(g == 2, 2, 0))))
    
    out[, (2*l - 1)] <- a1
    out[, (2*l)] <- a2
  }
  return(out)
}

#iterate through the 5 subsets of 200 randomly sampled SNPs
sets <- list(g_fixed_bal1, g_fixed_bal2, g_fixed_bal3, g_fixed_bal4, g_fixed_bal5)
for (s in seq_along(sets)) {
  g_obj <- sets[[s]]
  set_name <- paste0("set_",  s)
  cat("Processing", set_name, "\n")
  # convert genlight to matrix
  gl_mat <- as.matrix(g_obj)
  if (is.null(gl_mat)) {
    stop(paste("Matrix conversion failed for", set_name))
  }
  n_ind <- nrow(gl_mat)
  n_loci <- ncol(gl_mat)
  cat("Individuals:", n_ind, "Loci:", n_loci, "\n")
  pops <- pop(g_obj)
  # function to convert SNP genotypes
  geno <- gl_to_nh(gl_mat)
  # check dimensions
  if (ncol(geno) != 2 * n_loci) {
    stop(paste("Genotype matrix dimension mismatch in", set_name))
  }
  # check illegal genotypes
  bad <- apply(geno, 1, function(x)
    any(x[c(TRUE, FALSE)] == 0 & x[c(FALSE, TRUE)] != 0 |
          x[c(TRUE, FALSE)] != 0 & x[c(FALSE, TRUE)] == 0))
  
  cat("Illegal genotypes:", sum(bad), "\n")
  # assign z tags
  z <- rep("", n_ind)
  z[!is.na(pops) & pops == "JC"] <- "z0"
  z[!is.na(pops) & pops == "JA"] <- "z1"
  
  # convert to NewHybrids format
  geno_nh <- matrix(0, nrow = n_ind, ncol = n_loci)
  
  for (l in 1:n_loci) {
    
    a1 <- geno[, (2*l - 1)]
    a2 <- geno[, (2*l)]
    
    geno_nh[, l] <- ifelse(a1 == 0 | a2 == 0, 0,
                           a1 * 10 + a2)
  }
  # output file
  outfile <- paste0("./outputs/gbsNEWHYBRIDS/inputs/",set_name,".txt")
  
  con <- file(outfile, "w")
  
  writeLines(paste("NumIndivs", n_ind), con)
  writeLines(paste("NumLoci", n_loci), con)
  writeLines("Digits 1", con)
  writeLines("Format Lumped", con)
  writeLines("", con)
  
  for (ind in 1:n_ind) {
    
    geno_str <- paste(geno_nh[ind, ], collapse = " ")
    
    if (z[ind] == "") {
      line <- paste(ind, geno_str)
    } else {
      line <- paste(ind, z[ind], geno_str)
    }
    
    writeLines(line, con)
  }
  
  close(con)
  
  cat("File written:", outfile, "\n\n")
}

list.files("./outputs/gbsNEWHYBRIDS/inputs")
## check output
readLines("./outputs/gbsNEWHYBRIDS/inputs/set_4.txt")

##################################################################
## RUN NewHybrids in command line (see: ./misc/NEWHYBRIDSparams.txt)


########################################################################
#####  Import files from NewHybrids (all 3382 samples)  ################
########################################################################
ids <- indNames(g_fixed_bal)
par_pops<-pop(g_fixed_bal)
nhq1 <- read.delim("./outputs/gbsNEWHYBRIDS/run1/aa-PofZ.txt", header=TRUE)
nhq1$IndivName<-ids
sum(rownames(Q)==indNames(g_fixed_bal))
sum(rownames(Q)==ids)
nhq1<-cbind(nhq1, Q)
dim(nhq1)
colnames(nhq1)
head(nhq1)
rownames(nhq1)<-ids
nhq1$pop<-par_pops
nhq1<-rename(nhq1, c(
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

# run 2
nhq2 <- read.delim("./outputs/gbsNEWHYBRIDS/run2/aa-PofZ.txt", header=TRUE)
nhq2$IndivName<-ids
sum(rownames(Q)==indNames(g_fixed_bal))
sum(rownames(Q)==ids)
nhq2<-cbind(nhq2, Q)
rownames(nhq2)<-ids
nhq2$pop<-par_pops
nhq2<-rename(nhq2, c(
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

# run 3
nhq3 <- read.delim("./outputs/gbsNEWHYBRIDS/run3/aa-PofZ.txt", header=TRUE)
nhq3$IndivName<-ids
nhq3<-cbind(nhq3, Q)
dim(nhq3)
rownames(nhq3)<-ids
nhq3$pop<-par_pops
nhq3<-rename(nhq3, c(
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

nhq4 <- read.delim("./outputs/gbsNEWHYBRIDS/run4/aa-PofZ.txt", header=TRUE)
nhq4$IndivName<-ids
nhq4<-cbind(nhq4, Q)
dim(nhq4)
rownames(nhq4)<-ids
nhq4$pop<-par_pops
nhq4<-rename(nhq4, c(
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

runs <- list(nhq1, nhq2, nhq3, nhq4)
butternut_mat <- sapply(runs, function(x) x$NewHyb_JC)
heartnut_mat <- sapply(runs, function(x) x$NewHyb_JA)
butternutBC_mat <- sapply(runs, function(x) x$NewHyb_BCJC)
heartnutBC_mat <- sapply(runs, function(x) x$NewHyb_BCJA)
butternutBC2_mat <- sapply(runs, function(x) x$NewHyb_BC2JC)
heartnutBC2_mat <- sapply(runs, function(x) x$NewHyb_BC2JA)
butternutBC3_mat <- sapply(runs, function(x) x$NewHyb_BC3JC)
heartnutBC3_mat <- sapply(runs, function(x) x$NewHyb_BC3JA)
butternutBC4_mat <- sapply(runs, function(x) x$NewHyb_BC4JC)
heartnutBC4_mat <- sapply(runs, function(x) x$NewHyb_BC4JA)
F1_mat <- sapply(runs, function(x) x$NewHyb_F1)
F2_mat <- sapply(runs, function(x) x$NewHyb_F2)

cor(heartnut_mat)
cor(butternut_mat)
cor(butternutBC_mat)
cor(heartnutBC_mat)
cor(butternutBC2_mat)
cor(heartnutBC2_mat)
cor(butternutBC3_mat)
cor(heartnutBC3_mat)
cor(butternutBC4_mat)
cor(heartnutBC4_mat)
cor(F1_mat)
cor(F2_mat)

butternut_sd <- apply(butternut_mat, 1, sd)
summary(butternut_sd)
ids[butternut_sd > 0.1]

classes <- c(
  "NewHyb_JC",
  "NewHyb_JA",
  "NewHyb_F1",
  "NewHyb_F2",
  "NewHyb_BCJC",
  "NewHyb_BCJA", 
  "NewHyb_BC2JC", 
  "NewHyb_BC2JA", 
  "NewHyb_BC3JC", 
  "NewHyb_BC3JA", 
  "NewHyb_BC4JC", 
  "NewHyb_BC4JA", 
  "NewHyb_BC5JC", 
  "NewHyb_BC5JA", 
  "NewHyb_BC6JC", 
  "NewHyb_BC6JA", 
  "NewHyb_FxJC", 
  "NewHyb_BC2JCxBC3JC",
  "NewHyb_BC2JCxBC2JC", 
  "NewHyb_BCJCxBC3JC", 
  "NewHyb_BCJCxBC2JC", 
  "NewHyb_BCJCxBCJC", 
  "NewHyb_FxJA", 
  "NewHyb_BC2JAxBC3JA", 
  "NewHyb_BC2JAxBC2JA", 
  "NewHyb_BCJAxBC3JA", 
  "NewHyb_BCJAxBC2JA", 
  "NewHyb_BCJAxBCJA", 
  "NewHyb_BCxBC", 
  "NewHyb_BC2JAxBC3JC", 
  "NewHyb_BC3JAxBC2JC")
assignments <- lapply(runs, function(df){
  classes[max.col(df[,classes])]
})

assignments <- do.call(cbind, assignments)
threshold<-0.75
assignments_prob <- lapply(runs, function(df) {
  apply(df[, classes], 1, function(x) {
    if(max(x) >= threshold) {
      classes[which.max(x)]
    } else {
      NA  # no confident assignment
    }
  })
})
assignments_prob <- do.call(cbind, assignments_prob)
confident_counts <- apply(assignments_prob, 1, function(x) sum(!is.na(x)))
majority_class_prob <- apply(assignments_prob, 1, function(x) {
  x_conf <- x[!is.na(x)]
  if(length(x_conf) == 0) {
    NA  # no confident assignment in any run
  } else {
    ux <- unique(x_conf)
    ux[which.max(tabulate(match(x_conf, ux)))]
  }
})

called_counts <- mapply(function(row, maj){
  if(is.na(maj)) {
    NA
  } else {
    sum(row == maj, na.rm = TRUE)
  }
}, split(assignments_prob, row(assignments_prob)), majority_class_prob)
summary_df <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_ConfRunsCount = majority_class_prob,
  confident_called_runs = called_counts,      # runs supporting the called class
  confident_runs_total = confident_counts # runs with any class >0.95
)
head(summary_df)
prob_arrays <- lapply(runs, function(df) df[, classes])
prob_sum <- Reduce("+", prob_arrays)

top2 <- t(apply(prob_sum, 1, function(x) {
  ord <- order(x, decreasing = TRUE)
  c(
    classes[ord[1]],
    x[ord[1]],
    classes[ord[2]],
    x[ord[2]]
  )
}))

# build dataframe with proper types
summary_df1 <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_totalProb = top2[,1],
  summed_prob = as.numeric(top2[,2]),
  second_class = top2[,3],
  second_prob = as.numeric(top2[,4]),
  stringsAsFactors = FALSE
)

###
prob_arrays <- lapply(runs, function(df) df[, classes])

top2 <- t(sapply(seq_len(nrow(prob_arrays[[1]])), function(i) {
  
  # numeric matrix: rows = classes, cols = runs
  probs_i <- do.call(cbind, lapply(prob_arrays, function(mat) as.numeric(mat[i, ])))
  
  # summed probabilities (original behavior)
  sum_probs <- rowSums(probs_i)
  ord <- order(sum_probs, decreasing = TRUE)
  
  # extract per-run probabilities (rounded)
  top1_vals <- round(probs_i[ord[1], ], 3)
  top2_vals <- round(probs_i[ord[2], ], 3)
  
  c(
    classes[ord[1]],
    sum_probs[ord[1]],
    paste(top1_vals, collapse = "; "),
    classes[ord[2]],
    sum_probs[ord[2]],
    paste(top2_vals, collapse = "; ")
  )
}))

summary_df1 <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_totalProb = top2[,1],
  summed_prob = as.numeric(top2[,2]),   # summed probability
  majority_class_probs = top2[,3],      # per-run probs (rounded)
  second_class = top2[,4],
  second_prob = as.numeric(top2[,5]),   # summed probability (2nd)
  second_class_probs = top2[,6],        # per-run probs (rounded)
  stringsAsFactors = FALSE
)
###
rownames(summary_df1)<-summary_df1$Indiv
head(summary_df)
sum(rownames(summary_df1)==rownames(summary_df))
summary_df_all<-cbind(summary_df1, summary_df[,2:length(summary_df)])

summary_df_all$agree<-summary_df_all$majority_class_ConfRunsCount==summary_df_all$majority_class_totalProb
summary(as.factor(summary_df_all$majority_class_totalProb))
summary(as.factor(summary_df_all$majority_class_ConfRunsCount))
summary(summary_df_all$summed_prob)

colnames(Q_all)
colnames(summary_df_all)
summary_df_all<-rename(summary_df_all, "LabID"="Indiv")
summary_df_all2<-left_join(summary_df_all,Q_all,by="LabID")
summary_df_all2$NewHyb_FinalAssignment<-ifelse(summary_df_all2$confident_called_runs>2,summary_df_all2$majority_class_totalProb, ifelse(summary_df_all2$JC_struc>0.55, "complexBC_JC", ifelse(summary_df_all2$JC_struc<0.45,"complexBC_JA", "complexF")))
                                                   
write.csv(summary_df_all2, "genetics_summary.csv")


############################################################
############################################################
#archive?
#summarize newhybrids results for 37 categories
ids <- indNames(g_fixed_bal)
par_pops<-pop(g_fixed_bal)
list.files("./outputs/gbsNEWHYBRIDS/")
nhq1 <- read.delim("./outputs/gbsNEWHYBRIDS/37cats_run1/aa-PofZ.txt", header=TRUE)
nhq1$IndivName<-ids
sum(rownames(Q)==indNames(g_fixed_bal))
sum(rownames(Q)==ids)
nhq1<-cbind(nhq1, Q)
dim(nhq1)
colnames(nhq1)
head(nhq1)
rownames(nhq1)<-ids
nhq1$pop<-par_pops
colnames(nh)
nhq1<-rename(nhq1, c(
  "JC"="X1.000.0.000.0.000.0.000", 
  "JA"="X0.000.0.000.0.000.1.000",
  "F1"="X0.000.0.500.0.500.0.000", 
  "F2"="X0.250.0.250.0.250.0.250", 
  "BCJC"="X0.500.0.250.0.250.0.000", 
  "BCJA"="X0.000.0.250.0.250.0.500",
  "BC2JC"="X0.750.0.125.0.125.0.000",
  "BC2JA"="X0.000.0.125.0.125.0.750",
  "BC3JC"="X0.875.0.062.0.062.0.000",
  "BC3JA"="X0.000.0.062.0.062.0.875",
  "BC4JC"="X0.938.0.031.0.031.0.000",
  "BC4JA"="X0.000.0.031.0.031.0.938",
  "BC5JC"="X0.969.0.016.0.016.0.000",
  "BC5JA"="X0.000.0.016.0.016.0.969",
  "BC6JC"="X0.984.0.008.0.008.0.000",
  "BC6JA"="X0.000.0.008.0.008.0.984",
  "BC3JCxBC3JC"="X0.879.0.059.0.059.0.004",
  "BC2JCxBC3JC"="X0.820.0.086.0.086.0.008",
  "BC2JCxBC2JC"="X0.766.0.109.0.109.0.016",
  "BCJCxBC3JC"="X0.703.0.141.0.141.0.016",
  "BCJCxBC2JC"="X0.656.0.156.0.156.0.031",
  "BCJCxBCJC"="X0.562.0.188.0.188.0.062",
  "BC3JAxBC3JA"="X0.004.0.059.0.059.0.879",
  "BC2JAxBC3JA"="X0.008.0.086.0.086.0.820",
  "BC2JAxBC2JA"="X0.016.0.109.0.109.0.766",
  "BCJAxBC3JA"="X0.016.0.141.0.141.0.703",
  "BCJAxBC2JA"="X0.031.0.156.0.156.0.656",
  "BCJAxBCJA"="X0.062.0.188.0.188.0.562",
  "BCJCxBCJA"="X0.059.0.441.0.441.0.059",
  "BC2JAxBC3JC"="X0.117.0.414.0.414.0.055",
  "BC3JAxBC2JC"="X0.055.0.414.0.414.0.117",
  "F1xF2"="X0.125.0.375.0.375.0.125",
  "F2xBCJC"="X0.438.0.312.0.188.0.062",
  "F2xBCJA"="X0.062.0.188.0.312.0.438",
  "F1xBCJC"="X0.250.0.375.0.375.0.000",
  "F1xBCJA"="X0.000.0.375.0.375.0.250",
  "BCJCxBC2JAvs"="X0.281.0.219.0.219.0.281"))

# run 2
nhq2 <- read.delim("./outputs/gbsNEWHYBRIDS/37cats_run2/aa-PofZ.txt", header=TRUE)
nhq2$IndivName<-ids
sum(rownames(Q)==indNames(g_fixed_bal))
sum(rownames(Q)==ids)
nhq2<-cbind(nhq2, Q)
rownames(nhq2)<-ids
nhq2$pop<-par_pops
nhq2<-rename(nhq2, c(
  "JC"="X1.000.0.000.0.000.0.000", 
  "JA"="X0.000.0.000.0.000.1.000",
  "F1"="X0.000.0.500.0.500.0.000", 
  "F2"="X0.250.0.250.0.250.0.250", 
  "BCJC"="X0.500.0.250.0.250.0.000", 
  "BCJA"="X0.000.0.250.0.250.0.500",
  "BC2JC"="X0.750.0.125.0.125.0.000",
  "BC2JA"="X0.000.0.125.0.125.0.750",
  "BC3JC"="X0.875.0.062.0.062.0.000",
  "BC3JA"="X0.000.0.062.0.062.0.875",
  "BC4JC"="X0.938.0.031.0.031.0.000",
  "BC4JA"="X0.000.0.031.0.031.0.938",
  "BC5JC"="X0.969.0.016.0.016.0.000",
  "BC5JA"="X0.000.0.016.0.016.0.969",
  "BC6JC"="X0.984.0.008.0.008.0.000",
  "BC6JA"="X0.000.0.008.0.008.0.984",
  "BC3JCxBC3JC"="X0.879.0.059.0.059.0.004",
  "BC2JCxBC3JC"="X0.820.0.086.0.086.0.008",
  "BC2JCxBC2JC"="X0.766.0.109.0.109.0.016",
  "BCJCxBC3JC"="X0.703.0.141.0.141.0.016",
  "BCJCxBC2JC"="X0.656.0.156.0.156.0.031",
  "BCJCxBCJC"="X0.562.0.188.0.188.0.062",
  "BC3JAxBC3JA"="X0.004.0.059.0.059.0.879",
  "BC2JAxBC3JA"="X0.008.0.086.0.086.0.820",
  "BC2JAxBC2JA"="X0.016.0.109.0.109.0.766",
  "BCJAxBC3JA"="X0.016.0.141.0.141.0.703",
  "BCJAxBC2JA"="X0.031.0.156.0.156.0.656",
  "BCJAxBCJA"="X0.062.0.188.0.188.0.562",
  "BCJCxBCJA"="X0.059.0.441.0.441.0.059",
  "BC2JAxBC3JC"="X0.117.0.414.0.414.0.055",
  "BC3JAxBC2JC"="X0.055.0.414.0.414.0.117",
  "F1xF2"="X0.125.0.375.0.375.0.125",
  "F2xBCJC"="X0.438.0.312.0.188.0.062",
  "F2xBCJA"="X0.062.0.188.0.312.0.438",
  "F1xBCJC"="X0.250.0.375.0.375.0.000",
  "F1xBCJA"="X0.000.0.375.0.375.0.250",
  "BCJCxBC2JAvs"="X0.281.0.219.0.219.0.281"))


# run 3
nhq3 <- read.delim("./outputs/gbsNEWHYBRIDS/37cats_run3/aa-PofZ.txt", header=TRUE)
nhq3$IndivName<-ids
nhq3<-cbind(nhq3, Q)
dim(nhq3)
rownames(nhq3)<-ids
nhq3$pop<-par_pops
nhq3<-rename(nhq3, c(
  "JC"="X1.000.0.000.0.000.0.000", 
  "JA"="X0.000.0.000.0.000.1.000",
  "F1"="X0.000.0.500.0.500.0.000", 
  "F2"="X0.250.0.250.0.250.0.250", 
  "BCJC"="X0.500.0.250.0.250.0.000", 
  "BCJA"="X0.000.0.250.0.250.0.500",
  "BC2JC"="X0.750.0.125.0.125.0.000",
  "BC2JA"="X0.000.0.125.0.125.0.750",
  "BC3JC"="X0.875.0.062.0.062.0.000",
  "BC3JA"="X0.000.0.062.0.062.0.875",
  "BC4JC"="X0.938.0.031.0.031.0.000",
  "BC4JA"="X0.000.0.031.0.031.0.938",
  "BC5JC"="X0.969.0.016.0.016.0.000",
  "BC5JA"="X0.000.0.016.0.016.0.969",
  "BC6JC"="X0.984.0.008.0.008.0.000",
  "BC6JA"="X0.000.0.008.0.008.0.984",
  "BC3JCxBC3JC"="X0.879.0.059.0.059.0.004",
  "BC2JCxBC3JC"="X0.820.0.086.0.086.0.008",
  "BC2JCxBC2JC"="X0.766.0.109.0.109.0.016",
  "BCJCxBC3JC"="X0.703.0.141.0.141.0.016",
  "BCJCxBC2JC"="X0.656.0.156.0.156.0.031",
  "BCJCxBCJC"="X0.562.0.188.0.188.0.062",
  "BC3JAxBC3JA"="X0.004.0.059.0.059.0.879",
  "BC2JAxBC3JA"="X0.008.0.086.0.086.0.820",
  "BC2JAxBC2JA"="X0.016.0.109.0.109.0.766",
  "BCJAxBC3JA"="X0.016.0.141.0.141.0.703",
  "BCJAxBC2JA"="X0.031.0.156.0.156.0.656",
  "BCJAxBCJA"="X0.062.0.188.0.188.0.562",
  "BCJCxBCJA"="X0.059.0.441.0.441.0.059",
  "BC2JAxBC3JC"="X0.117.0.414.0.414.0.055",
  "BC3JAxBC2JC"="X0.055.0.414.0.414.0.117",
  "F1xF2"="X0.125.0.375.0.375.0.125",
  "F2xBCJC"="X0.438.0.312.0.188.0.062",
  "F2xBCJA"="X0.062.0.188.0.312.0.438",
  "F1xBCJC"="X0.250.0.375.0.375.0.000",
  "F1xBCJA"="X0.000.0.375.0.375.0.250",
  "BCJCxBC2JAvs"="X0.281.0.219.0.219.0.281"))


nhq4 <- read.delim("./outputs/gbsNEWHYBRIDS/37cats_run4/aa-PofZ.txt", header=TRUE)
nhq4$IndivName<-ids
nhq4<-cbind(nhq4, Q)
dim(nhq4)
rownames(nhq4)<-ids
nhq4$pop<-par_pops
nhq4<-rename(nhq4, c(
  "JC"="X1.000.0.000.0.000.0.000", 
  "JA"="X0.000.0.000.0.000.1.000",
  "F1"="X0.000.0.500.0.500.0.000", 
  "F2"="X0.250.0.250.0.250.0.250", 
  "BCJC"="X0.500.0.250.0.250.0.000", 
  "BCJA"="X0.000.0.250.0.250.0.500",
  "BC2JC"="X0.750.0.125.0.125.0.000",
  "BC2JA"="X0.000.0.125.0.125.0.750",
  "BC3JC"="X0.875.0.062.0.062.0.000",
  "BC3JA"="X0.000.0.062.0.062.0.875",
  "BC4JC"="X0.938.0.031.0.031.0.000",
  "BC4JA"="X0.000.0.031.0.031.0.938",
  "BC5JC"="X0.969.0.016.0.016.0.000",
  "BC5JA"="X0.000.0.016.0.016.0.969",
  "BC6JC"="X0.984.0.008.0.008.0.000",
  "BC6JA"="X0.000.0.008.0.008.0.984",
  "BC3JCxBC3JC"="X0.879.0.059.0.059.0.004",
  "BC2JCxBC3JC"="X0.820.0.086.0.086.0.008",
  "BC2JCxBC2JC"="X0.766.0.109.0.109.0.016",
  "BCJCxBC3JC"="X0.703.0.141.0.141.0.016",
  "BCJCxBC2JC"="X0.656.0.156.0.156.0.031",
  "BCJCxBCJC"="X0.562.0.188.0.188.0.062",
  "BC3JAxBC3JA"="X0.004.0.059.0.059.0.879",
  "BC2JAxBC3JA"="X0.008.0.086.0.086.0.820",
  "BC2JAxBC2JA"="X0.016.0.109.0.109.0.766",
  "BCJAxBC3JA"="X0.016.0.141.0.141.0.703",
  "BCJAxBC2JA"="X0.031.0.156.0.156.0.656",
  "BCJAxBCJA"="X0.062.0.188.0.188.0.562",
  "BCJCxBCJA"="X0.059.0.441.0.441.0.059",
  "BC2JAxBC3JC"="X0.117.0.414.0.414.0.055",
  "BC3JAxBC2JC"="X0.055.0.414.0.414.0.117",
  "F1xF2"="X0.125.0.375.0.375.0.125",
  "F2xBCJC"="X0.438.0.312.0.188.0.062",
  "F2xBCJA"="X0.062.0.188.0.312.0.438",
  "F1xBCJC"="X0.250.0.375.0.375.0.000",
  "F1xBCJA"="X0.000.0.375.0.375.0.250",
  "BCJCxBC2JAvs"="X0.281.0.219.0.219.0.281"))

nhq5 <- read.delim("./outputs/gbsNEWHYBRIDS/37cats_run5/aa-PofZ.txt", header=TRUE)
nhq5$IndivName<-ids
nhq5<-cbind(nhq5, Q)
dim(nhq5)
rownames(nhq5)<-ids
nhq5$pop<-par_pops
nhq5<-rename(nhq5, c(
  "JC"="X1.000.0.000.0.000.0.000", 
  "JA"="X0.000.0.000.0.000.1.000",
  "F1"="X0.000.0.500.0.500.0.000", 
  "F2"="X0.250.0.250.0.250.0.250", 
  "BCJC"="X0.500.0.250.0.250.0.000", 
  "BCJA"="X0.000.0.250.0.250.0.500",
  "BC2JC"="X0.750.0.125.0.125.0.000",
  "BC2JA"="X0.000.0.125.0.125.0.750",
  "BC3JC"="X0.875.0.062.0.062.0.000",
  "BC3JA"="X0.000.0.062.0.062.0.875",
  "BC4JC"="X0.938.0.031.0.031.0.000",
  "BC4JA"="X0.000.0.031.0.031.0.938",
  "BC5JC"="X0.969.0.016.0.016.0.000",
  "BC5JA"="X0.000.0.016.0.016.0.969",
  "BC6JC"="X0.984.0.008.0.008.0.000",
  "BC6JA"="X0.000.0.008.0.008.0.984",
  "BC3JCxBC3JC"="X0.879.0.059.0.059.0.004",
  "BC2JCxBC3JC"="X0.820.0.086.0.086.0.008",
  "BC2JCxBC2JC"="X0.766.0.109.0.109.0.016",
  "BCJCxBC3JC"="X0.703.0.141.0.141.0.016",
  "BCJCxBC2JC"="X0.656.0.156.0.156.0.031",
  "BCJCxBCJC"="X0.562.0.188.0.188.0.062",
  "BC3JAxBC3JA"="X0.004.0.059.0.059.0.879",
  "BC2JAxBC3JA"="X0.008.0.086.0.086.0.820",
  "BC2JAxBC2JA"="X0.016.0.109.0.109.0.766",
  "BCJAxBC3JA"="X0.016.0.141.0.141.0.703",
  "BCJAxBC2JA"="X0.031.0.156.0.156.0.656",
  "BCJAxBCJA"="X0.062.0.188.0.188.0.562",
  "BCJCxBCJA"="X0.059.0.441.0.441.0.059",
  "BC2JAxBC3JC"="X0.117.0.414.0.414.0.055",
  "BC3JAxBC2JC"="X0.055.0.414.0.414.0.117",
  "F1xF2"="X0.125.0.375.0.375.0.125",
  "F2xBCJC"="X0.438.0.312.0.188.0.062",
  "F2xBCJA"="X0.062.0.188.0.312.0.438",
  "F1xBCJC"="X0.250.0.375.0.375.0.000",
  "F1xBCJA"="X0.000.0.375.0.375.0.250",
  "BCJCxBC2JAvs"="X0.281.0.219.0.219.0.281"))

runs <- list(nhq1, nhq2, nhq3, nhq4, nhq5)

classes <- c(
  "JC",
  "JA",
  "F1",
  "F2",
  "BCJC",
  "BCJA",
  "BC2JC",
  "BC2JA",
  "BC3JC",
  "BC3JA",
  "BC4JC",
  "BC4JA",
  "BC5JC",
  "BC5JA",
  "BC6JC",
  "BC6JA",
  "BC3JCxBC3JC",
  "BC2JCxBC3JC",
  "BC2JCxBC2JC",
  "BCJCxBC3JC",
  "BCJCxBC2JC",
  "BCJCxBCJC",
  "BC3JAxBC3JA",
  "BC2JAxBC3JA",
  "BC2JAxBC2JA",
  "BCJAxBC3JA",
  "BCJAxBC2JA",
  "BCJAxBCJA",
  "BCJCxBCJA",
  "BC2JAxBC3JC",
  "BC3JAxBC2JC",
  "F1xF2",
  "F2xBCJC",
  "F2xBCJA",
  "F1xBCJC",
  "F1xBCJA",
  "BCJCxBC2JAvs")
assignments <- lapply(runs, function(df){
  classes[max.col(df[,classes])]
})

assignments <- do.call(cbind, assignments)
threshold<-0.75
assignments_prob <- lapply(runs, function(df) {
  apply(df[, classes], 1, function(x) {
    if(max(x) >= threshold) {
      classes[which.max(x)]
    } else {
      NA  # no confident assignment
    }
  })
})
assignments_prob <- do.call(cbind, assignments_prob)
confident_counts <- apply(assignments_prob, 1, function(x) sum(!is.na(x)))
majority_class_prob <- apply(assignments_prob, 1, function(x) {
  x_conf <- x[!is.na(x)]
  if(length(x_conf) == 0) {
    NA  # no confident assignment in any run
  } else {
    ux <- unique(x_conf)
    ux[which.max(tabulate(match(x_conf, ux)))]
  }
})

called_counts <- mapply(function(row, maj){
  if(is.na(maj)) {
    NA
  } else {
    sum(row == maj, na.rm = TRUE)
  }
}, split(assignments_prob, row(assignments_prob)), majority_class_prob)
summary_df <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_ConfRunsCount = majority_class_prob,
  confident_called_runs = called_counts,      # runs supporting the called class
  confident_runs_total = confident_counts # runs with any class >0.95
)
head(summary_df)
prob_arrays <- lapply(runs, function(df) df[, classes])
prob_sum <- Reduce("+", prob_arrays)

top2 <- t(apply(prob_sum, 1, function(x) {
  ord <- order(x, decreasing = TRUE)
  c(
    classes[ord[1]],
    x[ord[1]],
    classes[ord[2]],
    x[ord[2]]
  )
}))

# build dataframe with proper types
summary_df1 <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_totalProb = top2[,1],
  summed_prob = as.numeric(top2[,2]),
  second_class = top2[,3],
  second_prob = as.numeric(top2[,4]),
  stringsAsFactors = FALSE
)

###
prob_arrays <- lapply(runs, function(df) df[, classes])

top2 <- t(sapply(seq_len(nrow(prob_arrays[[1]])), function(i) {
  
  # numeric matrix: rows = classes, cols = runs
  probs_i <- do.call(cbind, lapply(prob_arrays, function(mat) as.numeric(mat[i, ])))
  
  # summed probabilities (original behavior)
  sum_probs <- rowSums(probs_i)
  ord <- order(sum_probs, decreasing = TRUE)
  
  # extract per-run probabilities (rounded)
  top1_vals <- round(probs_i[ord[1], ], 3)
  top2_vals <- round(probs_i[ord[2], ], 3)
  
  c(
    classes[ord[1]],
    sum_probs[ord[1]],
    paste(top1_vals, collapse = "; "),
    classes[ord[2]],
    sum_probs[ord[2]],
    paste(top2_vals, collapse = "; ")
  )
}))

summary_df1 <- data.frame(
  Indiv = runs[[1]]$IndivName,
  majority_class_totalProb = top2[,1],
  summed_prob = as.numeric(top2[,2]),   # summed probability
  majority_class_probs = top2[,3],      # per-run probs (rounded)
  second_class = top2[,4],
  second_prob = as.numeric(top2[,5]),   # summed probability (2nd)
  second_class_probs = top2[,6],        # per-run probs (rounded)
  stringsAsFactors = FALSE
)
###
rownames(summary_df1)<-summary_df1$Indiv
head(summary_df)
sum(rownames(summary_df1)==rownames(summary_df))
summary_df_all<-cbind(summary_df1, summary_df[,2:length(summary_df)])

summary_df_all$agree<-summary_df_all$majority_class_ConfRunsCount==summary_df_all$majority_class_totalProb
summary(as.factor(summary_df_all$majority_class_totalProb))
summary(as.factor(summary_df_all$majority_class_ConfRunsCount))
summary(summary_df_all$summed_prob)

colnames(Q_all)
colnames(summary_df_all)
summary_df_all<-rename(summary_df_all, "LabID"="Indiv")
summary_df_all2<-left_join(summary_df_all,Q_all,by="LabID")


summary_df_all2$NewHyb_FinalAssignment<-ifelse(summary_df_all2$confident_called_runs>2,summary_df_all2$majority_class_totalProb, ifelse(summary_df_all2$summed_prob>3,summary_df_all2$majority_class_totalProb, ifelse(summary_df_all2$JC_struc>0.55, "complexBC_JC", ifelse(summary_df_all2$JC_struc<0.45,"complexBC_JA", "complexBCF"))))
summary(as.factor(summary_df_all2$NewHyb_FinalAssignment))
write.csv(summary_df_all2, "genetics_summary_37cat.csv")

##################################################################
# generate expected ancestry proportions from NewHybrids categories
catgens<-read.csv("./summaries/NewHybrids_categories_generations.csv")
catgensnots<-catgens$Notation2[1:38]
ancestry <- function(x) {
  x <- gsub("\\s+", "", x)
  # base cases
  if (x == "JA") return(0)
  if (x == "JC") return(1)
  # remove outer parentheses only if they enclose the whole expression
  if (substr(x,1,1) == "(" && substr(x,nchar(x),nchar(x)) == ")") {
    depth <- 0
    enclosed <- TRUE
    for(i in seq_len(nchar(x))) {
      char <- substr(x,i,i)
      if(char=="(") depth <- depth + 1
      if(char==")") depth <- depth - 1
      # if depth reaches zero before the end, parentheses are not outer
      if(depth == 0 && i < nchar(x)) {
        enclosed <- FALSE
        break
      }
    }
    if(enclosed) {
      x <- substr(x,2,nchar(x)-1)
    }
  }
  # find top-level *
  depth <- 0
  split_pos <- NULL
  for(i in seq_len(nchar(x))) {
    char <- substr(x,i,i)
    if(char=="(") depth <- depth + 1
    if(char==")") depth <- depth - 1
    
    if(char=="*" && depth == 0) {
      split_pos <- i
      break
    }
  }
  if(is.null(split_pos)) {
    stop(paste("Cannot parse:", x))
  }
  left <- substr(x,1,split_pos-1)
  right <- substr(x,split_pos+1,nchar(x))
  return((ancestry(left)+ancestry(right))/2)
}
ances<-sapply(catgensnots, ancestry)
