### Load libraries, functions, and data ----
# Set directory where you have the data files as working dir.. 
WorkDir <- "R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Data"
setwd(WorkDir)

# Source in functions to be used. This also loads in all the libraries you need. Change the path to point R towards where you saved this script.
source("R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Scripts/bootcamp_functions.R")

# Set type of imputation method for use below, either naive or markov
impMethod <- "naive"


# Import phenotypic data
pheno <- read.csv("SoyNAM_Pheno_all.csv")


# Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "SoyNAM_Geno.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)


### Format, manipulate and filter genotype data ----

# Converting VCF file format to numerical matrix format that can be fit in statistical models
gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
gt2d_num <- apply(gt2d, 2, as.numeric)
#Adding row names back in
rownames(gt2d_num) <- rownames(gt2d)
geno_num <- t(gt2d_num)


# Filtering out markers on proportion of missing data
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.5)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num


# Filtering out individuals on proportion of missing data
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if (length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2 


# Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num3[, -ndx3] else geno_num4 <- geno_num3
  

# Match phenotypic data to marker data
ndx4 <- match(rownames(geno_num4), pheno$RIL)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

pheno2 <- pheno[ndx5, ]
geno_num5 <- geno_num4[-ndxNA, ]


# Impute genotype data using either naive imputation or Markov chain implemented in the NAM package
if (impMethod == "naive") geno_imp <- replaceNAwithMean(geno_num5)
if (impMethod == "markov") geno_imp <- markov(apply(geno_num5[, -1], 2, as.numeric))
if (impMethod == "markov") rownames(geno_imp_tst) <- rownames(geno_num5)

# Reduce the number of RILs in the dataset simply for the sake of saving time in computation for demonstration (we don't want to spend all of our time watching our computer work!)

ssNdx <- sample.int(n=dim(pheno2)[1], size=1000)
geno_imp_sub <- geno_imp[ssNdx, ]
pheno2_sub <- pheno2[ssNdx, ]


### Fit some genomic prediction models to the data ----


# Fit an RR-BLUP model using the rrBLUP package
rrModel <- mixed.solve(y=pheno2_sub$Seedsize, Z=geno_imp_sub)
mrk_effs_RR <- rrModel$u

# Use marker effects to calculate genomic estimated breeding values of individuals in training set by using . Here we are extracting the intercept and adding it back on.
int <- as.numeric(rrModel$beta)
gebv_rr <- int + geno_imp_sub%*%mrk_effs_RR


# Calculating a genomic relationship matrix using rrBLUP and fitting a G-BLUP model
G <- A.mat(geno_imp_sub)

gblupModel <- kin.blup(data=pheno2_sub, geno='RIL', pheno='Yield', K=G)
gblupGebv <- gblupModel$g


cor(gebv_rr, pheno2_sub$Seedsize)

# Compare GEBVs from ridge regression BLUP to G-BLUP
cor(rrGebv, gblupGebv)
plot(rrGebv, gblupGebv)



### Cross-validation analysis ----
# Now extend this to perform a 10-fold cross-validation analysis

# This works if my total sample size is divisible by 10. If not, need to subset so it is.

ndxShuf <- sample(1:dim(geno_imp_sub)[1], dim(geno_imp_sub)[1])

pheno_shuf <- pheno2_sub[ndxShuf, ]
geno_imp_sub_shuf <- geno_imp_sub[ndxShuf, ]

cnt <- 1:floor(length(ndxShuf)/10)

pred_stor <- vector(length=length(ndxShuf))

for (i in 1:10){
  pheno_trn <- pheno_shuf
  pheno_trn$Seedsize[cnt] <- NA
  
  rrModel <- mixed.solve(y=pheno_trn$Seedsize, Z=geno_imp_sub_shuf)
  mrkEffsRR <- rrModel$u
  
  # Use marker effects to calculate genomic estimated breeding values of individuals in training set by using . Here we are extracting the intercept and adding it back on.
  int <- as.numeric(rrModel$beta)
  gebv_rr <- int + geno_imp_sub_shuf%*%mrkEffsRR
  
  
  pred_stor[cnt] <- gebv_rr[cnt]
  
  cnt <- cnt + floor(length(ndxShuf)/10)
}

cor(pred_stor, pheno_shuf$Seedsize)



plot(pred_stor, pheno_shuf$Yield)











