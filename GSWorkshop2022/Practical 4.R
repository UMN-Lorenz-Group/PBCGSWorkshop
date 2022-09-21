### Load libraries, functions, and data ----
# Set directory where you have the data files as working dir.. 
WorkDir <- "R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Data"
setwd(WorkDir)

# Source in functions to be used. This also loads in all the libraries you need.
source("R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Scripts/bootcamp_functions.R")



# Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "SoyNAM_Geno.vcf"
genotypes_VCF <- read.table(infileVCF)
vcf <- read.vcfR(infileVCF, verbose = FALSE)


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
ndx <- which(mrkNA > 0.3)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num


# Filtering out individuals on proportion of missing data
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if (length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2 


# Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num3[, -ndx3] else geno_num4 <- geno_num3
  

# Import phenotypic data

pheno <- read.csv("SoyNAM_simPheno.csv")


### Manipulate, format, and impute data ----
# Match phenotypic data to marker data

ndx4 <- match(rownames(geno_num4), pheno$RIL)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

pheno2 <- pheno[ndx5, ]
geno_num5 <- geno_num4[-ndxNA, ]


# Impute using Naive imputation

geno_imp <- replaceNAwithMean(geno_num5)


# Reduce the number of RILs in the dataset simply for the sake of saving time in computation for demonstration (we don't want to spend all of our time watching our computer work!)

ssNdx <- sample.int(n=dim(pheno2)[1], size=5000)
geno_imp_sub <- geno_imp[ssNdx, ]
pheno2_sub <- pheno2[ssNdx, ]




# Implement multi-trait genomic prediction models ----

# Use the BGLR package to implement a multi-trait prediction model

ETA_MT_BRR <- list(list(X=geno_imp_sub, model='BRR', probIn=.10))
model_MT_BRR <- Multitrait(y=as.matrix(pheno2_sub[, 2:3]), ETA=ETA_MT_BRR, burnIn=1000, nIter=2000, verbose=FALSE)

gebvs_MT_BRR <- model_MT_BRR$ETAHat



# Multi-trait predictions, SpikeSlab in =BGLR

ETA_MT_SS <- list(list(X=geno_imp_sub, model='SpikeSlab', probIn=.10))
model_MT_SS <- Multitrait(y=as.matrix(pheno2_sub[, 2:3]), ETA=ETA_MT_SS, burnIn=1000, nIter=2000, verbose=FALSE)

gebvs_MT_SS <- model_MT_SS$ETAHat


cor(gebvs_MT_BRR[, 1], pheno2_sub$Trait1)
cor(gebvs_MT_BRR[, 2], pheno2_sub$Trait2)


# Single trait predictions

ETA_BB <- list(list(X=geno_imp_sub, model='BRR', probIn=.10))
modelBB <- BGLR(y=pheno2_sub$Trait1, ETA=ETA_BB, burnIn=500, nIter=2000, verbose=FALSE)
bb_gebvs <- modelBB$yHat

cor(bb_gebvs, pheno2_sub$Trait1)



# Multi-trait prediction using the G-BLUP model execuated in SOMMER

A.Tot <- A.mat(geno_imp_sub)

rownames(A.Tot) <- rownames(geno_imp_sub) 
colnames(A.Tot) <- rownames(geno_imp_sub) 

trait <- c("Yield_blup", "Pro_blup")

fm3 <- mmer(as.formula(paste("cbind(",paste(trait,collapse=","),")~1",sep="")),
            random=~vs(RIL,Gu=A.Tot),
            rcov=~units,
            data=pheno2_sub, verbose = TRUE)

fm3 <- mmer(cbind(Yield_blup, Pro_blup)~1,
            random=~vs(RIL,Gu=A.Tot),
            rcov=~units,
            data=pheno2_sub, verbose = TRUE)



# Perform 10-fold cross-validation analysis and test predictive ability ----

# This works if total sample size is divisible by 10. If not, need to subset so it is.
ndxShuf <- sample(1:dim(geno_imp_sub)[1], dim(geno_imp_sub)[1])


pheno_shuf <- pheno2_sub[ndxShuf, ]
geno_imp_sub_shuf <- geno_imp_sub[ndxShuf, ]

cnt <- 1:floor(length(ndxShuf)/10)

##Testing multi-trait model with cross-validation
pred_mt_stor <- vector(length=length(ndxShuf))

start_time <- Sys.time()

for (i in 1:10){
  pheno_trn <- pheno_shuf
  pheno_trn[cnt, 2] <- NA
  

  ETA_MT_BRR <- list(list(X=geno_imp_sub_shuf, model='BRR', probIn=.10))
  model_MT_BRR <- Multitrait(y=as.matrix(pheno_trn[, 2:3]), ETA=ETA_MT_BRR, burnIn=1000, nIter=2000, verbose=FALSE)
  
  
  gebvs_MT_BRR <- model_MT_BRR$ETAHat
  
  pred_mt_stor[cnt] <- gebvs_MT_BRR[cnt, 2]
  
  cnt <- cnt + floor(length(ndxShuf)/10)
}

end_time <- Sys.time()
end_time - start_time



# Test single-trait predictions for point of comparison

pred_st_stor <- vector(length=length(ndxShuf))
cnt <- 1:floor(length(ndxShuf)/10)

start_time <- Sys.time()

for (i in 1:10){
  pheno_trn <- pheno_shuf
  pheno_trn[cnt, 2] <- NA
  
  
  ETA_ST_BRR <- list(list(X=geno_imp_sub_shuf, model='BRR', probIn=.10))
  model_ST_BRR <- BGLR(y=pheno_trn[, 2], ETA=ETA_ST_BRR, burnIn=1000, nIter=2000, verbose=FALSE)
  
  
  gebvs_ST_BRR <- model_ST_BRR$yHat
  
  pred_st_stor[cnt] <- gebvs_ST_BRR[cnt]
  
  cnt <- cnt + floor(length(ndxShuf)/10)
}

end_time <- Sys.time()
end_time - start_time


# Estimate correlation between observed yield phenotypes and predicted yield from cross-validation

cor(pred_stor, pheno_shuf[, 3])
plot(pred_stor, pheno_shuf[, 3])

















