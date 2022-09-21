
#### Set directory where you have the data files as working dir.. 
WorkDir <- "H:/Minnesota/Teaching/Genomic selection workshops/Plant Breeding Center Data Bootcamp/Data"
setwd(WorkDir)

### Source in functions to be used. This also loads in all the libraries you need.
source("H:/Minnesota/Teaching/Genomic selection workshops/Plant Breeding Center Data Bootcamp/Scripts/bootcamp_functions.R")



### Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "SoyNAM_Geno.vcf"
genotypes_VCF <- read.table(infileVCF)
vcf <- read.vcfR(infileVCF, verbose = FALSE)


### Converting VCF file format to numerical matrix format that can be fit in statistical models
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


### Filtering out markers on proportion of missing data
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.3)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num


### Filtering out individuals on proportion of missing data
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if (length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2 


##Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num3[, -ndx3] else geno_num4 <- geno_num3
  

##Import phenotypic data

pheno <- read.csv("soynam_blups_yld.csv")


##Match phenotypic data to marker data

ndx4 <- match(rownames(geno_num4), pheno$RIL)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

pheno2 <- pheno[ndx5, ]
geno_num5 <- geno_num4[-ndxNA, ]


##Impute using Naive imputation

geno_imp <- replaceNAwithMean(geno_num5)


##Reduce the number of RILs in the dataset simply for the sake of saving time in computation for demonstration (we don't want to spend all of our time watching our computer work!)

ssNdx <- sample.int(n=dim(pheno2)[1], size=5000)
geno_imp_sub <- geno_imp[ssNdx, ]
pheno2_sub <- pheno2[ssNdx, ]


##Fit an RR-BLUP model using the rrBLUP package
rrModel <- mixed.solve(y=pheno2_sub$Yield_blup, Z=geno_imp_sub)
mrkEffsRR <- rrModel$u

##Use marker effects to calculate genomic estimated breeding values of individuals in training set by using . Here we are extracting the intercept and adding it back on.
int <- as.numeric(rrModel$beta)
rrGebv <- int + geno_imp_sub%*%mrkEffsRR


##Calculating G mat using rrBLUP and BGLR and comparing
K <- A.mat(geno_imp_sub)

gblupModel <- kin.blup(data=pheno2_sub, geno='RIL', pheno='Yield_blup', K=K)
gblupGebv <- gblupModel$g


##Compare GEBVs from ridge regression BLUP to G-BLUP
cor(rrGebv, gblupGebv)
plot(rrGebv, gblupGebv)


#############################################################################
#### Start new code for Practical 2 #########################################
#############################################################################

##Use the BGLR package to fit a BayesB model

ETA_BB <- list(list(X=geno_imp_sub, model='BayesB', probIn=.10))

modelBB <- BGLR(y=pheno2_sub$Yield_blup, ETA=ETA_BB, burnIn=500, nIter=2000, verbose=FALSE)

#Pull out GEBVs from model object
bbGebvs <- modelBB$yHat
#Obtain predictions
predict(modelBB)
#Summarize and plot results
summary(modelBB)
plot(modelBB)


##Remove some data to perform a validation analysis
##Use line coding to identify RILs by family. 
fam <- gsub("...$", "", rownames(geno_imp_sub))
ndxFam <- which(fam=="DS11-64")

pheno2_sub_trn <- pheno2_sub

pheno2_sub_trn$Yield_blup[ndxFam] <- NA


ETA_BB <- list(list(X=geno_imp_sub, model='BayesB', probIn=.10))

modelBB <- BGLR(y=pheno2_sub_trn$Yield_blup, ETA=ETA_BB, burnIn=500, nIter=2000, verbose=FALSE)
bbGebvsVal <- modelBB$yHat

##Correlate predictions of RILs left out of the analysis, with predictions
cor(pheno2_sub$Yield_blup[ndxFam],  bbGebvsVal[ndxFam])
plot(pheno2_sub$Yield_blup[ndxFam],  bbGebvsVal[ndxFam])


##Now extend this to perform a 10-fold cross-validation analysis




##This works if my total sample size is divisible by 10. If not, need to subset so it is.
ndxShuf <- sample(1:dim(geno_imp_sub)[1], dim(geno_imp_sub)[1])


pheno_shuf <- pheno2_sub[ndxShuf, ]
geno_imp_sub_shuf <- geno_imp_sub[ndxShuf, ]

cnt <- 1:floor(length(ndxShuf)/10)

predStor <- vector(length=length(ndxShuf))

for (i in 1:10){
  phenoTrn <- pheno_shuf
  phenoTrn[cnt, 2] <- NA
  
  ETA <- list(list(X=geno_imp_sub_shuf, model='BayesB', probIn=.10))
  modelBB <- BGLR(y=phenoTrn$Yield_blup, ETA=ETA, burnIn=500, nIter=2000, verbose=FALSE)
  
  bbGebvs <- modelBB$yHat
  
  predStor[cnt] <- bbGebvs[cnt]
  
  cnt <- cnt + floor(length(ndxShuf)/10)
}

cor(predStor, pheno_shuf[, 2])



















