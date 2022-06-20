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

ssNdx <- sample.int(n=dim(pheno2)[1], size=5000)
geno_imp_sub <- geno_imp[ssNdx, ]
pheno2_sub <- pheno2[ssNdx, ]



# New code for practical 2

### BGLR model fitting ----
# Use the BGLR package to fit various types of models. BRR = Bayesian ridge regression, BL = Bayes LASSO, BayesA, BayesB, BayesC 

# Remove some data to perform a validation analysis
# Use line coding to identify RILs by family. 
fam <- gsub("...$", "", rownames(geno_imp_sub))
ndxFam <- which(fam=="DS11-64")

pheno2_sub_trn <- pheno2_sub

pheno2_sub_trn$Seedsize[ndxFam] <- NA

G <- A.mat(geno_imp_sub)

ETA <- list(list(K=NULL, X=geno_imp_sub, model='BayesB', probIn=.10))

model_bglr <- BGLR(y=pheno2_sub_trn$Seedsize, ETA=ETA, burnIn=500, nIter=2000, verbose=FALSE)
gebv_bglr <- model_bglr$yHat

# Extract marker effect predictions from model object. Try different models, changing the name of the object storing the effect (e.g,, "bhat_brr") and plot them against one another on a scatter plot.

bhat <- model_bglr$ETA[[1]]$b

# Here is a way to make a trace plot
plot(bhat_brr^2, ylab='Estimated squared marker effect', type='o')


##Correlate predictions of RILs left out of the analysis, with predictions
cor(pheno2_sub$Seedsize[ndxFam],  gebv_bglr[ndxFam])
plot(pheno2_sub$Seedsize[ndxFam],  gebv_bglr[ndxFam])

# Fit a multi-kernel model using BGLR to treat some large-effect QTL as fixed effects, and remaining QTL as random effects. QTL here were previously declared significant using a GWAS analysis. SNP positions of QTL were 1926, 829, 683, 678.
qtl <- c(1926, 829, 683, 678)

ETA_mk <- list(list(X=geno_imp_sub[, qtl], model='FIXED', probIn=.10), list(K=G, X=geno_imp_sub[, -qtl], model='RKHS', probIn=.10))

model_bglr_mk <- BGLR(y=pheno2_sub_trn$Seedsize, ETA=ETA_mk, burnIn=500, nIter=2000, verbose=FALSE)

gebv_bglr_mk <- model_bglr_mk$yHat

cor(pheno2_sub$Seedsize[ndxFam],  gebv_bglr_mk[ndxFam])
plot(pheno2_sub$Seedsize[ndxFam],  gebv_bglr_mk[ndxFam])





