
#### Set directory where you have the data files as working dir.. 
WorkDir <- "C:/Users/ivanv/Desktop/GS_Workshop"
setwd(WorkDir)

##Source in functions to be used
source("C:/Users/ivanv/Desktop/UMN_GIT/GPSoy/SoyGen2/App/GS_Pipeline_Jan_2022_FnsApp.R")
#source("C:/Users/ivanv/Desktop/GS_Workshop/bootcamp_functions.R")

##Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "SoyNAM_Geno.vcf"
genotypes_VCF <- read.table(infileVCF)
vcf <- read.vcfR(infileVCF, verbose = FALSE)




##Converting VCF file format to numerical matrix format. Final genotype matrix is geno_num

gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
rownames(gt2) <- rownames(gt)

gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))

gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)

gt2d_num <- apply(gt2d,2,as.numeric)
geno_num <- t(gt2d_num)

# write.csv(gt2d_num, "geno_numerical.csv")
# geno_num <- read.csv("geno_numerical.csv", row.names=1)
# colNames <- scan("geno_numerical.csv", what='character', sep = ",", nlines=1)[-1]
# colnames(geno_num) <- colNames


#table(geno_num)/sum(table(geno_num))
# geno_num
# 0        0.5          1 
# 0.46856114 0.06415616 0.46728270

##Filter markers on % missing
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.2)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num

##Filter individuals on % missing
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

 if(length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2


##Filter markers based on MAF
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num2[, -ndx3] else geno_num4 <- geno_num3
  

##Import phenotypic data

#pheno <- read.csv("soynam_blups_yld.csv")

pheno <- read.csv("SoyNAM_Pheno.csv")
geno_num4_x <- cbind(rownames(geno_num4),geno_num4)
Data <- merge(geno_num4_x,pheno,by="strain",all=TRUE)


YldNA_Indices <- which(is.na(Data$yield))
Data_Sub <- Data[-YldNA_Indices,]


pheno_sub <- Data_Sub[,phenoIndices]
geno_num5 <- Data_Sub[uniqueStrainIndices,genoIndices]

genoStrain <- unique(geno_num4_x[,"strain"])
genoStrainIndices <- which(Data_Sub[,"strain"] %in% genoStrain)

uniqueGenoStrainIndices <- 
phenoIndices <- c(1,1430:1444)
genoIndices <- c(1:1429)
pheno_sub <- Data_Sub[genoStrainIndices,phenoIndices]


geno_num4b <- Data_Sub[genoStrainIndices,genoIndices]
uniqueStrainIndices<- which(!duplicated(geno_num4b[,"strain"]))
geno_num5 <- geno_num4b[uniqueStrainIndices,]

## Check 
length(which(pheno_sub[,"strain"] %in% geno_num5[,"strain"]))

##Match phenotypic data to marker data
#ndx4 <- match(rownames(geno_num4), pheno$RIL)
#ndx4 <- match(rownames(geno_num4), pheno$strain)
# ndxNA <- which(is.na(ndx4))
# if(length(ndxNA)>0) ndx5 <- ndx4[-ndxNA] 
# 
# pheno_sub <- pheno[ndx5, ]
# geno_num5 <- geno_num4[-ndxNA, ]

##Impute using Naive imputation

#geno_imp <- replaceNAwithMean(geno_num5)
anyNA(geno_num5)
geno_imp <- markov(apply(geno_num5[,-1],2,as.numeric))
anyNA(geno_imp)

##Fit several models using BGLR

#BayesB
ETA_BB <- list(list(X=geno_imp, model='BayesB', probIn=.10))
#modelBB <- BGLR(y=pheno_sub$Yield_blup, ETA=ETA, burnIn=1000, nIter=2000, verbose=FALSE)
modelBB <- BGLR(y=pheno_sub$yield, ETA=ETA_BB, burnIn=1000, nIter=2000, verbose=FALSE)

bbGebvs <- modelBB$yHat

predict(modelBB)


summary(modelBB)
plot(modelBB)

yldCol <- which(colnames(pheno_sub) %in% "yield")
colnames(pheno_sub)[13] <- "Yield_blup" 

##BayesRR
ETA_BRR <- list(list(X=geno_imp, model='BRR', probIn=.10))
modelBRR <- BGLR(y=pheno_sub$Yield_blup, ETA=ETA_BRR, burnIn=1000, nIter=2000, verbose=FALSE)

predict(modelBRR)

##BayesC
ETA_BC <- list(list(X=geno_imp, model='BayesC', probIn=.10))
modelBC <- BGLR(y=pheno_sub$Yield_blup, ETA=ETA_BC, burnIn=1000, nIter=2000, verbose=FALSE)
predict(modelBC)
plot(predict(modelBB), predict(modelBC))


##Calculating G mat using rrBLUP and BGLR and comparing
K_rr <- A.mat(geno_imp)
X <- scale(geno_imp, center=TRUE, scale=TRUE)

#D <- (as.matrix(dist(X, method='euclidean'))^2)/dim(geno_imp)[2]

##G-BLUP using RKHS model in BGLR
ETA_GB <- list(list(K=K_rr, model='RKHS'))
modelGB <- BGLR(y=pheno_sub$Yield_blup, ETA=ETA_GB, burnIn=2000, nIter=12000, verbose=FALSE)

##Try fitting G-BLUP using rrBLUP package

gblupModel <- kin.blup(data=pheno_sub, geno='RIL', pheno='Yield_blup', K=K_rr)
gblupGebv <- gblupModel$g




##Cross validation

##To do this simply here, I am going to remove eight RILs so we have a sample size evenly divisible by 10
ndxShuf <- sample(1:dim(geno_imp)[1], dim(geno_imp)[1])

pheno_shuf <- pheno_sub[ndxShuf, ]
geno_imp_shuf <- geno_imp[ndxShuf, ]

cnt <- 1:floor(length(ndxShuf)/10)

predVec <- vector(length=length(ndxShuf))

for (i in 1:10){
  phenoTrain <- pheno_shuf
  phenoTrain[cnt, 2] <- NA
  
  ETA <- list(list(X=geno_imp_shuf, model='BayesB', probIn=.10))
  modelBB <- BGLR(y=phenoTrain$Yield, ETA=ETA, burnIn=1000, nIter=2000, verbose=FALSE)
  
  bbGebvs <- modelBB$yHat
  
  predVec[cnt] <- bbGebvs[cnt]
  
  cnt <- cnt + floor(length(ndxShuf)/10)
}

PA_BB <- cor(predVec, pheno_shuf[, 2])

##Training population optimization

##Perform eigenvalue decomposition for clustering
svdSoy <- svd(K_rr, nu=50, nv=50)
PC50Soy <- K_rr%*%svdSoy$v
rownames(PC50Soy) <- rownames(K_rr)
distSoy <- dist(PC50Soy)
treeSoy <- cutree(hclust(distSoy), k=4)
plot(PC50Soy[,1], PC50Soy[,2], col=treeSoy, pch=as.character(treeSoy), xlab="pc1", ylab="pc2")

##Select cluster #3 as target set
trgtPop <- rownames(PC50Soy)[treeSoy==3]
length(trgtPop)

##Define candidates as all others not in target set
candidates <- setdiff(rownames(PC50Soy), trgtPop)

phenoTrgt <- pheno[match(trgtPop, pheno_sub$RIL), ]
genoTrgt <- geno_imp[match(trgtPop, pheno_sub$RIL), ]

phenoCand <- pheno[match(candidates, pheno_sub$RIL), ]
genoCand <- geno_imp[match(candidates, pheno_sub$RIL), ]


##Now, use an optimization algorithm to select the most informative training set
paOptVec <- vector(length=10)

for (i in 1:10){
  train1 <- GenAlgForSubsetSelection(P=PC50Soy,Candidates=candidates, Test=trgtPop, ntoselect=25, mc.cores=4)
  
  
  ndx <- match(train1$'Solution with rank 1', pheno_sub$RIL)
  
  phenoOpt <- pheno_sub[ndx, ]
  genoOpt <- geno_imp[ndx, ]
  
  
##Now estimate marker effects using rrBLUP package and make predictions in test set

  rrModelOpt <- mixed.solve(y=phenoOpt$Yield_blup, Z=genoOpt)
  mrkEffsRR_opt <- rrModelOpt$u
  
##Use marker effects to calculate genomic estimated breeding values of individuals in training set. Here we are extracting the intercept and adding it back on.
  int <- as.numeric(rrModelOpt$beta)
  rrGebvOpt <- int + genoCand%*%mrkEffsRR_opt
  paOptVec[i] <- cor(rrGebvOpt, phenoCand[, 2])
}

### GxE using SOMMER/BMTME

library(sommer)
data(DT_example)
DT <- DT_example
A <- A_example 

ansSingle <- mmer(Yield~1,
                  random= ~ vs(Name, Gu=A),
                  rcov= ~ units,
                  data=DT, verbose = FALSE)
summary(ansSingle)

### CS
E <- diag(length(unique(DT$Env)))
rownames(E) <- colnames(E) <- unique(DT$Env)
EA <- kronecker(E,A, make.dimnames = TRUE)
ansCS <- mmer(Yield~Env,
              random= ~ vs(Name, Gu=A) + vs(Env:Name, Gu=EA),
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ansCS)





### Single Env SoyNAM

library(sommer)
data(DT_example)
DT <- pheno_sub

K_rr <- A.mat(geno_imp)
colnames(K_rr) <-rownames(geno_imp)
rownames(K_rr) <- rownames(geno_imp)
A <- K_rr


ansSingle <- mmer(Yield_blup~1,
                  random= ~ vs(strain, Gu=A),
                  rcov= ~ units,
                  data=DT, verbose = FALSE)
summary(ansSingle)


### Multi-Environ / CS 


yldCol <- which(colnames(pheno_sub) %in% "yield")
colnames(pheno_sub)[13] <- "Yield_blup" 

env_sub <-  names(which(table(pheno_sub[,"environ"])>5100)[1:3])

env_sub_indices <- which(pheno_sub[,"environ"] %in% env_sub)
DT <- pheno_sub[env_sub_indices,]

DT$environ <- as.factor(DT$environ)

rownames(geno_imp) <- geno_num5[,"strain"]
env_geno_sub_indices <- which(rownames(geno_imp) %in% unique(DT[,"strain"]))
geno_imp_sub <- geno_imp[env_geno_sub_indices,]
K_rr <- A.mat(geno_imp_sub)
colnames(K_rr) <-rownames(geno_imp_sub)
rownames(K_rr) <- rownames(geno_imp_sub)
A <- K_rr

### Limit no_of_environments to create EA matrix of reasonable mem size.. 

A_Sub <- A[1:500,1:500]
DT_Sub <- DT[which(DT[,"strain"] %in% rownames(A_Sub)),]

E <- diag(length(unique(DT$environ)))
rownames(E) <- colnames(E) <- unique(DT$environ)
DT_Sub$environ <- as.factor(DT_Sub$environ)



### Same set of strains in each of the environments 

rmStrains <- names(which(table(DT_Sub[,"strain"]) <3))
DT_Sub1 <- DT_Sub[-which(DT_Sub[,"strain"] %in% rmStrains),]

E <- diag(length(unique(DT_Sub1$environ)))
rownames(E) <- colnames(E) <- unique(DT_Sub1$environ)
EA <- kronecker(E,A_Sub1, make.dimnames = TRUE)

ansCS <- mmer(Yield_blup~environ,
              random= ~ vs(strain, Gu=A_Sub1) + vs(environ:strain, Gu=EA),
              rcov= ~ units,
              data=DT_Sub1, verbose = FALSE)
summary(ansCS)


DT_Sub1_Tb <- as_tibble(DT_Sub1)
DT_Sub1_Strain <- DT_Sub1_Tb %>% group_by(strain)
DT_Sub1_Strain_Comb <- (DT_Sub1_Strain) %>% summarise( 
    Yield_blup =mean(Yield_blup,na.rm=TRUE))



DT_Sub1_Tb <- as_tibble(DT_Sub1)
DT_Sub1_env <- DT_Sub1_Tb %>% group_by(environ)
DT_Sub1_env_Comb <- (DT_Sub1_env) %>% summarise( 
  Yield_blup =mean(Yield_blup,na.rm=TRUE))




















  
  
