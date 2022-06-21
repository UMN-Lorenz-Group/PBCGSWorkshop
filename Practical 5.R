##### Practical 5 - Modeling GxE in Genomic Prediction

WorkDir <- getwd()
setwd(WorkDir)

##Source in functions to be used
source("R_Functions/GS_Pipeline_Jan_2022_FnsApp.R")
gc()



##Load in genotype data. Use package vcfR to read in and work with vcf file.
infileVCF <- "Data/SoyNAM_Geno.vcf"
genotypes_VCF <- read.table(infileVCF)
vcf <- read.vcfR(infileVCF, verbose = FALSE)
vcf

gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
rownames(gt2) <- rownames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)

gt2d_num<- apply(gt2d,2,as.numeric)
rownames(gt2d_num)<- rownames(gt2d)
geno_num <- t(gt2d_num)
dim(geno_num)
rm(list=grep("gt2",ls(),value=TRUE))


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
  
dim(geno_num4)


pheno <- read.csv("Data/SoyNAM_Pheno.csv")

geno_num4_x <- cbind(rownames(geno_num4),geno_num4)
colnames(geno_num4_x)[1]<- "strain"

### Check strain names have same format in pheno and geno 
pheno[,1] <- gsub("[-.]","",pheno[,1])
geno_num4_x[,1] <- gsub("[-.]","",geno_num4_x[,1])

## Merge Geno and Pheno Data
Data <- merge(geno_num4_x,pheno,by="strain",all=TRUE)

## Remove with missing yiled_blup values 

YldNA_Indices <- which(is.na(Data$yield))
if(length(YldNA_Indices) >0){Data_Sub <- Data[-YldNA_Indices,]}else{Data_Sub <- Data}


genoStrain <- unique(as.character(geno_num4_x[,"strain"]))

genoStrainIndices <- which(Data_Sub[,"strain"] %in% genoStrain)
length(genoStrainIndices)
genoIndices <- grep("ss",colnames(geno_num4_x))
initGenoIndx <- genoIndices[1]
finalGenoIndx <- genoIndices[length(genoIndices)]
phenoIndices <- c(1,c((finalGenoIndx+1):ncol(Data_Sub)))

pheno_sub <- Data_Sub[genoStrainIndices,phenoIndices]
geno_num4b <- Data_Sub[genoStrainIndices,c(1,genoIndices)]


uniqueStrainIndices<- which(!duplicated(geno_num4b[,"strain"]))

if(length(uniqueStrainIndices)>0) {geno_num5 <- geno_num4b[uniqueStrainIndices,]}else{geno_num5 <- geno_num4b}

dim(geno_num5)

rm(geno_num4b)
rm(geno_num4)
rm(geno_num3)
rm(geno_num2)

### set 'yield' colname to 'Yield_blup'

yldCol <- which(colnames(pheno_sub) %in% "yield")
colnames(pheno_sub)[yldCol] <- "Yield_blup" 



### Select 3 environs with largest number of evaluations (lines)  

env_sub <-  names(which(table(pheno_sub[,"environ"])>5100)[1:3])

env_sub_indices <- which(pheno_sub[,"environ"] %in% env_sub)

## Subset Data and Geno tables 
DT <- pheno_sub[env_sub_indices,]

DT$environ <- as.factor(DT$environ)

dim(DT)

#### Impute genotable using markov function from 'NAM' package 

geno_imp <- markov(apply(geno_num5[,-1],2,as.numeric))
rownames(geno_imp) <- geno_num5[,"strain"]
dim(geno_imp)

### 
env_geno_sub_indices <- which(rownames(geno_imp) %in% unique(DT[,"strain"]))
geno_imp_sub <- geno_imp[env_geno_sub_indices,]

dim(geno_imp_sub)

K_rr <- A.mat(geno_imp_sub)
colnames(K_rr) <-rownames(geno_imp_sub)
rownames(K_rr) <- rownames(geno_imp_sub)
A <- K_rr
dim(A)



  
A_Sub <- A[1:500,1:500]
DT_Sub <- DT[which(DT[,"strain"] %in% rownames(A_Sub)),]

E <- diag(length(unique(DT$environ)))
rownames(E) <- colnames(E) <- unique(DT$environ)
dim(E)

### Same set of strains in each of the environments 

rmStrains <- names(which(table(DT_Sub[,"strain"]) <3))
DT_Sub1 <- DT_Sub[-which(DT_Sub[,"strain"] %in% rmStrains),]

A_Sub1 <- A_Sub[-which(rownames(A_Sub) %in% rmStrains),-which(rownames(A_Sub) %in% rmStrains)]
dim(A_Sub1)

#### Exercise 1 : 


#### 1a) Main Effects Model 

fitMain <- mmer(Yield_blup~environ-1,
                random=~vs(strain,Gu=A_Sub1),
                rcov=~units,
                data=DT_Sub1,verbose=FALSE)
summary(fitMain)



m <- model.matrix(~ environ-1 ,data=DT_Sub1)
m_beta <- m %*% as.numeric(fitMain$Beta[,3]) 
PredMain <- m_beta+fitMain$U$`u:strain`$Yield_blup
cor(PredMain,DT_Sub1[,"Yield_blup"]) 


#### 1b) Compound Symmetry Var-Covar 
E <- diag(length(unique(DT_Sub1$environ)))
rownames(E) <- colnames(E) <- unique(DT_Sub1$environ)

EA <- kronecker(E,A_Sub1, make.dimnames = TRUE)
DT_Sub1$environ <- as.factor(DT_Sub1$environ)
DT_Sub1$strain <- as.factor(DT_Sub1$strain)

fitCS <- mmer(Yield_blup~environ-1,
              random= ~ vs(strain, Gu=A_Sub1) + vs(environ:strain, Gu=EA),
              rcov= ~ units,
              data=DT_Sub1, verbose = FALSE)
summary(fitCS)


m <- model.matrix(~ environ-1 ,data=DT_Sub1)
m_beta <- m %*% as.numeric(fitCS$Beta[,3]) 
PredCS <- m_beta+fitCS$U$`u:environ:strain`$Yield_blup
cor(PredCS,DT_Sub1[,"Yield_blup"]) 


#### 1c) Compound Symmetry Var-Covar with heterogeneous variance


fitCSDG <- mmer(Yield_blup~environ-1,
                random=~vs(strain,Gu=A_Sub1) +vs(ds(environ),strain,Gu=A_Sub1),
                rcov=~units,
                data=DT_Sub1,verbose=FALSE) 

summary(fitCSDG)


m2 <- cbind(c(rep(1,nrow(DT_Sub1)/3),rep(0,2*nrow(DT_Sub1)/3)),c(rep(0,nrow(DT_Sub1)/3),rep(1,nrow(DT_Sub1)/3),rep(0,nrow(DT_Sub1)/3)),
c(rep(0,nrow(DT_Sub1)/3),rep(0,nrow(DT_Sub1)/3),rep(1,nrow(DT_Sub1)/3)))

m_beta <- m2 %*% as.numeric(fitCSDG$Beta[,3]) 
length(m_beta)
m_env_strain <- do.call(cbind,lapply(fitCSDG$U,function(x) x$Yield_blup))
dim(m_env_strain)
envStrain_blup <-c(m_env_strain[,2:4])                              
                  
strain_blup <- rep(fitCSDG$U$`u:strain`$Yield_blup,3)
length(strain_blup)

PredCSDG <- m_beta+strain_blup+envStrain_blup

indES <-  sort.int(as.numeric(DT_Sub1[,"environ"]),decreasing=FALSE,index.return=TRUE)[[2]]

cor(PredCSDG,DT_Sub1[indES,"Yield_blup"]) 


#### 1d) Unstructured Var-Covar with heterogeneous environment specific variance and covariance 

fitUS <- mmer(Yield_blup~environ-1,
                random=~vs(us(environ),strain,Gu=A_Sub1),
                rcov=~units,
                data=DT_Sub1,verbose=FALSE) 
summary(fitUS)

envNames <- levels(factor(DT_Sub1$environ))
print(envNames)
print(names(fitUS$U))
env1Ind <- c(1,3,6)
U_envStrain <- list()
PredUS <- list()
  for(i in 1:length(envNames)){
       envInd <-  grep(envNames[i],names(fitUS$U))
       envIndNames <-  grep(envNames[i],names(fitUS$U),value=TRUE)
       U_envStrain[[i]] <-  as.numeric(fitUS$U[[env1Ind[i]]]$Yield_blup)
       for(j  in 2:length(envInd)){ 
         indJ <- envInd[j]
         b <- cbind(names(fitUS$U[[indJ]]$Yield_blup),fitUS$U[[indJ]]$Yield_blup)
         colnames(b) <- c("strain","Yield_blup")
         b_group <- as_tibble(b) %>% group_by(strain)
         YldBlup_group <- b_group %>% summarise(Yield_blup = sum(as.numeric(Yield_blup)))
         U_envStrain[[i]] <- U_envStrain[[i]] +YldBlup_group[,2]
       }
      
      PredUS[[i]] <- c(U_envStrain[[i]] + fitUS$Beta[i,3])
     }

indES <-  sort.int(as.numeric(DT_Sub1[,"environ"]),decreasing=FALSE,index.return=TRUE)[[2]]  
cor(unlist(PredUS),DT_Sub1[indES,"Yield_blup"]) 


##### Exercise 2: Predict GEBVs of untested and tested genotypes in untested and tested environments

### 2a) Tested genotypes in untested environments

### Remove lines from IA2013 and train the model using IA2012 and IL2013 only and predict 
### performance of lines for IA2013 (untested environ) and compare accuracy with model 
### incorporating data from IA2013 in the training model 

tstIndices1 <- which(DT_Sub1[,"environ"] %in% "IA_2013") 

DT_Sub1A <- DT_Sub1
DT_Sub1A[tstIndices1 ,"Yield_blup"] <- NA
#DT_Sub1A[tstIndices1 ,"environ"] <- NA

dim(DT_Sub1A)


E <- diag(length(unique(DT_Sub1A$environ)))
rownames(E) <- colnames(E) <- unique(DT_Sub1A$environ)
dim(E)

EA <- kronecker(E,A_Sub1, make.dimnames = TRUE)
DT_Sub1$Aenviron <- as.factor(DT_Sub1A$environ)
DT_Sub1A$strain <- as.factor(DT_Sub1A$strain) 

dim(EA)

fitUS1A <- mmer(Yield_blup~environ-1,
                random=~vs(us(environ),strain,Gu=A_Sub1),
                rcov=~units,
                data=DT_Sub1A,verbose=FALSE) 
summary(fitUS1A)

names(fitUS1A$U)
envNames <- levels(factor(DT_Sub1A$environ))
env1Ind <- c(1,3,6)
U_envStrain <- list()
PredUS1A <- list()
  for(i in 1:length(envNames)){
       envInd <-  grep(envNames[i],names(fitUS1A$U))
       U_envStrain[[i]] <-  as.numeric(fitUS1A$U[[env1Ind[i]]]$Yield_blup)
     for(j  in 2:length(envInd)){ 
         indJ <- envInd[j]
         b <- cbind(names(fitUS1A$U[[indJ]]$Yield_blup),fitUS1A$U[[indJ]]$Yield_blup)
         colnames(b) <- c("strain","Yield_blup")
         b_group <- as_tibble(b) %>% group_by(strain)
         YldBlup_group <- b_group %>% summarise(Yield_blup = sum(as.numeric(Yield_blup)))
         U_envStrain[[i]] <- U_envStrain[[i]] + YldBlup_group[,2] 
       } 
        PredUS1A[[i]] <- c(U_envStrain[[i]] + fitUS$Beta[i,3])
     }
    
length(unlist(PredUS1A[[2]]))

length(tstIndices1)

cor(unlist(PredUS1A[[2]]),DT_Sub1[tstIndices1,"Yield_blup"]) 


#### 2b) Untested genotypes in tested environments

### Subset Data to generate untested genotypes 

set.seed(125)
tstStrain <- sample(unique(DT_Sub1[,"strain"]),0.2*length(unique(DT_Sub1[,"strain"])))
length(tstStrain)
tstIndices2 <- which(DT_Sub1[,"strain"] %in% tstStrain)
DT_Sub1B <- DT_Sub1
DT_Sub1B[tstIndices2 ,"Yeild_blup"] <- NA
dim(DT_Sub1B)

fitCSDG1B <- mmer(Yield_blup~environ-1,
                random=~vs(strain,Gu=A_Sub1) +vs(ds(environ),strain,Gu=A_Sub1),
                rcov=~units,
                data=DT_Sub1B,verbose=FALSE) 

summary(fitCSDG1B)

m2 <- cbind(c(rep(1,nrow(DT_Sub1B)/3),rep(0,2*nrow(DT_Sub1B)/3)),c(rep(0,nrow(DT_Sub1B)/3),rep(1,nrow(DT_Sub1B)/3),rep(0,nrow(DT_Sub1B)/3)),
c(rep(0,nrow(DT_Sub1B)/3),rep(0,nrow(DT_Sub1B)/3),rep(1,nrow(DT_Sub1B)/3)))

m_beta <- m2 %*% as.numeric(fitCSDG1B$Beta[,3]) 
length(m_beta)
m_env_strain <- do.call(cbind,lapply(fitCSDG1B$U,function(x) x$Yield_blup))
dim(m_env_strain)
envStrain_blup <-c(m_env_strain[,2:4])                              
                  
strain_blup <- rep(fitCSDG1B$U$`u:strain`$Yield_blup,3)
length(strain_blup)

PredCSDG1B<- m_beta+strain_blup+envStrain_blup

indES <-  sort.int(as.numeric(DT_Sub1[,"environ"]),decreasing=FALSE,index.return=TRUE)[[2]]

cor(PredCSDG1B,DT_Sub1[indES,"Yield_blup"]) 

