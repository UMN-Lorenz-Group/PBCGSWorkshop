### Load libraries, functions, and data ----
# Start with a clean working environment by removing all objects stored in R's memory
rm(list=ls(all=TRUE))

# Source in functions to be used. This also loads in all the libraries you need.
source("R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Scripts/bootcamp_functions.R")


# Set directory where you have the data files as working dir.. 
WorkDir <- "R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Data"
setwd(WorkDir)



# Import phenotypes from wheat dataset included as example in BGLR package
data(wheat)
pheno <- wheat.Y
geno <- wheat.X

pheno_means <- rowMeans(pheno)
pheno2 <- cbind(pheno,  pheno_means)
colnames(pheno2) <- c("Env1", "Env2", "Env3", "Env4", "Mean")


rownames(geno) <- rownames(pheno2)



##################################
# Start new code for Practical 3 #
##################################

#Set number of principal components for dimension reduction to select TP
pcNum <- 25

G <- wheat.A

##Perform eigenvalue decomposition for clustering
svdSoyNam <- svd(G, nu=pcNum, nv=pcNum)
PCSoyNam <- G%*%svdSoyNam$v
rownames(PCSoyNam)<-rownames(G)


DistMaize <- dist(PCSoyNam)
TreeMaize <- cutree(hclust(DistMaize), k=6)
plot(PCSoyNam[,1], PCSoyNam[,2], col=TreeMaize, pch=as.character(TreeMaize), xlab="pc1", ylab="pc2")

##Select cluster #5 as target set
TrgtPop <- rownames(PCSoyNam)[TreeMaize==5]
length(TrgtPop)

##Define candidates as all others not in target set
Candidates <- setdiff(rownames(PCSoyNam), TrgtPop)


phenoTrgt <- pheno2[match(TrgtPop, rownames(pheno2)), ]
genoTrgt <- geno[match(TrgtPop, rownames(geno)), ]

phenoCand <- pheno2[match(Candidates, rownames(pheno2)), ]
genoCand <- geno[match(Candidates, rownames(geno)), ]


############################################

##Randomly select a population of size N from candidate set to predict target population. Use this loop to repeat this 10 times
N <- 250


paRandStor <- vector(length=10)

for (i in 1:10){
  
  rnd <- sample(1:dim(phenoCand)[1], N)
  phenoCandRnd <- phenoCand[rnd, ]
  genoCandRnd <- genoCand[rnd, ]
  
  ##Train RR-BLUP model using randomly selected training set and predict target population.
  rrModelRnd <- mixed.solve(y=phenoCandRnd[, 1], Z=genoCandRnd)
  
  mrkEffsRR_rnd <- rrModelRnd$u
  
  ##Use marker effects to calculate genomic estimated breeding values of individuals in training set. Here we are extracting the intercept and adding it back on.
  int <- as.numeric(rrModelRnd$beta)
  rrGebvTrgt_Rnd <- int + genoTrgt%*%mrkEffsRR_rnd
  
  paRandStor[i] <- cor(rrGebvTrgt_Rnd, phenoTrgt[, 1])
  
}



##Now, use an optimization algorithm to select the most informative training set
paOptStor <- vector(length=10)

start_time <- Sys.time()


for (i in 1:10){
  trainOpt <- GenAlgForSubsetSelection(P=PCSoyNam, Candidates=rownames(genoCand), Test=rownames(genoTrgt), ntoselect=N, mc.cores=4, niterations=50, errorstat="CDMEAN")
  
  
  ndx <- match(trainOpt$'Solution with rank 1', rownames(phenoCand))
  
  phenoOpt <- phenoCand[ndx, ]
  genoOpt <- genoCand[ndx, ]
  
  
  ##Now estimate marker effects using rrBLUP package and make predictions in test set
  rrModelOpt <- mixed.solve(y=phenoOpt[, 1], Z=genoOpt)
  
  mrkEffsRR_opt <- rrModelOpt$u
  
  ##Use marker effects to calculate genomic estimated breeding values of individuals in training set. Here we are extracting the intercept and adding it back on.
  int <- as.numeric(rrModelOpt$beta)
  rrGebvOpt <- int + genoTrgt%*%mrkEffsRR_opt
  
  paOptStor[i] <- cor(rrGebvOpt, phenoTrgt[, 1])
  print(paste("completed repetition", i, sep=" "))
}

end_time <- Sys.time()

end_time - start_time


mean(paRandStor)
mean(paOptStor)

paMat <- cbind(paRandStor, paOptStor)
colnames(paMat) <- c("Rand TP Pred Acc", "Opt TP Pred Acc")

boxplot(paMat)









