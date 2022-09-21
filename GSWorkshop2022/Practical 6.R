### Load libraries, functions, and data ----
# Start with a clean working environment by removing all objects stored in R's memory
rm(list=ls(all=TRUE))

# Source in functions to be used. This also loads in all the libraries you need.
source("R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Scripts/bootcamp_functions.R")


# Set directory where you have the data files as working dir.. 
WorkDir <- "R:/cfans_agro_lore0149_lorenzlabresearch/GS_Workshop/Plant Breeding Center Data Bootcamp/Data"
setwd(WorkDir)


# Import phenotypes
pheno <- read.csv("barley_cross_pred_pheno.csv")


# Load marker data
geno <- read.csv('barley_cross_pred_geno.csv', row.names=1)
mrk_names <- scan("barley_cross_pred_geno.csv", what='character', sep=',', nlines=1)[-1]
colnames(geno) <- mrk_names

# Load genetic map for variance prediction below
map <- read.csv(file = "genoForMap2.csv", na.strings = c("NA", "#N/A"))


### Manipulate, format, and impute data ----

# Randomly select 1000 markers to make computations faster
ranNum <- sample.int(dim(geno)[2], 1000)
geno2 <- geno[, ranNum]


# Impute missing data using naive imputation
geno_imp <- replaceNAwithMean(geno2)
geno_imp <- round(geno_imp, 0)
geno_imp[which(geno_imp==0)] <- -1


# Write file back out, and read it back in, setting header as F so that column names end up being first row. 
write.csv(geno_imp, 'geno_imp.csv')
geno_imp2 <- read.csv('geno_imp.csv', header=F)


# Match up markers in map with marker data file
mrkname <- geno_imp2[1, ][-1]
ndx <- match(map$mrk, mrkname)
ndxNa <- which(is.na(ndx))
map2 <- map[-ndxNa, ]
ndxNa2 <- match(mrkname, map2$mrk)
geno_imp3 <- geno_imp2[, -(which(is.na(ndxNa2))+1)]
map3 <- map2[order(map2$chr, map2$pos), -1]
ndxOrd <- match(map3$mrk, geno_imp3[1, ])
geno_imp4 <- geno_imp3[, c(1, ndxOrd)]


### Identify parents execute PopVar function ----

# Identify a set of parents to predict cross combinations for. Use set of 100 arbitrarily from the middle of the set to save computation time and easier data handling.

parents <- pheno$line_name[101:200]

cross_table <- t(combn(parents, 2))
colnames(cross_table) <- c("Par1", "Par2")

cross_table <- as.data.frame(cross_table)


# Call PopVar function, deterministic version

pop_predict_out <- pop.predict2(G.in = geno_imp4, y.in = pheno, map.in = map3, crossing.table = cross_table)

write.csv(pop_predict_out, "pop_predict_out.csv")

