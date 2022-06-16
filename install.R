#install.packages("rlang")
install.packages("dplyr",dependencies=TRUE)
install.packages("rrBLUP")
install.packages("vcfR")

install.packages("NAM")
install.packages("bWGR")
install.packages("STPGA")
install.packages("BGLR")
install.packages("SOMMER")

install.packages("BiocManager")
BiocManager::install(version = "3.14")
 
options(repos = BiocManager::repositories())
library(BiocManager)


install.packages("devtools")
library(devtools)

install.packages("rJava")
#devtools::install_bitbucket(
#		repo = "bucklerlab/rTASSEL",
#		ref = "master",
#		build_vignettes = FALSE
#) 




