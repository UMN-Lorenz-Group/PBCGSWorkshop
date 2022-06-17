##Functions for GS bootcamp course



#easy impute, replacing missing with marker mean
replaceNAwithMean <- function(mat){
  replaceNAwithMeanVec <- function(vec){
    mu <- mean(vec, na.rm=TRUE)
    vec[is.na(vec)] <- mu
    return(vec)
  }
  return(apply(mat, 2, replaceNAwithMeanVec))
}