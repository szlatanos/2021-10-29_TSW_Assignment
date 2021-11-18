make.matrix <- function(z)
{
  zd <- as.Date(z[,1])
  zz <- as.matrix(z[,2:NCOL(z)])
  zz <- apply(zz, 2, as.numeric)
  rownames(zz) <- as.character(zd)
  colnames(zz) <- colnames(z)[2:NCOL(z)]
  return(zz)
}

remove.NAs <- function(z1, z2)
{
  na1 <- na2 <- {}
  
  for(i in 1:NROW(z1))
  {
    if(sum(is.na(z1[i,]))>0){
      na1 <- c(na1, i)
    }
    if(sum(is.na(z2[i,]))>0){
      na2 <- c(na2, i)
    }
  }
  if(is.null(na1)==TRUE){ na1 <- NA}
  if(is.null(na2)==TRUE){ na2 <- NA}
  
  allna <- na.omit(c(na1, na2))
  allna <- unique(allna)
  
  if(sum(is.na(allna)==TRUE)==0){
    allna <- 0
  }
  
  newX <- z1[-allna,]
  newY <- as.matrix(z2[-allna,])
  return(list(x=newX, y=newY))
}

# Four VAR
remove.NAs.4v <- function(z1, z2, z3, z4)
{
  na1 <- na2 <- na3 <- na4 <- {}
  
  for(i in 1:NROW(z1))
  {
    if(sum(is.na(z1[i,]))>0){
      na1 <- c(na1, i)
    }
    if(sum(is.na(z2[i,]))>0){
      na2 <- c(na2, i)
    }
    if(sum(is.na(z3[i,]))>0){
      na3 <- c(na3, i)
    }
    if(sum(is.na(z4[i,]))>0){
      na4 <- c(na4, i)
    }
  }
  if(is.null(na1)==TRUE){ na1 <- NA}
  if(is.null(na2)==TRUE){ na2 <- NA}
  if(is.null(na3)==TRUE){ na3 <- NA}
  if(is.null(na4)==TRUE){ na4 <- NA}
  
  allna <- na.omit(c(na1, na2, na3, na4))
  allna <- unique(allna)
  
  if(sum(is.na(allna))>0){
    allna <- 0
  }
  
  newX1 <- z1[-allna,]
  newX2 <- z2[-allna,]
  newX3 <- z3[-allna,]
  newY <- as.matrix(z4[-allna,])
  return(list(x1=newX1, x2=newX2, x3=newX3, y=newY))
}
