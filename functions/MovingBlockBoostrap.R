# There seems to be some sort of persistence
# therefore, we cannot use iid bootstrap.
# So, we have to use a block bootstrap.
# We will use the Moving Block Bootstrap
#
# We will use a custom function, that we have to load first:
MBB <- function(x, b)
{
  xboot <- NULL
  t <- NROW(x)
  k <- NROW(xboot)
  while(k<=(t+5)){
    length <- b
    point <- runif(1, 1, t-b)
    xstart <- point;
    xend <- point+length
    x_new <- x[(xstart+1):(xend)]
    xboot <- c(xboot, x_new)
    k <- NROW(xboot)
  }
  xboot <- xboot[1:t]
  return(xboot)
}