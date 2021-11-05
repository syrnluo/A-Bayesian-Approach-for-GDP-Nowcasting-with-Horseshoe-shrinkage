#install.packages("MCMCpack")
#install.packages("rARPACK")
#install.packages("pscl")
# install.packages("magic")
# install.packages("MBESS")
#install.packages("ggpubr")

library(dplyr)
#library(tidyr)
library(magic)
library(MASS)
library(coda)
library(MCMCpack)
library(rARPACK)
library(Matrix)
library(Metrics)
#library(MBESS)
library(truncnorm)
library(bayesAB)
#library(matrixNormal)
library(LaplacesDemon)
#library(tidyverse)
library(mltools)
#library(ggpubr)
library(parallel)



###functions###

dbind <- function(m, a, rev=F)
{
  # returns the matrix
  # ( m 0 ) or ( a 0 )
  # ( 0 a ) ( 0 m ) if rev = T
  
  if ((is.vector(m) & length(m) > 1) | (is.vector(a) & length(a) > 1)) stop("Both m and a must be matrix or scalar.")
  
  matrix.index <- function(m, row, col) {(col-1)*nrow(m)+row}
  
  if (!is.matrix(m)) m <- matrix(m, 1, 1)
  if (!is.matrix(a)) a <- matrix(a, 1, 1)
  
  if (rev)
  {
    tmp <- m
    m <- a
    a <- tmp
  }
  
  m.dim <- dim(m)
  a.dim <- dim(a)
  
  if (all(m.dim==0))
  {
    m <- a
  }
  else if (any(a.dim>0))
  {
    new.dim <- m.dim + a.dim
    
    new <- matrix(0, nrow=new.dim[1], ncol=new.dim[2])
    j <- as.vector(matrix.index(new, row(m), col(m)))
    new[j] <- m
    j <- as.vector(matrix.index(new, m.dim[1] + row(a), m.dim[2] + col(a)))
    new[j] <- a
    m <- new
  }
  
  m
} # end of dbind


