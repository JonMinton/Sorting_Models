# Input vignette: 

# To do 
#   1) [R] Create list object two deep 
#   2) [C++] load the object
#   3) [C++] load elements of object as subobjects of different types

rm(list=ls())
setwd("X:/RcppTest/vignettes/")

require("Rcpp")

sourceCpp("Input.cpp")


RListObject <- list(
  A = c(1, 6, 3, 2),
  B = c("hello", "there"),
  C = list(
    D = c(6, 2),
    E = c(4.2, 5.3)
    )
  )

Load(RListObject)

#LoadVec(RListObject)

b <- new(Base)

b$set(4)
b$get()

b$loadvec(RListObject)



m <- matrix(c(1,2,3,4), nrow=2)

b$loadmatrix(m)

b$displaymatrix()



b$returnlist()






