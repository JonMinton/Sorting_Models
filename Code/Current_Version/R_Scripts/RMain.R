###################################################################################################################################
# Housekeeping
###################################################################################################################################

rm(list=ls())
#setwd("X:/PROJECTS/Sorting Models/SortingModelRcpp/")


setwd("~/Google Drive/PROJECTS/Sorting Models/SortingModelRcpp/")

###################################################################################################################################
# Dependencies
###################################################################################################################################

require(compiler)
require("Rcpp")

sourceCpp("CppMain_06.cpp")


###################################################################################################################################
# R Functions
###################################################################################################################################

MakeMovingCostMatrix <- function(
  proprents,
  medval,
  vendCost=0.03,
  moveCostNum=2910, # moving cost numerator
  moveCostDenom=24.556, # moving cost denominator
  moveOutsideCost=528.71 # cost of moving out of country      
){
  
  # Error checking
  if (length(proprents)!=length(medval)) stop("proprents and medval of different lengths")
  N <- length(proprents)
  
  out <- matrix(nrow=N+1, ncol=N+1)
  
  for (i in 1:(N+1)){
    for (j in 1:(N+1)){      
      # If not moving
      if (i==j){
        out[i,j] <- 0
      } 
      # If not moving out of the county
      if ((i != j) & (i != (N+1)) & (j != (N+1))){
        out[i,j] <- moveCostNum + vendCost *(
          (1 - proprents[i]) * medval[i] + (1 - proprents[j]) * medval[j]
        ) / moveCostDenom
      }
      # If moving out of the county
      if ((i != j) & ((i==(N+1)) | (j==(N+1)))){
        out[i,j] <- moveOutsideCost    
      }  
      
    }
  }
  return (out) 
}



MakeMovingCostMatrix <- cmpfun(MakeMovingCostMatrix)


###############
#cccccccccccccccccccccccccccccccc
#c Stay-put % by race:          c
#c    - white_l = 0.4174        c
#c    - black_l = 0.3254        c
#c    - asian_l = 0.2808        c
#c    - hisp_l = 0.2892          c
#c    - other_l = 0.2687        c
#c                              c
#c    - white_m = 0.4082        c
#c    - black_m = 0.3933        c

propstay.white_low <- 0.4174
propstay.black_low <- 0.3254
propstay.asian_low <- 0.2808
propstay.hisp_low <- 0.2892
propstay.other_low <- 0.2687

propstay.white_middle <- 0.4082
propstay.black_middle <- 0.3933

###################################################################################################################################
# Input Data
###################################################################################################################################
Data0 <- read.csv("Data/InputsCombined.csv")
Data1 <- read.csv("Data/stage2_dec_vc_sheet1.csv")
Data2 <- read.csv("Data/stage2_dec_vc_sheet2.csv")

require(reshape)
#Data0 <- rename(Data0, c(area="id"))
names(Data0)[1] <- "id"

Data01 <- merge(Data0, Data1, by="id", all.x=T)
Data012 <- data.frame(Data01, Data2)

countsold <- {
  Data012$white.l00 + Data012$white.m00 + Data012$white.h00 +
    Data012$black.l00 + Data012$black.m00 + Data012$black.h00 + 
    Data012$hisp.l00 +  Data012$hisp.m00 + Data012$hisp.h00 + 
    Data012$asian.l00 + Data012$asian.m00 + Data012$asian.h00
}

countsnew <- {
  Data012$white.l09 + Data012$white.m09 + Data012$white.h09 +
    Data012$black.l09 + Data012$black.m09 + Data012$black.h09 + 
    Data012$hisp.l09 +  Data012$hisp.m09 + Data012$hisp.h09 + 
    Data012$asian.l09 + Data012$asian.m09 + Data012$asian.h09
}

countsold.white_low <- Data012$white.l00
countsnew.white_low <- Data012$white.l09


correction  <- 0.0001

countsold <- countsold + correction
countsnew <- countsnew + correction

proprents = Data012$prent
medval = Data012$medval


############################################################################################################################

proprents <- proprents
medval <- medval
countsnew <- countsnew
countsold <- countsold


movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)

###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################


In <- list(
  data = list(
    stayerprop= 0.4174,
    counts = list(
      oldcounts=countsold,
      newcounts=countsnew      
      )
  ),
  movingcostmatrix = movingcostmatrix,
  utils = list(
    deltas = rep(0, length(countsold) + 1) # Added one for outside region
  ),
  params = list(
    tol_delta = 10^-2,
    tol_mu = 10^-1,
    mu_upper = 0.50,
    mu_lower = 0,
    maxit_outer = 1000,
    maxit_inner = 100000,
    use_logmu=F
    ),
  dbg = list(
    dbgflag= T,
    verboseflag=T

  )
)

rm(movingcostmatrix)

x <- new(Contraction, In)


x$run()

x$extractests() -> outputs

rm(x)




###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################
movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)


In2 <- list(
    data = list(
        stayerprop= 0.4174,
        counts = list(
            oldcounts=countsold,
            newcounts=countsnew      
        )
    ),
    movingcostmatrix = movingcostmatrix,
    utils = list(
        deltas = outputs$deltas
    ),
    params = list(
        tol_delta = 10^-5,
        tol_mu = 10^-3,
        mu_upper = 1.10 * outputs$mu,
        mu_lower = 0.90 * outputs$mu,
        maxit_outer = 100000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

x2 <- new(Contraction, In2)
x2$run()
x2$extractests() -> outputs2

rm(x2)


#save(outputs2, file="Contraction_Mapping_refined__14_04_14.Rdata")

#######################

###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################
movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)


In3 <- list(
    data = list(
        stayerprop= 0.4174,
        counts = list(
            oldcounts=countsold,
            newcounts=countsnew      
        )
    ),
    movingcostmatrix = movingcostmatrix,
    utils = list(
        deltas = outputs2$deltas
    ),
    params = list(
        tol_delta = 10^-8,
        tol_mu = 10^-6,
        mu_upper = 1.5 * outputs2$mu,
        mu_lower = 0.5 * outputs2$mu,
        maxit_outer = 10000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

x3 <- new(Contraction, In3)
x3$run()
x3$extractests() -> outputs3

rm(x3)

# save(outputs, outputs2, outputs3, file="Contraction_Mapping_Three_Stage__2014_04_16.rData")

###########################################################################################
###########################################################################################
 # White lower only


###################################################################################################################################
# Housekeeping
###################################################################################################################################

rm(list=ls())
#setwd("X:/PROJECTS/Sorting Models/SortingModelRcpp/")


setwd("~/Google Drive/PROJECTS/Sorting Models/SortingModelRcpp/")

filename <- paste0("logs/SortingSessionOutput_", Sys.time(), ".txt")

file.create(
    filename
)
sink()
sink(filename, append=TRUE, split=TRUE)
###################################################################################################################################
# Dependencies
###################################################################################################################################

require(compiler)
require("Rcpp")

sourceCpp("CppMain_06.cpp")


###################################################################################################################################
# R Functions
###################################################################################################################################

MakeMovingCostMatrix <- function(
    proprents,
    medval,
    vendCost=0.03,
    moveCostNum=2910, # moving cost numerator
    moveCostDenom=24.556, # moving cost denominator
    moveOutsideCost=528.71 # cost of moving out of country      
){
    
    # Error checking
    if (length(proprents)!=length(medval)) stop("proprents and medval of different lengths")
    N <- length(proprents)
    
    out <- matrix(nrow=N+1, ncol=N+1)
    
    for (i in 1:(N+1)){
        for (j in 1:(N+1)){      
            # If not moving
            if (i==j){
                out[i,j] <- 0
            } 
            # If not moving out of the county
            if ((i != j) & (i != (N+1)) & (j != (N+1))){
                out[i,j] <- moveCostNum + vendCost *(
                    (1 - proprents[i]) * medval[i] + (1 - proprents[j]) * medval[j]
                ) / moveCostDenom
            }
            # If moving out of the county
            if ((i != j) & ((i==(N+1)) | (j==(N+1)))){
                out[i,j] <- moveOutsideCost    
            }  
            
        }
    }
    return (out) 
}



MakeMovingCostMatrix <- cmpfun(MakeMovingCostMatrix)


###############
#cccccccccccccccccccccccccccccccc
#c Stay-put % by race:          c
#c    - white_l = 0.4174        c
#c    - black_l = 0.3254        c
#c    - asian_l = 0.2808        c
#c    - hisp_l = 0.2892          c
#c    - other_l = 0.2687        c
#c                              c
#c    - white_m = 0.4082        c
#c    - black_m = 0.3933        c

propstay.white_low <- 0.4174
propstay.black_low <- 0.3254
propstay.asian_low <- 0.2808
propstay.hisp_low <- 0.2892
propstay.other_low <- 0.2687

propstay.white_middle <- 0.4082
propstay.black_middle <- 0.3933

###################################################################################################################################
# Input Data
###################################################################################################################################
Data0 <- read.csv("Data/InputsCombined.csv")
Data1 <- read.csv("Data/stage2_dec_vc_sheet1.csv")
Data2 <- read.csv("Data/stage2_dec_vc_sheet2.csv")

#require(reshape)
#Data0 <- rename(Data0, c(area="id"))
names(Data0)[1] <- "id"


Data01 <- merge(Data0, Data1, by="id", all.x=T)
Data012 <- data.frame(Data01, Data2)

countsold <- {
    Data012$white.l00 + Data012$white.m00 + Data012$white.h00 +
        Data012$black.l00 + Data012$black.m00 + Data012$black.h00 + 
        Data012$hisp.l00 +  Data012$hisp.m00 + Data012$hisp.h00 + 
        Data012$asian.l00 + Data012$asian.m00 + Data012$asian.h00
}

countsnew <- {
    Data012$white.l09 + Data012$white.m09 + Data012$white.h09 +
        Data012$black.l09 + Data012$black.m09 + Data012$black.h09 + 
        Data012$hisp.l09 +  Data012$hisp.m09 + Data012$hisp.h09 + 
        Data012$asian.l09 + Data012$asian.m09 + Data012$asian.h09
}

countsold.white_low <- Data012$white.l00
countsnew.white_low <- Data012$white.l09


correction  <- 0.0001

countsold <- countsold + correction
countsnew <- countsnew + correction

countsold.white_low <- countsold.white_low + correction
countsnew.white_low <- countsnew.white_low + correction

proprents = Data012$prent
medval = Data012$medval


############################################################################################################################

proprents <- proprents
medval <- medval
countsnew <- countsnew
countsold <- countsold


movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)

###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################


In <- list(
    data = list(
        stayerprop= propstay.white_low,
        counts = list(
            oldcounts=countsold.white_low,
            newcounts=countsnew.white_low      
        )
    ),
    movingcostmatrix = movingcostmatrix,
    utils = list(
        deltas = rep(0, length(countsold) + 1) # Added one for outside region
    ),
    params = list(
        tol_delta = 10^-1,
        tol_mu = 10^-1,
        mu_upper = 2,
        mu_lower = 0,
        maxit_outer = 10000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

#rm(movingcostmatrix)

x <- new(Contraction, In)


x$run()

x$extractests() -> outputs

rm(x)




###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################
#movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)


In2 <- list(
    data = list(
        stayerprop= 0.4174,
        counts = list(
            oldcounts=countsold.white_low,
            newcounts=countsnew.white_low)        
    ),
    movingcostmatrix = movingcostmatrix,
    utils = list(
        deltas = outputs$deltas
    ),
    params = list(
        tol_delta = 10^-8,
        tol_mu = 10^-6,
        mu_upper = 2,
        mu_lower = 0,
        maxit_outer = 100000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

x2 <- new(Contraction, In2)
x2$run()
x2$extractests() -> outputs2

rm(x2)


#save(outputs2, file="Contraction_Mapping_refined__14_04_14.Rdata")

#######################

###################################################################################################################################
# Input structure for Rcpp
###################################################################################################################################
movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)


In3 <- list(
    data = list(
        stayerprop= 0.4174,
        counts = list(
            oldcounts=countsold,
            newcounts=countsnew      
        )
    ),
    movingcostmatrix = movingcostmatrix,
    utils = list(
        deltas = outputs2$deltas
    ),
    params = list(
        tol_delta = 10^-8,
        tol_mu = 10^-6,
        mu_upper = 1.5 * outputs2$mu,
        mu_lower = 0.5 * outputs2$mu,
        maxit_outer = 10000,
        maxit_inner = 100000,
        use_logmu=F
    ),
    dbg = list(
        dbgflag= T,
        verboseflag=T
        
    )
)

x3 <- new(Contraction, In3)
x3$run()
x3$extractests() -> outputs3

rm(x3)

# save(outputs, outputs2, outputs3, file="Contraction_Mapping_Three_Stage__2014_04_16.rData")

####################################################################################################

####################################################################################################
rm(list=ls())
#setwd("X:/PROJECTS/Sorting Models/SortingModelRcpp/")


setwd("~/Google Drive/PROJECTS/Sorting Models/SortingModelRcpp/")
#load(file="Contraction_Mapping_Three_Stage__2014_04_16.rData")
load("Results_White_lower.Rdata")

# Quantile regression in R

require(quantreg)

deltas_to_use <- outputs_final$deltas
attach(Data012)

qreg <- rq(outputs_final$deltas ~ c(vcrime, 0), tau=0.5)

detach(Data012)



