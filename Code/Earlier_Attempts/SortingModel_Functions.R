# Functions
require(compiler)

CalcMovingCostR <- function(fromArea, toArea,
                           percentRents,
                           medval,
                           vendCost=0.03,
                           moveCostNum=2910, # moving cost numerator
                           moveCostDenom=24.556, # moving cost denominator
                           moveOutsideCost=528.71 # cost of moving out of country                         
){
  
  # PROCESSING
  N <- length(percentRents)  
  
  # If not moving, no cost of move
  if (fromArea==toArea){
    output <- 0
  } 
  # If not moving out of the county
  if ((fromArea != toArea) & (fromArea != N) & (toArea != N)){
    output <- (1 - percentRents[fromArea]) * medval[fromArea] + (1 - percentRents[toArea]) * medval[toArea]
    output <- (moveCostNum + vendCost * output)/moveCostDenom    
  }
  
  # If moving out of the county
  if ((fromArea != toArea) &((fromArea==N) | (toArea==N))){
    output <- moveOutsideCost    
  }  
  
  # ERROR CATCHING
  if (is.na(output)){
    stop("CalcMovingCost: output is NA")
  }
  
  # OUTPUT
  return(output)      
}

CalcMovingCostR <- cmpfun(CalcMovingCostR) # compiles function

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


CalcDenomR <- function(
  areak, 
  deltas, 
  mu, 
  percentRents, 
  medval
){
  
  # PROCESSING 
  output <- 0
  for (l in 1:length(deltas)){
    output <- output + exp(deltas[l] - deltas[areak] - mu * CalcMovingCost(fromArea=areak, toArea=l, percentRents=percentRents, medval=medval))
  }
  
  # ERROR CATCHING
  
  if (is.infinite(output)){
    stop("CalcDenom: output is infinite")
  }
  
  if (output==0){
    stop("CalcDenom: output is zero")
  }
  
  if (is.na(output)){
    stop("CalcDenom: output is NA")
  }
  
  # RETURN 
  return(output)
}

CalcDenomR<- cmpfun(CalcDenomR)



CalcForceAttractR <- function(
  toArea, 
  fromArea, 
  mu, 
  deltas, 
  percentRents, 
  medval,
  verbose=F
  ){
  
  # PROCESSING
  output <- exp(deltas[toArea] - deltas[fromArea] - mu * CalcMovingCost(fromArea=fromArea, toArea=toArea, percentRents, medval=medval))
  
  # ERROR CATCHING
  if(is.infinite(output)){
    stop("CalcForceAttract: output infinite")
  }
  
  if (is.na(output)){
    stop("CalcForceAttract: output NA")
  }  
  
  # OUTPUT
  return(output)
}

CalcForceAttractR <- cmpfun(CalcForceAttractR)


CalcDenomR <- function(
  areak, 
  deltas, 
  mu, 
  percentRents, 
  medval
  ){
  
  # PROCESSING 
  output <- 0
  for (l in 1:length(deltas)){
    output <- output + exp(deltas[l] - deltas[areak] - mu * CalcMovingCost(fromArea=areak, toArea=l, percentRents=percentRents, medval=medval))
  }
  
  # ERROR CATCHING
  
  if (is.infinite(output)){
    stop("CalcDenom: output is infinite")
  }
  
  if (output==0){
    stop("CalcDenom: output is zero")
  }
  
  if (is.na(output)){
    stop("CalcDenom: output is NA")
  }
  
  # RETURN 
  return(output)
}

CalcDenomR<- cmpfun(CalcDenomR)


CalcPredictedShareR <- function(
  fromArea, 
  toArea, 
  deltas, 
  mu, 
  percentRents, 
  medval
  ){
  forceAttract <- CalcForceAttract(fromArea=fromArea, toArea=toArea, deltas=deltas, mu=mu, percentRents=percentRents, medval=medval)
  denom <- CalcDenom(areak=fromArea, deltas=deltas, mu=mu, percentRents=percentRents, medval=medval)
  
  output <- forceAttract/denom
  return(output)
}

CalcPredictedShareR <- cmpfun(CalcPredictedShareR)


CalcPredictedStayerProp <- function(
  popT2, 
  deltas, 
  mu, 
  percentRents, 
  medval,
  verbose=F
  ){
  
  # ERROR CATCHING
  
  if (any((!is.finite(popT2)))){
    stop("popT2 is not finite")
  }
  if (any(!(percentRents >= 0 & percentRents <= 1))){
    stop("some percentRents elements are not proportions")
  }
  
  
  # PROCESSING
  N <- length(deltas)
  output <- 0
  for (i in 1:N){
    tmp <- CalcPredictedShare(fromArea=i, toArea=i, deltas=deltas, mu=mu, percentRents=percentRents, medval=medval)
    tmp <- tmp * popT2[i]
    output <- output + tmp  
  }
  
  # ERROR CATCHING
  if (is.na(output)){
    stop("output is NA")
  }
  
  if (is.infinite(output)){
    stop("output is inf")
  }
  
  # OUTPUT
  output <- output/sum(popT2)
  return(output)
}


CalcPredictedStayerProp <- cmpfun(CalcPredictedStayerProp)



# INPUTS:
#   deltas :  candidate utility for each area
#   mu :      (race specific) general utility
#   oldShares : shares of residence in earlier time period
#   calcMovingCosts: function for calculating moving costs of going from area 1 to area 2

# INTERNAL FUNCTIONS
#   CalcDenom: denominator for areak 

# OUTPUTS: 
#   predictedNewShares : predicted shares in each area in later time period given deltas


PredictNewShares <- function(
  deltas, 
  mu, 
  oldCounts, 
  percentRents,
  medval,
  verbose=F
  ){
  
  
  if(verbose==T){
    print("Entered PredictNewShares")
  }

  # ERROR CHECKING
  
  if (any(oldCounts < 0)){
    stop("Some oldCounts are negative")
  }
  
  if (any(!(percentRents >= 0 & percentRents <= 1))){
    stop("some percentRents elements are not proportions")
  }
  
  # PROCESSING
  
  Nareas <- length(deltas)
  
  oldShares <- oldCounts/sum(oldCounts)

  predictedNewShares <- rep(0, Nareas)
  for (j in 1:Nareas){
    
    if (verbose==T){
      if ((j %%1)==0){
        cat("o")
      }
      if ((j %%10)==0){
        cat("|")
      }
      if ((j %%100)==0){
        cat("\n")
      }
    }
    
    for (k in 1:Nareas){
      denom <- CalcDenom(
        areak=k, 
        deltas=deltas, 
        mu=mu, 
        percentRents=percentRents, 
        medval=medval
        )
      
      forceAttract <- CalcForceAttract(
        toArea=j, 
        fromArea=k, 
        mu=mu, 
        deltas=deltas, 
        percentRents=percentRents, 
        medval=medval
        )
      
      tmp <- forceAttract / denom
      tmp <- tmp * oldShares[k]
      predictedNewShares[j] <- predictedNewShares[j] + tmp
    }
  }
  
  if (verbose){
    print("Exiting PredictNewShares")
  }
  
    return(predictedNewShares)    

}


PredictNewShares <- cmpfun(PredictNewShares)




# INPUTS:
#   initialDeltas: first guess of area utilities
#   currentMu: overall utility estimate 
#   tolerance: how far out the discrepency between any element of realNewShares and PredictedNewShares can be
#   realOldShares: proportion of people in each area in earlier time
#   realNewShares: proportion of people in each area in later time

# OUTPUT:
#   currentDeltas: New estimates for area utilities


DeltaContractionMapping <- function(
  initialDeltas,
  currentMu,
  tolerance,
  realOldCounts,
  realNewCounts,
  percentRents,
  medval,
  verbose=F
  ){
  
  if(verbose){
    print("Entered DeltaContractionMapping")
  }
  
  currentDeltas <- initialDeltas
  predictedOldShares <- Inf
  
  
  previousDeltas <- currentDeltas + tolerance * 100 # to make sure this object exists and fails first condition
  
  # This may need changing again later, because there is not a closed system
  realOldShares <- realOldCounts/sum(realOldCounts)
  realNewShares <- realNewCounts/sum(realNewCounts)
  
  inner <- 0 
  
  repeat
  {
    # BUG CATCHING
    if (any(is.na(previousDeltas))){
      stop("Some elements in previousDeltas are NA")
    }
    
    if (any(is.na(currentDeltas))){
      stop("Some elements in currentDeltas are NA")
    }
    
    tmp <- abs(previousDeltas-currentDeltas) 
    tmp <- tmp < tolerance
        
    if (all(tmp))
    {
      
      if (verbose){
        print("exiting DeltaContractionMapping loop: condition met")
      }
      
      break
    } 
    else 
    {
      
      if (verbose){
        tmp <- abs(previousDeltas - currentDeltas)
        print("repeating DeltaContractionMapping loop: condition not met")
        cat("Discrepancy: mean\t", mean(tmp), "\tmax:\t", max(tmp),"\n")
        cat("inner increased: ", inner, "\n")
      }
      
      inner <- inner + 1   
      
      
      predictedNewShares <- PredictNewShares(
        deltas=currentDeltas, 
        mu=currentMu, 
        oldCounts=realOldCounts, 
        percentRents=percentRents, 
        medval=medval,
        verbose=verbose)     
      
      updatedDeltas <- currentDeltas + log(realNewShares) - log(predictedNewShares)
      updatedDeltas <- updatedDeltas - mean(updatedDeltas)
      
      tmp <- which(abs(updatedDeltas - currentDeltas) >= tolerance)
      
      if (verbose){
        cat("Changing ", length(tmp), " elements of delta vector\n")  
      }
      
      previousDeltas <- currentDeltas
      currentDeltas[tmp] <- updatedDeltas[tmp] # change the deltas where the discrepency is greater than the tolerance
    }
  }
  
  return(currentDeltas)  
}


DeltaContractionMapping <- cmpfun(DeltaContractionMapping)


# INPUTS:
#   lowerLimit: starting lower limit for mu
#   upperLimit: starting upper limit for mu
#   muTolerance: tolerance for when to stop trying to optimise mu
#   maxIt: maximum number of iterations for trying to optimise mu


FindMu <- function(
  
  initLowerLimit=-0.15,
  initUpperLimit=0.10,
  muTol=0.001,
  deltaTol=0.001,
  maxIt=1000,

  initDeltas=NA,
  percentRents,
  medval,  
  realStayerProp,
  oldCounts,
  newCounts,
  verbose=F, 
  interim=T
  ){
  
  # ERROR CHECKING
  if (initLowerLimit > initUpperLimit){
    stop("initLowerLimit greater than initUpperLimit")
  }
  
  if (any(percentRents > 1) | any(percentRents < 0)){
    stop("percentRents not a valid proportion")
  }
  
  
  
  if (verbose) {
    print("Starting FindMu")
  }
  
  
  
  N <- length(medval)
  iteration <- 0
  
  
  upperLim <- initUpperLimit
  lowerLim <- initLowerLimit
  
  
  curMu <- (upperLim + lowerLim)/2

  if (verbose){
    cat("upperLim:\t", upperLim, "\tlowerLim:\t", lowerLim, "\tcurMu:\t", curMu, "\n")
    print("Starting first contraction mapping")
  }
  
  newDeltas <- DeltaContractionMapping(
    initialDeltas=initDeltas,
    currentMu=curMu,
    tolerance=deltaTol,
    realOldCounts=oldCounts,
    realNewCounts=newCounts,
    percentRents=percentRents,
    medval=medval,
    verbose=verbose
    )
  
  if (verbose){
    print("Finished first contraction mapping")
    print("About to enter CalcPredictedStayerProp")
  }
  
  predictedStayerProp <- CalcPredictedStayerProp(popT2=newCounts, deltas=newDeltas, mu=curMu, percentRents=percentRents, medval=medval)
  difStayerProp <- realStayerProp - predictedStayerProp
  
  if (verbose){
    cat("StayerProp. Real:\t", realStayerProp, "\tPredicted:\t",predictedStayerProp, "\tDif:\t", difStayerProp, "\n")
    print("About to enter while loop")

  }
  while( (iteration < maxIt) & (abs(difStayerProp) > muTol)){

    if (verbose){
      print("In while loop within FindMu")
      cat("FindMu iteration:\t", iteration, "\n")
    }    
    
    if (difStayerProp > 0){
      lowerLim <- curMu
      if (verbose){
        print("difStayerProp is positive")
      }
      
    } else if (difStayerProp < 0){
      upperLim <- curMu
      if (verbose){
        print("difStayerProp is negative")
      }
    } else {
      stop("difStayerProp neither positive nor negative")
    }
    
    if (verbose){
      cat("curMu is ", curMu, "\n")
    }
    
    curMu <- (lowerLim + upperLim)/2
    
    if (verbose){
      cat("curMu is now ", curMu, "\n")
    }
    
    prevDeltas <- newDeltas
    
    if (verbose){
      print("Entering DeltaContractionMapping for a subsequent time")
    }
    
    newDeltas <- DeltaContractionMapping(
      initialDeltas=prevDeltas,
      currentMu=curMu,
      tolerance=deltaTol,
      realOldCounts=oldCounts,
      realNewCounts=newCounts,
      percentRents=percentRents,
      medval=medval,
      verbose=verbose
    )
    
    if (verbose){
      print("Exited DeltaContractionMapping a subsequent time")
      print("Entering CalcPredictedStayerProp")
    }
    
    predictedStayerProp <- CalcPredictedStayerProp(
      popT2=newCounts, 
      deltas=newDeltas, 
      mu=curMu, 
      percentRents=percentRents, 
      medval=medval)
    
    difStayerProp <- realStayerProp - predictedStayerProp
    if (verbose){
      cat("difStayerProp now:\t", difStayerProp, "\n")
    }
    iteration <- iteration + 1
    
    if (verbose){
      cat("iteration now ", iteration, "\n")
      print("Reached end of while loop")
    }
    
  }
  # Now need to calculate different between predicted and actual new shares
  
  output <- list(
    mu=curMu,
    deltas=newDeltas
    )
  
  return(output)
}

# 
# 
# 
# ############################################################
# ############################################################
# 
#   while(abs(mudiff)>mutol & k1 < maxit){
#     
#     
#     # muguess = (upper+lower)/2.d0
#     muguess <<- (upper+lower)/2
#     muguesstracker<<- c(muguesstracker, muguess)
#     # CALL fcn1(xguess,muguess,x,mudiff)
#     cat("fcn1 called from outer\n")
#     fcn1(xguess, muguess, mudiff)
#     
#     # if (mudiff.gt.0.d0) lower = muguess
#     # if (mudiff.le.0.d0) upper = muguess  
#     if(mudiff > 0){
#       cat("lower changed from", lower, "to", muguess, "\n")
#       lower <<- muguess 
#     } else {
#       cat("upper changed from", upper, "to", muguess, "\n")
#       
#       #     cat("upper changed from ", lower, " to ", muguess",\n")
#       upper <<- muguess
#     }
#   
#   repeat
#   {
#     if (all((abs(updatedDeltas-currentDeltas)) < tolerance))
#     {
#       break
#     } 
#     else 
#     {
#       predictedNewShares <- PredictNewShares(deltas=currentDeltas, mu=currentMu, oldShares=realOldShares)     
#       
#       updatedDeltas <- currentDeltas + log(realNewShares) - log(predictedNewShares)
#       updatedDeltas <- updatedDeltas - mean(updatedDeltas)
#       
#       tmp <- which(abs(updatedDeltas - currentDeltas) >= tolerance)
#       currentDeltas[tmp] <- updatedDeltas[tmp] # change the deltas where the discrepency is greater than the tolerance
#     }
#   }
#   
#   return(NULL)
# }
# 
#   # Let's think about this: find mu MUST call DeltaContractionMapping
#  #      which calls predict new shares
#   #         which calls CalcMovingCosts
#   
#   
#   #What are the inputs FindMu need?
#   
#   # DATA
#   #   Popn in 2000 (T1)
#   #   Popn in 2007 (T2) 
#   #   Percent stayers
#   
#   # OTHER INPUTS
#   # Lower limit for MuSearch
#   # Upper limit for MuSearch
#   # 
#   
#   
#   currentMu <- (upperLimit + lowerLimit)/2
#   fcn1called <- fcn1called + 1
#   
#   
# 
#   
#   mu <- muguess
#   
#   totpop <- sum(pop2)
#   
#   # c
#   cat("subfcn1 called from fcn1\n")
#   subfcn1()
#   
#   #   mudiff = stayp - pstay
#   mudiff <<- stayp - pstay
#   #  cat("\n")
#   return(currentMu)
#   
# }
# 
# subfcn1 <- function(){
#   subfcn1called <<- subfcn1called + 1 # keep track of the number of times this function has been called
#   
#   # if the object 'test' does not exist, create it and give it an infinite value (So the while loop will run)
#   if (!exists("test")){ test <<- Inf} 
#   
#   # while any element exceeds the tolerance
#   while(any(test > tolc)){
#     
#     # Track the number of times the while loop has run for
#     
#     trk <<- trk+1
#     # every 10th run, print |; every other run, print .
#     if (trk%%10==0) cat("|") else cat(".")
#     # every 100th run, print the average of the discrepencies, then start a new line
#     if (trk%%100==0) cat(" ", mean(test),"\n")
#     
#     #create an empty matrix to hold utilty estimates of going from area k to area j
#     util <<- matrix(NA, nrow=n1, ncol=n1)
#     for (j in 1:n1){
#       for (k in 1:n1){
#         mc <- GetMovingCost(from.area=k, to.area=j, N=n1, percent.rents=rperc, medval=medval)      
#         util[k, j] <<- x1[k] - x1[j] - mu*mc        
#       }
#     }
# 
# 
#     
#     # for each column in exputil, calculate the sum, and return as an n+1 length vector
#     denom <- apply(exputil,2, sum)
# 
#     rhs <- rep(0, n1)
#     for (i in 1:n1){
#       for (j in 1:n1){
#         # the numerator is the exp(util) of going from area j to area i
#         numer <- exputil[i,j]
#         # rhs for area i is sum of exp(utils) of all areas j to i, divided by sum of exp utils for area j, multiplied
#         # by popn of area js at time t=1
#         rhs[i] <- rhs[i] + pop1[j]*(numer/denom[j]) 
#       }
#     }
#     
#     # stayer percentages are sums of  diagonals of exputilities (U11, U22, U33 etc), multiplied populations of each area j in time 1, divided by denom of area j
#     stay <<- sum(pop1* diag(exputil)/denom)
#     #... then divided by totpop-change in population from 2001 to 2009
#     stay <<- stay/(totpop-pop2[n1])
# 
#     # predicted stayer %: sum of the diagonals of exputil matrix, multiplied by the number of areal units, divided by vector of area denominators
#     pstay <<- sum(n * diag(exputil)/denom)
#     
#     share <<- pop2/totpop
#     pshare <<- rhs/totpop
#     
#     
#     xf <<- x1 + log(share) - log(pshare)
#     xavg <<- mean(xf)
#     xf <<- xf - xavg
#     
#     # enddo
#     # c     
#     
#     test <<- abs(xf - x1)
#     
#     x1 <<- xf
#   }
# }
# #SUBROUTINE fcn1(x1,muguess,xf,mudiff)
# fcn1 <- function(x1, muguess, mudiff){
#   fcn1called <<- fcn1called + 1
#   
#   
#   #c
#   #ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   #c
#   #USE numerical_libraries
#   #c
#   #PARAMETER (tolc=0.000001d0)
#   #PARAMETER (n=1921,n1=n+1)
#   #INTEGER i,j
#   # DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
#   # &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
#   # &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
#   # &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
#   # &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
#   # &   rperc(n),mc
#   # c
#   # common /data1/ pop1,pop2,medval,rperc
#   # common /data2/ stayp,pstay
#   # common /distance/ dist
#   # c
#   # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   # c
#   # c     
#   # mu = muguess
#   mu <<- muguess
# 
#   
# 
#   
#   totpop <<- sum(pop2)
#   
#   # c
#   cat("subfcn1 called from fcn1\n")
#   subfcn1()
#   
#   #   mudiff = stayp - pstay
#   mudiff <<- stayp - pstay
# #  print(mudiff)
#   mudifftracker <<- c(mudifftracker, mudiff) # look for convergence
# #  cat("\n")
#   return(NULL)
#   
#   # END  
# }
# 
# # Check all values are within tolerance
# # if all values are within tolerance, end function, and 
# # return mudiff
# 
# # c
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # 
# # }
# 
# 
# ################################################
# #ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# #c
# subfcn1.old <- function(){
#   
#   trk <<- trk+1
#   if (trk%%10==0) cat("|") else cat(".")
#   if (trk%%100==0) cat("\n")
#   
#   # 300  continue
#   # do j = 1,n1
#   # do k = 1,n1
#   for (j in 1:n1){
#     for (k in 1:n1){
#       
#       # mc = 0.d0
#       mc <- 0
#       
#       #     if ((j.ne.k).and.(j.ne.n1).and.(k.ne.n1)) 
#       #       mc = (2910.d0+(0.03d0*(1.d0-rperc(j))*medval(j))+(0.03d0*(1.d0-rperc(k))*medval(k)))/24.556d0
#       if ((j != k) & (j != n1) & (k!=n1)){
#         mc <- (2910 + (0.03 * (1 - rperc[j]) * medval[j]) + (0.03 * (1 - rperc[k])*medval[k]))/24.556
#       }
#       
#       #     if ((j.ne.k).and.((j.eq.n1).or.(k.eq.n1))) 
#       #       mc = 528.71d0
#       if ((j!=k) & ((j ==n1) | (k ==n1))){      
#         mc <- 528.71
#       }
#       
#       #    util(k,j) = x1(k)-x1(j)-mu*mc # utility of going from area k to area j    
#       util[k, j] <<- x1[k] - x1[j] - mu*mc
#       
#       # enddo
#       # enddo    
#     }
#   }
#   
#   # 
#   #do j = 1,n1
#   for (j in 1:n1){
#     #   denom(j) = 0.d0
#     denom[j] <- 0  
#     
#     #   do k = 1,n1
#     for (k in 1:n1){
#       # denom(j) = denom(j)+dexp(util(k,j)) # exponential of util
#       denom[j] <- denom[j] + exp(util[k,j])
#       # enddo
#       # enddo      
#     }
#   }
#   
#   # denominator for area j is sum of exponentials for area utilities 
#   # for all other js  
#   
#   # do i = 1,n1
#   
#   for (i in 1:n1){
#     #   rhs(i) = 0.d0
#     rhs[i] <- 0
#     #   do j = 1,n1
#     for (j in 1:n1){
#       # numer = dexp(util(i,j))
#       numer <- exp(util[i,j])
#       # rhs(i) = rhs(i)+(numer/denom(j))*pop1(j)
#       rhs[i] <- rhs[i] + (numer/denom[j]) * pop1[j]
#       #     enddo
#       #     enddo
#       
#     }
#   }
#   
#   # 
#   # stay = 0.d0
#   stay <- 0
#   
#   # do i = 1,n
#   for (i in 1:n){
#     #   stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
#     stay <- stay + (exp(util[i,j]))/denom[i]*pop1[i]  
#     #   enddo  
#   }
#   
#   # stay = stay/(totpop-pop2(n1))
#   stay <- stay/(totpop-pop2[n1])
#   # c
#   # pstay = 0.d0
#   pstay <<- 0
#   # do i = 1,n
#   for (i in 1:n){
#     #   pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
#     pstay <<- pstay + (exp(util[i,j])/denom[i])/(1 * n)  
#     #   enddo
#     
#   }
#   # c
#   # do i = 1,n1
#   share <<- rep(NA, n1)
#   pshare <<- rep(NA, n1)
#   for (i in 1:n1){
#     #   share(i) = pop2(i)/totpop
#     share[i] <<- pop2[i]/totpop
#     #   pshare(i) = rhs(i)/totpop
#     pshare[i] <<- rhs[i]/totpop
#   }
#   
#   # enddo
#   # c
#   # xavg = 0.d0
#   xavg <<- 0
#   xf <<- rep(NA, n1)
#   # DO j = 1,n1
#   for (j in 1:n1){
#     # Contraction mapping : adding difference in log shares
#     # xf: updated value of x
#     # x1 : previous value of x
#     # normalised by subtracting mean from x
#     
#     # xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
#     xf[j] <<- x1[j] + log(share[j]) - log(pshare[j])
#     
#     #   c         write(*,*) j,share(j),pshare(j)
#     
#     #   xavg = xavg+(xf(j)/(1.d0*n1))
#     xavg <<- xavg + xf[j]/(1*n1)
#     
#     #   ENDDO  
#   }
#   
#   
#   
#   # 225  format(i10,f15.6,f15.6)
#   # c
#   # temp = xavg 
#   # normalisation (needs to be done)
#   temp <<- xavg
#   
#   # DO j = 1,n1
#   for (j in 1:n1){
#     #   xf(j) = xf(j)-temp # normalised updated x estimate
#     xf[j] <<- xf[j] - temp
#     
#   }
#   # enddo
#   # c     
#   test <<- rep(NA, n1)
#   # DO j = 1,n1
#   for (j in 1:n1){
#     #   test(j) = dabs(xf(j)-x1(j))
#     test[j] <<- abs(xf[j] - x1[j])
#     #   IF (test(j).gt.tolc) then
#     if (test[j] > tolc){
#       #     DO k = 1,n1
#       for (k in 1:n1){
#         # x1(k) = xf(k)
#         x1[k] <<- xf[k]
#         # ENDDO
#       }
#       # go back to subfunction
#       #     GO TO 300
#       
#       subfcn1.old()
#       # ENDIF
#     }
#     # ENDDO
#   }
# }
# #SUBROUTINE fcn1(x1,muguess,xf,mudiff)
# fcn1.old <- function(x1, muguess, xf, mudiff){
#   
#   
#   
#   #c
#   #ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   #c
#   #USE numerical_libraries
#   #c
#   #PARAMETER (tolc=0.000001d0)
#   #PARAMETER (n=1921,n1=n+1)
#   #INTEGER i,j
#   # DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
#   # &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
#   # &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
#   # &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
#   # &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
#   # &   rperc(n),mc
#   # c
#   # common /data1/ pop1,pop2,medval,rperc
#   # common /data2/ stayp,pstay
#   # common /distance/ dist
#   # c
#   # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   # c
#   # c     
#   # mu = muguess
#   mu <<- muguess
#   x1 <<- x1
#   
#   denom <<- rep(NA, n1)
#   rhs <<- rep(NA, n1)
#   
#   util <<- matrix(NA, nrow=n1, ncol=n1)
#   
#   # totpop = 0.d0
#   totpop <<- 0
#   for (i in 1:n1){
#     
#     # DO i = 1,n1
#     # totpop = totpop+pop2(i)
#     totpop <<- totpop + pop2[i]
#     # Total in the second period
#     # Shares are shares based on sum of pop2
#     # ENDDO
#   }
#   # c
#   subfcn1.old()
#   
#   #   mudiff = stayp - pstay
#   mudiff <- stayp - pstay
#   print(mudiff)
#   mudifftracker <<- c(mudifftracker, mudiff) # look for convergence
#   cat("\n")
#   return(mudiff)
#   
#   # END  
# }
# 
# # Check all values are within tolerance
# # if all values are within tolerance, end function, and 
# # return mudiff
# 
# # c
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # 
# # }
# 
# # For producing simulated area counts given the available count data, but also incorporating
# # a noninformative prior of 0.5 to all cell counts
# MakeDirichletFreqs <- function(InputData, M, prior=0.5, as.freqs=F){
#   require(MCMCpack)
#   N <- dim(InputData)[1]  
#   tmp <- vector("list", N)
#   if (as.freqs==T){
#     area.totals <- apply(InputData,1, sum) + prior
#   }
#   
#   # by default, the input to rdirichlet is a vector of the frequencies in each area
#   # the output from rdirichlet is a dataframe where each row is a different estimate 
#   # for that area. 
#   # This will require rearranging later
#   for (i in 1:N){
#     tmp[[i]] <- rdirichlet(M, t(InputData[i,] + prior) )
#   }
#   
#   # This code rearranges the output tmp from above
#   # so that the final output from the function is
#   # M sets of estimated frequencies, where each row 
#   # indicates a different area
#   
#   Output <- vector("list", M)
#   for (i in 1:M){
#     Output[[i]] <- t(sapply(tmp, function(x) x[i,] ))
#     if (as.freqs==T){
#       # n.b. double check this works as intended 
#       Output[[i]] <- Output[[i]] * area.totals
#     }
#   }
#   
#   return(Output)
# }
# 
# 
# # 
# # 
# # GuessMu <- function(
# #   #  lower = -0.15d0
# #   #  upper = 0.1d0
# #   lower.start= -0.15,
# #   upper.start=  0.10,
# #   maxit=10^3,
# #   mutol=10^-8,
# #   n.areas
# #   
# # ){
# #   
# #   # c
# #   #       do j = 1,n1
# #   #          xguess(j) = 0.d0
# #   #       enddo  
# #   xguess <- rep(0, n.areas)
# #   upper <- upper.start
# #   lower <- lower.start
# #   for (i in 1:maxit){
# #     muguess <- (upper + lower) / 2
# #     mudiff <- CalcImplications(xguess, muguess, x, mudiff)
# #     
# #     
# #   }
# #   
# # }
# # # c
# # # c
# # #       do k1 = 1,maxit
# # #          muguess = (upper+lower)/2.d0
# # #          CALL fcn1(xguess,muguess,x,mudiff)
# # #          if (mudiff.gt.0.d0) lower = muguess
# # #          if (mudiff.le.0.d0) upper = muguess
# # #          do j = 1,n1
# # #             xguess(j) = x(j)
# # #          enddo
# # #          write(*,40) k1,muguess,mudiff,lower,upper
# # #  40      format(i6,5f9.5)
# # #          if(dabs(mudiff).le.mutol) go to 120
# # #       enddo
# # #  120  continue
# # #       do j = 1,n
# # #          delta(j) = x(j)
# # #       enddo
# # # c
# # #       do j = 1,n
# # #          write(30,*) delta(j)
# # # c         write(*,*) j,delta(j)
# # #       enddo
# # # c
# # #       write(*,*) 'delta(n1) =',x(n1)
# # #       write(*,*) 'pstay     =',pstay
# # #       END
# # # c
