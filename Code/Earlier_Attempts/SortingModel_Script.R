rm(list=ls())
setwd("X:/Manuscripts/01 Early Draft/Sorting model with UK Data/R model/")

require(Rcpp)
sourceCpp("SortingModel_Functions.cpp")

source("SortingModel_Functions.R")

# 
# # Start with a 'toy' example:
# 
# countsT1 <- c(4, 6, 2)
# countsT2 <- c(6, 3, 3)
# 
# stayerProp <- 0.6
# 
# medval <- rep(1000, 3)
# propRent <- c(0.5, 0.3, 0.7)
# 
# 
# 
# # EVERYTHING apart from calibration parameters is derived from the above
# 
# sharesT1 <- countsT1/sum(countsT1)
# sharesT2 <- countsT2/sum(countsT2)
# 
# 

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

require(compiler)

MakeMovingCostMatrix <- cmpfun(MakeMovingCostMatrix)



###################################################################################################################################
# Input Data
###################################################################################################################################
Data0 <- read.csv("inputs/InputsCombined.csv")
Data1 <- read.csv("inputs/stage2_dec_vc_sheet1.csv")
Data2 <- read.csv("inputs/stage2_dec_vc_sheet2.csv")

require(reshape)
Data0 <- rename(Data0, c(area="id"))


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



proprents = Data012$prent
medval = Data012$medval




#Data012 <- Data012[1:300,] # for testing algorithm
movingcostmatrix <- MakeMovingCostMatrix(proprents,medval)

In <- list(
  counts = list(
    oldcounts=countsold,
    newcounts=countsnew
  ),
  movingcostmatrix = movingcostmatrix,
  utils = list(
    mu = 0.5,
    deltas = rep(0, length(countsold))
  ),
  dbg = list(
    dbgflag= T,
    verboseflag=T
    )  
)




##################################################################################################################################
###################################################################################################################################
# Global variables
###################################################################################################################################


##############################################################################################

# cccccccccccccccccccccccccccccccc
# c Stay-put % by race:          c
# c    - white_l = 0.4174        c
# c    - black_l = 0.3254        c
# c    - asian_l = 0.2808        c
# c    - hisp_l = 0.2892          c
# c    - other_l = 0.2687        c
# c                              c
# c    - white_m = 0.4082        c
# c    - black_m = 0.3933        c
# c    - asian_m =         c
# c    - hisp_m =          c
# c    - other_m =         c
# c                              c
# c    - white_h =         c
# c    - black_h =         c
# c    - asian_h =         c
# c    - hisp_h =          c
# c    - other_h =         c
# cccccccccccccccccccccccccccccccc
# stayp = 0.4174d0

# Stayer % 
stayp <- 0.4174

# mutol = 0.00000001d0

# 1t: t1 (first time period)
# 2t: t2 (second time period)
  n <- dim(Data012)[1]
  n1 <- n+1
# with 'patch' added
#  pop1t <- rep(NA, n)
#  pop2t <- rep(NA, n)
#  patch <- 0.001
# 'prior' as patch
patch <- 0.5
#  white_l00 <- Data012$white.l00
#  white_l09 <- Data012$white.l09


# do i = 1,n
pop1t <- Data012$white.l00[1:20] + patch
pop2t <- Data012$white.l09[1:20] + patch
rperc <- Data012$prent[1:20]
medval <- Data012$medval[1:20]


output <- StageOne(
  pop1t, 
  pop2t,
  rperc,
  medval,
  stayp,
  T)


######################################
system.time(
  output <- FindMu(
    initLowerLimit=-0.15, 
    initUpperLimit=0.1, 
    percentRents=rperc, 
    initDeltas=rep(0, length(pop1t)),
    medval=medval, 
    realStayerProp=stayp, 
    oldCounts=pop1t, 
    newCounts=pop2t,
    verbose=T)
)





  initLowerLimit=-0.15, 
  initUpperLimit=0.1, 
  percentRents=, 
  initDeltas=rep(0, length(pop1t)),
  medval=, 
  realStayerProp=, 
  oldCounts=pop1t, 
  newCounts=pop2t,
  verbose=T)
# USING 20 observations


# USING ALL R Functions
# VERBOSE =T
# user  system elapsed 
# 41.26    0.20   41.17 


# AS ABOVE, BUT CONVERTING GetMovingCost to C
# VERBOSE = T
# user  system elapsed 
# 12.68    0.04   12.71 


# Now with GetMovingCost, CalcDenom, and CalcForceAttract in  C
# VERBOSE = T
# user  system elapsed 
# 1.22    0.0.   1.20




# pop00 <- sum(pop1t)
# pop09 <- sum(pop2t)
# 
# dpop <- pop09-pop00
# 
# pop1 <- pop1t
# pop2 <- pop2t
# 
# 
# AreaData <- data.frame(medval=Data012$medval, percentRents=...
# 
# # pop1(n1) = 2.d0*dabs(dpop)
# pop1 <- c(pop1, 
#           2 * abs(dpop)
#         )
# 
# pop2 <- c(
#   pop2,
#   pop1[n1]-dpop
#   )
# 
# 
# xguess <- rep(0, n1)
# 
# x1 <- rep(0, n1)
# # c
# # c
# # lower = -0.15d0
# lower <- -60
# # upper = 0.1d0
# upper <- 40
# 
# maxit <- 1000
# 
# #mutol <- 0.00000001
# #tolc <- 0.000001
# 
# mutol <- 0.001
# tolc <- 0.001
# 
# mudiff <-  Inf # to allow while loop to work

# mudifftracker <- c()  # use for checking convergence in later functions

#debug(fcn1)
#debug(subfcn1)

# trk <<- 0 # Tracker variable
# k1 <<- 1
# # do k1 = 1,maxit
# 
# fcn1called <<- 0
# subfcn1called <<- 0
# muguesstracker <<- c()
# 
# debug(subfcn1)
# while(abs(mudiff)>mutol & k1 < maxit){
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
#       } else {
#       cat("upper changed from", upper, "to", muguess, "\n")
# 
#  #     cat("upper changed from ", lower, " to ", muguess",\n")
#       upper <<- muguess
#   }
#     
# #    cat("mudiff: ", mudiff, "\n")
#  
#     
#     # do j = 1,n1
# #     for (j in 1:(n1)){
# #       # xguess(j) = x(j)
# #       xguess[j] <<- x1[j]
# #       # enddo  
# #     }
#     xguess <<- x1
#     
#     # write(*,40) k1,muguess,mudiff,lower,upper
#     # 40      format(i6,5f9.5)
#     
#     # if(dabs(mudiff).le.mutol) go to 120
#     # ###### REPLACED WITH A WHILE LOOP
#     # enddo 
#     k1 <<- k1 + 1
#     browser()
# }
# 
# # 120  continue
# 
# 
# delta <- x
# 
# 
# # NONE OF THE FOLLOWING NEEDED
# 
# # # c
# # # do j = 1,n
# # for (j in 1:n){
# #   
# #   
# # }
# # # write(30,*) delta(j)
# # # c         write(*,*) j,delta(j)
# # # enddo
# # # c
# # # write(*,*) 'delta(n1) =',x(n1)
# # # write(*,*) 'pstay     =',pstay
# # # END
# # # c
# # 
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############################################################################################
# 
# # N.reps # number of replications (Dirichlet method)
# # # Proportion of stayers in 2000, by group
# # white.l <- 0.4174
# # black.l <- 0.3254
# # asian.l <- 0.2808
# # hisp.l <- 0.2892
# # other.l <- 0.2687
# # 
# # # Proportion of stayers in 2009, by group
# # white.m <- 0.4082
# # black.m <- 0.3933
# # 
# # # stayp <- 0.4174
# # mutol <- 10^-8
# # maxit <- 10^3 # the maximum number of tries for convergence
# # 
# # count.correction <- 10^-3
# # 
# # moving.cost.init <- 0
# # 
# # moving.cost.out <- 528.71 # cost of moving out of LA
# 
# 
# 
# 
# # cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# #       PROGRAM estimate
# # c
# # cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# #       USE numerical_libraries
# # c   
# #       PARAMETER (n=1921,n1=n+1,maxit=1000,nn=(n-1)*(n-1))
# 
# # Derived variables
# 
# N <- n + 1
# nn <- (n-1)*(n-1)
# 
# #       INTEGER i,j,tract(n),id1(nn),id2(nn),id3(n)
# 
# # nn implies the objects below may be better handled as matrices rather than vectors
# id1 <- matrix(NA, nrow=n, ncol=n)
# id2 <- id1
# id3 <- id2
# 
# # Let's try to understand the data first
# 
# # l: low income
# # m: medium income
# # h: high income
# 
# attach(Data)
# 
# white.00 <- white.l00 + white.m00 + white.h00
# white.09 <- white.l09 + white.m09 + white.h09
# 
# black.00 <- black.l00 + black.m00 + black.h00
# black.09 <- black.l09 + black.m09 + black.h09
# 
# asian.00 <- asian.l00 + asian.m00 + asian.h00
# asian.09 <- asian.l09 + asian.m09 + asian.h09
# 
# hisp.00 <- hisp.l00 + hisp.m00 + hisp.h00
# hisp.09 <- hisp.l09 + hisp.m09 + hisp.h09
# 
# other.00 <- other.l00 + other.m00 + other.h00
# other.09 <- other.l09 + other.m09 + other.h09
# 
# total.00 <- white.00 + black.00 + asian.00 + hisp.00 + other.00
# total.09 <- white.09 + black.09 + asian.09 + hisp.09 + other.09
# 
# # totals
# total.l00 <- white.l00 + black.l00 + asian.l00 + hisp.l00 + other.l00
# total.m00 <- white.m00 + black.m00 + asian.m00 + hisp.m00 + other.m00
# total.h00 <- white.h00 + black.h00 + asian.h00 + hisp.h00 + other.h00
# 
# total.l09 <- white.l09 + black.l09 + asian.l09 + hisp.l09 + other.l09
# total.m09 <- white.m09 + black.m09 + asian.m09 + hisp.m09 + other.m09
# total.h09 <- white.h09 + black.h09 + asian.h09 + hisp.h09 + other.h09
# 
# # percentages
# white.l00.pc <- white.l00 / total.l00
# white.m00.pc <- white.m00 / total.m00
# white.h00.pc <- white.h00 / total.h00
# 
# white.l09.pc <- white.l09 / total.l09
# white.m09.pc <- white.m09 / total.m09
# white.h09.pc <- white.h09 / total.h09
# 
# black.l00.pc <- black.l00 / total.l00
# black.m00.pc <- black.m00 / total.m00
# black.h00.pc <- black.h00 / total.h00
# 
# black.l09.pc <- black.l09 / total.l09
# black.m09.pc <- black.m09 / total.m09
# black.h09.pc <- black.h09 / total.h09
# 
# asian.l00.pc <- asian.l00 / total.l00
# asian.m00.pc <- asian.m00 / total.m00
# asian.h00.pc <- asian.h00 / total.h00
# 
# asian.l09.pc <- asian.l09 / total.l09
# asian.m09.pc <- asian.m09 / total.m09
# asian.h09.pc <- asian.h09 / total.h09
# 
# hisp.l00.pc <- hisp.l00 / total.l00
# hisp.m00.pc <- hisp.m00 / total.m00
# hisp.h00.pc <- hisp.h00 / total.h00
# 
# hisp.l09.pc <- hisp.l09 / total.l09
# hisp.m09.pc <- hisp.m09 / total.m09
# hisp.h09.pc <- hisp.h09 / total.h09
# 
# other.l00.pc <- other.l00 / total.l00
# other.m00.pc <- other.m00 / total.m00
# other.h00.pc <- other.h00 / total.h00
# 
# other.l09.pc <- other.l09 / total.l09
# other.m09.pc <- other.m09 / total.m09
# other.h09.pc <- other.h09 / total.h09
# 
# # Total population sizes by group
# 
# white.totpop.l00 <- sum(white.l00)
# black.totpop.l00 <- sum(black.l00)
# asian.totpop.l00 <- sum(asian.l00)
# hisp.totpop.l00 <- sum(hisp.l00)
# other.totpop.l00 <- sum(other.l00)
# 
# white.totpop.l09 <- sum(white.l09)
# black.totpop.l09 <- sum(black.l09)
# asian.totpop.l09 <- sum(asian.l09)
# hisp.totpop.l09 <- sum(hisp.l09)
# other.totpop.l09 <- sum(other.l09)
# 
# white.totpop.m00 <- sum(white.m00)
# black.totpop.m00 <- sum(black.m00)
# asian.totpop.m00 <- sum(asian.m00)
# hisp.totpop.m00 <- sum(hisp.m00)
# other.totpop.m00 <- sum(other.m00)
# 
# white.totpop.m09 <- sum(white.m09)
# black.totpop.m09 <- sum(black.m09)
# asian.totpop.m09 <- sum(asian.m09)
# hisp.totpop.m09 <- sum(hisp.m09)
# other.totpop.m09 <- sum(other.m09)
# 
# white.totpop.h00 <- sum(white.h00)
# black.totpop.h00 <- sum(black.h00)
# asian.totpop.h00 <- sum(asian.h00)
# hisp.totpop.h00 <- sum(hisp.h00)
# other.totpop.h00 <- sum(other.h00)
# 
# white.totpop.h09 <- sum(white.h09)
# black.totpop.h09 <- sum(black.h09)
# asian.totpop.h09 <- sum(asian.h09)
# hisp.totpop.h09 <- sum(hisp.h09)
# other.totpop.h09 <- sum(other.h09)
# 
# # Differences in populations over time
# 
# white.l.dif <- white.l09 - white.l00
# black.l.dif <- black.l09 - black.l00
# asian.l.dif <- asian.l09 - asian.l00
# hisp.l.dif  <- hisp.l09 - hisp.l00
# other.l.dif <- other.l09 - other.l00
# 
# white.m.dif <- white.m09 - white.m00
# black.m.dif <- black.m09 - black.m00
# asian.m.dif <- asian.m09 - asian.m00
# hisp.m.dif  <- hisp.m09 - hisp.m00
# other.m.dif <- other.m09 - other.m00
# 
# white.h.dif <- white.h09 - white.h00
# black.h.dif <- black.h09 - black.h00
# asian.h.dif <- asian.h09 - asian.h00
# hisp.h.dif  <- hisp.h09 - hisp.h00
# other.h.dif <- other.h09 - other.h00
# 
# # Differences in populations over time as a % of 2000 size
# 
# white.l.dif.pc <- white.l.dif / white.l00
# black.l.dif.pc <- black.l.dif / black.l00
# asian.l.dif.pc <- asian.l.dif / asian.l00
# hisp.l.dif.pc <- hisp.l.dif / hisp.l00
# other.l.dif.pc <- other.l.dif / other.l00
# 
# white.m.dif.pc <- white.m.dif / white.m00
# black.m.dif.pc <- black.m.dif / black.m00
# asian.m.dif.pc <- asian.m.dif / asian.m00
# hisp.m.dif.pc <- hisp.m.dif / hisp.m00
# other.m.dif.pc <- other.m.dif / other.m00
# 
# white.h.dif.pc <- white.h.dif / white.h00
# black.h.dif.pc <- black.h.dif / black.h00
# asian.h.dif.pc <- asian.h.dif / asian.h00
# hisp.h.dif.pc <- hisp.h.dif / hisp.h00
# other.h.dif.pc <- other.h.dif / other.h00
# 
# 
# # have following equations:
# 
# #condition to maximise likelihood of:
#  
# # Choose area utilities
# # utility multiplier of moving cost
# # such that the discrepency between
# #   sigma.2007.real
# # and 
# #   sigma.2007.predicted
# # is minimised
# 
# # where
# 
# # Functions: 
# 
# 
# # area index: which area pop share is being predicted for
# # util.area: current estimate of area utilities (vector)
# # util.ind: current estimate of individual utility assoc with moving (scalar)
# # popshares.earlier: earlier 
# 
# 
# CompareRealPredicted()
# 
# 
# 
# 
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# #       SUBROUTINE fcn1(x1,muguess,xf,mudiff)
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# #       USE numerical_libraries
# # c
# #       PARAMETER (tolc=0.000001d0)
# #       PARAMETER (n=1921,n1=n+1)
# #       INTEGER i,j
# #       DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
# #      &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
# #      &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
# #      &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
# #      &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
# #      &   rperc(n),mc
# # c
# #       common /data1/ pop1,pop2,medval,rperc
# #       common /data2/ stayp,pstay
# #       common /distance/ dist
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# # c     
# fcn1 <- function(x1, muguess, xf, mudiff){
#   mu <- muguess  
#   totpop <- sum(pop2)
# 
#   
# 
#   
# # #   c
# # #   300  continue
# # 
# # 
# # #   
# #   
# #   if ((j.ne.k).and.((j.eq.n1).or.(k.eq.n1))) mc = 528.71d0
#   
# #   util(k,j) = x1(k)-x1(j)-mu*mc
# #   enddo
# #   enddo
# #   do j = 1,n1
# #   denom(j) = 0.d0
# #   do k = 1,n1
# #   denom(j) = denom(j)+dexp(util(k,j))
# #   enddo
# #   enddo
# #   do i = 1,n1
# #   rhs(i) = 0.d0
# #   do j = 1,n1
# #   numer = dexp(util(i,j))
# #   rhs(i) = rhs(i)+(numer/denom(j))*
# #     &               pop1(j)
# #   enddo
# #   enddo
# #   stay = 0.d0
# #   do i = 1,n
# #   stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
# #   enddo
# #   stay = stay/(totpop-pop2(n1))
# #   c
# #   pstay = 0.d0
# #   do i = 1,n
# #   pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
# #   enddo
# #   c
# #   do i = 1,n1
# #   share(i) = pop2(i)/totpop
# #   pshare(i) = rhs(i)/totpop
# #   enddo
# #   c
# #   xavg = 0.d0
# #   DO j = 1,n1
# #   xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
# #   c         write(*,*) j,share(j),pshare(j)
# #   xavg = xavg+(xf(j)/(1.d0*n1))
# #   ENDDO
# #   225  format(i10,f15.6,f15.6)
# #   c
# #   temp = xavg
# #   DO j = 1,n1
# #   xf(j) = xf(j)-temp
# #   enddo
# #   c      
# #   DO j = 1,n1
# #   test(j) = dabs(xf(j)-x1(j))
# #   IF (test(j).gt.tolc) then
# #   DO k = 1,n1
# #   x1(k) = xf(k)
# #   ENDDO
# #   GO TO 300
# #   ENDIF
# #   ENDDO
# #   mudiff = stayp - pstay
# #   c
# #   END
# #   c
# #   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# #   
# #   ##########################
# }
# 
# # cccccccccccccccccccccccccccccccc
# # c Stay-put % by race:          c
# # c    - white_l = 0.4174        c
# # c    - black_l = 0.3254        c
# # c    - asian_l = 0.2808        c
# # c    - hisp_l = 0.2892          c
# # c    - other_l = 0.2687        c
# # c                              c
# # c    - white_m = 0.4082        c
# # c    - black_m = 0.3933        c
# # c    - asian_m =         c
# # c    - hisp_m =          c
# # c    - other_m =         c
# # c                              c
# # c    - white_h =         c
# # c    - black_h =         c
# # c    - asian_h =         c
# # c    - hisp_h =          c
# # c    - other_h =         c
# # cccccccccccccccccccccccccccccccc
# # stayp = 0.4174d0
# # mutol = 0.00000001d0
# 
# # 1t: t1 (first time period)
# # 2t: t2 (second time period)
# 
# # with 'patch' added
# 
# # do i = 1,n
# # pop1t(i) = white_l00(i)+0.001d0
# # pop2t(i) = white_l09(i)+0.001d0
# # pop00 = pop00+pop1t(i)
# # pop09 = pop09+pop2t(i)
# # enddo
# # 555  format(i9,3f21.7)
# # dpop = pop09-pop00
# # write(*,*) 'dpop =',dpop
# # do i = 1,n
# # pop1(i) = pop1t(i)
# # pop2(i) = pop2t(i)
# # enddo
# # pop1(n1) = 2.d0*dabs(dpop)
# # pop2(n1) = pop1(n1)-dpop
# # c
# # do j = 1,n1
# # xguess(j) = 0.d0
# # enddo
# # c
# # c
# # lower = -0.15d0
# # upper = 0.1d0
# # do k1 = 1,maxit
# # muguess = (upper+lower)/2.d0
# # CALL fcn1(xguess,muguess,x,mudiff)
# # if (mudiff.gt.0.d0) lower = muguess
# # if (mudiff.le.0.d0) upper = muguess
# # do j = 1,n1
# # xguess(j) = x(j)
# # enddo
# # write(*,40) k1,muguess,mudiff,lower,upper
# # 40      format(i6,5f9.5)
# # if(dabs(mudiff).le.mutol) go to 120
# # enddo
# # 120  continue
# # do j = 1,n
# # delta(j) = x(j)
# # enddo
# # c
# # do j = 1,n
# # write(30,*) delta(j)
# # c         write(*,*) j,delta(j)
# # enddo
# # c
# # write(*,*) 'delta(n1) =',x(n1)
# # write(*,*) 'pstay     =',pstay
# # END
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# # SUBROUTINE fcn1(x1,muguess,xf,mudiff)
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# # USE numerical_libraries
# # c
# # PARAMETER (tolc=0.000001d0)
# # PARAMETER (n=1921,n1=n+1)
# # INTEGER i,j
# # DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
# # &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
# # &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
# # &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
# # &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
# # &   rperc(n),mc
# # c
# # common /data1/ pop1,pop2,medval,rperc
# # common /data2/ stayp,pstay
# # common /distance/ dist
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# # c
# # c     
# # mu = muguess
# # totpop = 0.d0
# # DO i = 1,n1
# # totpop = totpop+pop2(i)
# # ENDDO
# # c
# # 300  continue
# # do j = 1,n1
# # do k = 1,n1
# # mc = 0.d0
# # if ((j.ne.k).and.(j.ne.n1).and.
# #     &         (k.ne.n1)) mc = (2910.d0+(0.03d0*(1.d0-
# #                                                   &         rperc(j))*medval(j))+(0.03d0*(1.d0-
# #                                                                                             &         rperc(k))*medval(k)))/24.556d0
# # if ((j.ne.k).and.((j.eq.n1).or.
# #                   &         (k.eq.n1))) mc = 528.71d0
# # util(k,j) = x1(k)-x1(j)-mu*mc
# # enddo
# # enddo
# # do j = 1,n1
# # denom(j) = 0.d0
# # do k = 1,n1
# # denom(j) = denom(j)+dexp(util(k,j))
# # enddo
# # enddo
# # do i = 1,n1
# # rhs(i) = 0.d0
# # do j = 1,n1
# # numer = dexp(util(i,j))
# # rhs(i) = rhs(i)+(numer/denom(j))*
# #   &               pop1(j)
# # enddo
# # enddo
# # stay = 0.d0
# # do i = 1,n
# # stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
# # enddo
# # stay = stay/(totpop-pop2(n1))
# # c
# # pstay = 0.d0
# # do i = 1,n
# # pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
# # enddo
# # c
# # do i = 1,n1
# # share(i) = pop2(i)/totpop
# # pshare(i) = rhs(i)/totpop
# # enddo
# # c
# # xavg = 0.d0
# # DO j = 1,n1
# # xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
# # c         write(*,*) j,share(j),pshare(j)
# # xavg = xavg+(xf(j)/(1.d0*n1))
# # ENDDO
# # 225  format(i10,f15.6,f15.6)
# # c
# # temp = xavg
# # DO j = 1,n1
# # xf(j) = xf(j)-temp
# # enddo
# # c      
# # DO j = 1,n1
# # test(j) = dabs(xf(j)-x1(j))
# # IF (test(j).gt.tolc) then
# # DO k = 1,n1
# # x1(k) = xf(k)
# # ENDDO
# # GO TO 300
# # ENDIF
# # ENDDO
# # mudiff = stayp - pstay
# # c
# # END
# # c
# # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
