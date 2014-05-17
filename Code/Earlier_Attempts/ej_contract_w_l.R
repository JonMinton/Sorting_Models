rm(list=ls())
setwd("X:/Manuscripts/01 Early Draft/Sorting model with UK Data/R model/")

# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
#       PROGRAM estimate
# c
# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
#       USE numerical_libraries
# c   
#       PARAMETER (n=1921,n1=n+1,maxit=1000,nn=(n-1)*(n-1))

n <- 1921
n1 <- n + 1
maxit <- 1000
nn <- (n-1)*(n-1)

#       INTEGER i,j,tract(n),id1(nn),id2(nn),id3(n)
i <- c()
j <- c()
tract <- vector("numeric", n)
id1 <- vector("numeric", nn)
id2 <- vector("numeric", nn)
id3 <- vector("numeric", nn)

#       DOUBLE PRECISION pop1(n1),pop2(n1),pop1t(n),
#      &   pop2t(n),pop00,pop09,dpop,xguess(n1),x(n1),
#      &   muguess,mudiff,stayp,upper,lower,mutol,
#      &   dist(n1,n1),asian09(n),black09(n),hisp09(n),
#      &   white09(n),asian00(n),delta(n),pstay,rperc(n),
#      &   black00(n),hisp00(n),white00(n),medval(n),
#      &   white_l00(n),white_m00(n),white_h00(n),
#      &   black_l00(n),black_m00(n),black_h00(n),
#      &   asian_l00(n),asian_m00(n),asian_h00(n),
#      &   hisp_l00(n),hisp_m00(n),hisp_h00(n),
#      &   other_l00(n),other_m00(n),other_h00(n),
#      &   white_l09(n),white_m09(n),white_h09(n),
#      &   black_l09(n),black_m09(n),black_h09(n),
#      &   asian_l09(n),asian_m09(n),asian_h09(n),
#      &   hisp_l09(n),hisp_m09(n),hisp_h09(n),
#      &   other_l09(n),other_m09(n),other_h09(n)

pop1 <- vector("numeric", n1)
pop2 <- vector("numeric", n1)
pop1t <- vector("numeric", n)
pop2t <- vector("numeric", n)
pop00 <- c()
pop09 <- c()
stayp <- c()
upper <- c()
xguess <- vector("numeric", n1)
x <- vector("numeric", n1)
muguess <- c()
mudiff <- c()
stayp <- c()
upper <- c()
lower <- c()
mutol <- c()
dist <- matrix(nrow=n1, ncol=n1)
asian09 <- vector("numeric", n)
black09 <- vector("numeric" n)
hisp09 <- vector("numeric", n)
white09 <- vector("numeric", n)
asian00 <- vector("numeric", n)
delta <- vector("numeric", n)
pstay <- c()
rperc <- vector("numeric", n)
black00 <- vector("numeric", n)
hisp00 <- vector("numeric", n)
white00 <- vector("numeric", n)
medval <- vector("numeric", n)
white_l00 <-vector("numeric", n)
white_m00 <- vector("numeric", n)
white_h00 <- vector("numeric", n)
black_l00 <- vector("numeric", n)
black_m00 <- vector("numeric", n)
black_h00 <- vector("numeric", n)
asian_l00 <- vector("numeric", n)
asian_m00 <- vector("numeric", n)
asian_h00 <- vector("numeric", n)
hisp_l00 <- vector("numeric", n)
hisp_m00 <- vector("numeric", n)
hisp_h00 <- vector("numeric", n)
other_l00 <- vector("numeric", n)
other_m00 <- vector("numeric", n)
other_h00 <- vector("numeric", n)
white_l09 <-vector("numeric", n)
white_m09 <- vector("numeric", n)
white_h09 <- vector("numeric", n)
black_l09 <- vector("numeric", n)
black_m09 <- vector("numeric", n)
black_h09 <- vector("numeric", n)
asian_l09 <- vector("numeric", n)
asian_m09 <- vector("numeric", n)
asian_h09 <- vector("numeric", n)
hisp_l09 <- vector("numeric", n)
hisp_m09 <- vector("numeric", n)
hisp_h09 <- vector("numeric", n)
other_l09 <- vector("numeric", n)
other_m09 <- vector("numeric", n)
other_h09 <- vector("numeric", n)

# c
#       common /data1/ pop1,pop2,medval,rperc
#       common /data2/ stayp,pstay
#       common /distance/ dist
# c
#       EXTERNAL fcn1
# c
# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
#       OPEN(10,file='pop_data.txt')
#       REWIND(10)
#       DO j = 1,n
#          READ(10,*) white_l00(j),white_m00(j),white_h00(j),
#      &        black_l00(j),black_m00(j),black_h00(j),
#      &        asian_l00(j),asian_m00(j),asian_h00(j),
#      &        hisp_l00(j),hisp_m00(j),hisp_h00(j),
#      &        other_l00(j),other_m00(j),other_h00(j),
#      &        white_l09(j),white_m09(j),white_h09(j),
#      &        black_l09(j),black_m09(j),black_h09(j),
#      &        asian_l09(j),asian_m09(j),asian_h09(j),
#      &        hisp_l09(j),hisp_m09(j),hisp_h09(j),
#      &        other_l09(j),other_m09(j),other_h09(j)
#       ENDDO
# c

Popdata <- read.delim("inputs/pop_data.txt")

names(Popdata) <- c(
  "white_l00" ,"white_m00" ,"white_h00",
  "black_l00", "black_m00", "black_h00",
  "asian_l00", "asian_m00", "asian_h00",
  "hisp_l00","hisp_m00","hisp_h00",
  "other_l00","other_m00","other_h00",
  "white_l09","white_m09","white_h09",
  "black_l09","black_m09","black_h09",
  "asian_l09","asian_m09","asian_h09",
  "hisp_l09","hisp_m09","hisp_h09",
  "other_l09","other_m09","other_h09"
)



  
#       OPEN(15,file='medval.txt')
#       REWIND(15)
#       DO j = 1,n
#          READ(15,*) medval(j)
#       enddo
# c

Medval <- read.delim("inputs/pop_data.txt")

#       OPEN(18,file='perc_renter.txt')
#       REWIND(18)
#       DO j = 1,n
#          READ(18,*) id3(j),rperc(j)
#       enddo
# c
#       do i = 1,n1
#          do j = 1,n1
#             dist(i,j) = 0.d0
#          enddo
#       enddo

Percrenter <- read.delim("inputs/pop_data.txt")


#       open(20,file='dist_matrix.txt')
#       rewind(20)
#       do i = 1,nn
#          read(20,*) id1(i),id2(i),dist(id1(i),id2(i))
#       enddo
#       do i = 1,n
#          dist(i,n1) = maxd
#       enddo
#       dist(n1,n1) = 0.d0

#Distmatrix <- read.delim("inputs/dist_matrix.txt")

# c
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
white_l <- 0.4174
black_l <- 0.3254
asian_l <- 0.2808
hisp_l <- 0.2892
other_l <- 0.2687

white_m <- 0.4082
black_m <- 0.3933

# c
#       open(30,file='delta_w_l.txt')
#       rewind(30)
#       write(30,*) 'delta_w_l'
#       stayp = 0.4174d0
#       mutol = 0.00000001d0
# c

stayp <- 0.4174
mutol <- 0.00000001

# cccccccccccccccccccccccccc
# c CHANGE POPULATION RACE c
# cccccccccccccccccccccccccc
# c
#       do i = 1,n
#          pop1t(i) = white_l00(i)+0.001d0
#          pop2t(i) = white_l09(i)+0.001d0
#          pop00 = pop00+pop1t(i)
#          pop09 = pop09+pop2t(i)
#       enddo

pop1t <- Popdata$white_l00 + 0.001
pop2t <- Popdata$white_l09 + 0.001

pop00 <- pop1t
pop09 <- pop2t

#  555  format(i9,3f21.7)
#       dpop = pop09-pop00
#       write(*,*) 'dpop =',dpop
#       do i = 1,n
#          pop1(i) = pop1t(i)
#          pop2(i) = pop2t(i)
#       enddo
#       pop1(n1) = 2.d0*dabs(dpop)
#       pop2(n1) = pop1(n1)-dpop

# dabs : absolute difference?

dpop <- pop09 - pop00
pop1[-n1] <- pop1t
pop2[-n1] <- pop2t

pop1[n1] <- 


# c
#       do j = 1,n1
#          xguess(j) = 0.d0
#       enddo
# c
# c
#       lower = -0.15d0
#       upper = 0.1d0
#       do k1 = 1,maxit
#          muguess = (upper+lower)/2.d0
#          CALL fcn1(xguess,muguess,x,mudiff)
#          if (mudiff.gt.0.d0) lower = muguess
#          if (mudiff.le.0.d0) upper = muguess
#          do j = 1,n1
#             xguess(j) = x(j)
#          enddo
#          write(*,40) k1,muguess,mudiff,lower,upper
#  40      format(i6,5f9.5)
#          if(dabs(mudiff).le.mutol) go to 120
#       enddo
#  120  continue
#       do j = 1,n
#          delta(j) = x(j)
#       enddo
# c
#       do j = 1,n
#          write(30,*) delta(j)
# c         write(*,*) j,delta(j)
#       enddo
# c
#       write(*,*) 'delta(n1) =',x(n1)
#       write(*,*) 'pstay     =',pstay
#       END
# c


# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
#       SUBROUTINE fcn1(x1,muguess,xf,mudiff)
# c
# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
#       USE numerical_libraries
# c
#       PARAMETER (tolc=0.000001d0)
#       PARAMETER (n=1921,n1=n+1)
#       INTEGER i,j
#       DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
#      &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
#      &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
#      &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
#      &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
#      &   rperc(n),mc
# c
#       common /data1/ pop1,pop2,medval,rperc
#       common /data2/ stayp,pstay
#       common /distance/ dist
# c
# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c
# c     
fcn1 <- function(x1, muguess, xf, mudiff){
#   mu = muguess
  mu <- muguess
  
  #   totpop = 0.d0
  #   DO i = 1,n1
  #   totpop = totpop+pop2(i)
  #   ENDDO
  totpop <- sum(pop2)
  
# #   c
# #   300  continue
# #   do j = 1,n1
# #   do k = 1,n1
# #   mc = 0.d0
#   if ((j.ne.k).and.(j.ne.n1).and.
#       &         (k.ne.n1)) mc = (2910.d0+(0.03d0*(1.d0-
#                                                     &         rperc(j))*medval(j))+(0.03d0*(1.d0-
#                                                                                               &         rperc(k))*medval(k)))/24.556d0
#   if ((j.ne.k).and.((j.eq.n1).or.
#                     &         (k.eq.n1))) mc = 528.71d0
#   util(k,j) = x1(k)-x1(j)-mu*mc
#   enddo
#   enddo
#   do j = 1,n1
#   denom(j) = 0.d0
#   do k = 1,n1
#   denom(j) = denom(j)+dexp(util(k,j))
#   enddo
#   enddo
#   do i = 1,n1
#   rhs(i) = 0.d0
#   do j = 1,n1
#   numer = dexp(util(i,j))
#   rhs(i) = rhs(i)+(numer/denom(j))*
#     &               pop1(j)
#   enddo
#   enddo
#   stay = 0.d0
#   do i = 1,n
#   stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
#   enddo
#   stay = stay/(totpop-pop2(n1))
#   c
#   pstay = 0.d0
#   do i = 1,n
#   pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
#   enddo
#   c
#   do i = 1,n1
#   share(i) = pop2(i)/totpop
#   pshare(i) = rhs(i)/totpop
#   enddo
#   c
#   xavg = 0.d0
#   DO j = 1,n1
#   xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
#   c         write(*,*) j,share(j),pshare(j)
#   xavg = xavg+(xf(j)/(1.d0*n1))
#   ENDDO
#   225  format(i10,f15.6,f15.6)
#   c
#   temp = xavg
#   DO j = 1,n1
#   xf(j) = xf(j)-temp
#   enddo
#   c      
#   DO j = 1,n1
#   test(j) = dabs(xf(j)-x1(j))
#   IF (test(j).gt.tolc) then
#   DO k = 1,n1
#   x1(k) = xf(k)
#   ENDDO
#   GO TO 300
#   ENDIF
#   ENDDO
#   mudiff = stayp - pstay
#   c
#   END
#   c
#   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#   
#   
}
