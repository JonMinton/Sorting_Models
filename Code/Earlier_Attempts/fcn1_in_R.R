#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#c
#SUBROUTINE fcn1(x1,muguess,xf,mudiff)
fcn1 <- function(x1, muguess, xf, mudiff){
  subfcn1 <- function(){
    # 300  continue
    # do j = 1,n1
    # do k = 1,n1
    for (j in 1:n1){
      for (k in 1:n1){
        
        # mc = 0.d0
        mc <- 0
        
        #     if ((j.ne.k).and.(j.ne.n1).and.(k.ne.n1)) 
        #       mc = (2910.d0+(0.03d0*(1.d0-rperc(j))*medval(j))+(0.03d0*(1.d0-rperc(k))*medval(k)))/24.556d0
        if (j != k) & (j != n1) & (k!=n1){
          mc <- (2910 + (0.03 * (1 - rperc[j]) * medval[j]) + (0.03 * (1 - rperc[k])*medval[k]))/24.556
        }
        
        #     if ((j.ne.k).and.((j.eq.n1).or.(k.eq.n1))) 
        #       mc = 528.71d0
        if ((j!=k) & ((j ==n1) | (k ==n1))){      
          mc <- 528.71
        }
        
        #    util(k,j) = x1(k)-x1(j)-mu*mc # utility of going from area k to area j    
        util[k, j] <<- x1[k] - x1[j] - mu*mc
        
        # enddo
        # enddo    
      }
    }
    
    # 
    #do j = 1,n1
    for (j in 1:n1){
      #   denom(j) = 0.d0
      denom[j] <- 0  
      
      #   do k = 1,n1
      for (k in 1:n1){
        # denom(j) = denom(j)+dexp(util(k,j)) # exponential of util
        denom[j] <- denom[j] + exp(util[k,j])
        # enddo
        # enddo      
      }
    }
    
    # denominator for area j is sum of exponentials for area utilities 
    # for all other js  
    
    # do i = 1,n1
    
    for (i in 1:n1){
      #   rhs(i) = 0.d0
      rhs[i] <- 0
      #   do j = 1,n1
      for (j in 1:n1){
        # numer = dexp(util(i,j))
        numer <- exp(util[i,j])
        # rhs(i) = rhs(i)+(numer/denom(j))*pop1(j)
        rhs[i] <- rhs(i) + (numer/denom[j]) * pop1[j]
        #     enddo
        #     enddo
        
      }
    }
    
    # 
    # stay = 0.d0
    stay <- 0
    
    # do i = 1,n
    for (i in 1:n){
      #   stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
      stay <- stay + (exp(util[i,j])/denom[i]*pop1[i]  
                      #   enddo  
    }
    
    # stay = stay/(totpop-pop2(n1))
    stay <- stay/(totpop-pop2[n1])
    # c
    # pstay = 0.d0
    pstay <- 0
    # do i = 1,n
    for (i in 1:n){
      #   pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
      pstay <- pstay + (exp(util[i,j])/denom[i])/(1 * n)  
      #   enddo
      
    }
    # c
    # do i = 1,n1
    for (i in 1:n1){
      #   share(i) = pop2(i)/totpop
      share[i] <<- pop2[i]/totpop
      #   pshare(i) = rhs(i)/totpop
      pshare[i] <<- rhs[i]/totpop
    }
    
    # enddo
    # c
    # xavg = 0.d0
    xavg <<- 0
    xf <<- rep(NA, n1)
    # DO j = 1,n1
    for (j in 1:n1){
      # Contraction mapping : adding difference in log shares
      # xf: updated value of x
      # x1 : previous value of x
      # normalised by subtracting mean from x
      
      # xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
      xf[j] <<- x1[j] + log(share[j]) - log(pshare[j])
      
      #   c         write(*,*) j,share(j),pshare(j)
      
      #   xavg = xavg+(xf(j)/(1.d0*n1))
      xavg <<- xavg + xf[j]/(1*n1)
      
      #   ENDDO  
    }
    
    
    
    # 225  format(i10,f15.6,f15.6)
    # c
    # temp = xavg 
    # normalisation (needs to be done)
    temp <<- xavg
    
    # DO j = 1,n1
    for (j in 1:n1){
      #   xf(j) = xf(j)-temp # normalised updated x estimate
      xf[j] <<- xf[j] - temp
      
    }
    # enddo
    # c     
    test <<- rep(NA, n1)
    # DO j = 1,n1
    for (j in 1:n1){
      #   test(j) = dabs(xf(j)-x1(j))
      test[j] <<- abs(xf[j] - x1[j])
      #   IF (test(j).gt.tolc) then
      if (test[j] > tolc){
        #     DO k = 1,n1
        for (j in 1:n1){
          # x1(k) = xf(k)
          x1[k] <<- xf[k]
          # ENDDO
        }
        # go back to subfunction
        #     GO TO 300
        
        subfcn1()
        # ENDIF
      }
      # ENDDO
    }
  }
  
  
  #c
  #ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  #c
  #USE numerical_libraries
  #c
  #PARAMETER (tolc=0.000001d0)
  #PARAMETER (n=1921,n1=n+1)
  #INTEGER i,j
  # DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
  # &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
  # &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
  # &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
  # &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
  # &   rperc(n),mc
  # c
  # common /data1/ pop1,pop2,medval,rperc
  # common /data2/ stayp,pstay
  # common /distance/ dist
  # c
  # ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  # c
  # c     
  # mu = muguess
  mu <<- muguess
  
  denom <<- rep(NA, n1)
  rhs <<- rep(NA, n1)
  
  util <<- matrix(NA, nrow=n1, ncol=n1)
  
  # totpop = 0.d0
  totpop <<- 0
  for (i in 1:n1){
    
    # DO i = 1,n1
    # totpop = totpop+pop2(i)
    totpop <<- totpop + pop2[i]
    # Total in the second period
    # Shares are shares based on sum of pop2
    # ENDDO
  }
  # c
  subfcn1()
  
  #   mudiff = stayp - pstay
  mudiff <- stayp - pstay
  return(mudiff)
  
  # END  
}

# Check all values are within tolerance
# if all values are within tolerance, end function, and 
# return mudiff

# c
# c
# ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# 
# }