cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      PROGRAM estimate
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      USE numerical_libraries
c   
      PARAMETER (n=1921,n1=n+1,maxit=1000,nn=(n-1)*(n-1))
      INTEGER i,j,tract(n),id1(nn),id2(nn),id3(n)
      DOUBLE PRECISION pop1(n1),pop2(n1),pop1t(n),
     &   pop2t(n),pop00,pop09,dpop,xguess(n1),x(n1),
     &   muguess,mudiff,stayp,upper,lower,mutol,
     &   dist(n1,n1),asian09(n),black09(n),hisp09(n),
     &   white09(n),asian00(n),delta(n),pstay,rperc(n),
     &   black00(n),hisp00(n),white00(n),medval(n),
     &   white_l00(n),white_m00(n),white_h00(n),
     &   black_l00(n),black_m00(n),black_h00(n),
     &   asian_l00(n),asian_m00(n),asian_h00(n),
     &   hisp_l00(n),hisp_m00(n),hisp_h00(n),
     &   other_l00(n),other_m00(n),other_h00(n),
     &   white_l09(n),white_m09(n),white_h09(n),
     &   black_l09(n),black_m09(n),black_h09(n),
     &   asian_l09(n),asian_m09(n),asian_h09(n),
     &   hisp_l09(n),hisp_m09(n),hisp_h09(n),
     &   other_l09(n),other_m09(n),other_h09(n)
c
      common /data1/ pop1,pop2,medval,rperc
      common /data2/ stayp,pstay
      common /distance/ dist
c
      EXTERNAL fcn1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      OPEN(10,file='pop_data.txt')
      REWIND(10)
      DO j = 1,n
         READ(10,*) white_l00(j),white_m00(j),white_h00(j),
     &        black_l00(j),black_m00(j),black_h00(j),
     &        asian_l00(j),asian_m00(j),asian_h00(j),
     &        hisp_l00(j),hisp_m00(j),hisp_h00(j),
     &        other_l00(j),other_m00(j),other_h00(j),
     &        white_l09(j),white_m09(j),white_h09(j),
     &        black_l09(j),black_m09(j),black_h09(j),
     &        asian_l09(j),asian_m09(j),asian_h09(j),
     &        hisp_l09(j),hisp_m09(j),hisp_h09(j),
     &        other_l09(j),other_m09(j),other_h09(j)
      ENDDO
c
      OPEN(15,file='medval.txt')
      REWIND(15)
      DO j = 1,n
         READ(15,*) medval(j)
      enddo
c
      OPEN(18,file='perc_renter.txt')
      REWIND(18)
      DO j = 1,n
         READ(18,*) id3(j),rperc(j)
      enddo
c
      do i = 1,n1
         do j = 1,n1
            dist(i,j) = 0.d0
         enddo
      enddo
      open(20,file='dist_matrix.txt')
      rewind(20)
      do i = 1,nn
         read(20,*) id1(i),id2(i),dist(id1(i),id2(i))
      enddo
      do i = 1,n
         dist(i,n1) = maxd
      enddo
      dist(n1,n1) = 0.d0
c
cccccccccccccccccccccccccccccccc
c Stay-put % by race:          c
c    - white_l = 0.4174        c
c    - black_l = 0.3254        c
c    - asian_l = 0.2808        c
c    - hisp_l = 0.2892          c
c    - other_l = 0.2687        c
c                              c
c    - white_m = 0.4082        c
c    - black_m = 0.3933        c
c    - asian_m =         c
c    - hisp_m =          c
c    - other_m =         c
c                              c
c    - white_h =         c
c    - black_h =         c
c    - asian_h =         c
c    - hisp_h =          c
c    - other_h =         c
cccccccccccccccccccccccccccccccc
c
      open(30,file='delta_w_l.txt')
      rewind(30)
      write(30,*) 'delta_w_l'
      stayp = 0.4174d0
      mutol = 0.00000001d0
c
cccccccccccccccccccccccccc
c CHANGE POPULATION RACE c
cccccccccccccccccccccccccc
c
      do i = 1,n
         pop1t(i) = white_l00(i)+0.001d0
         pop2t(i) = white_l09(i)+0.001d0
         pop00 = pop00+pop1t(i)
         pop09 = pop09+pop2t(i)
      enddo
 555  format(i9,3f21.7)
      dpop = pop09-pop00
      write(*,*) 'dpop =',dpop
      do i = 1,n
         pop1(i) = pop1t(i)
         pop2(i) = pop2t(i)
      enddo
      pop1(n1) = 2.d0*dabs(dpop)
      pop2(n1) = pop1(n1)-dpop
c
      do j = 1,n1
         xguess(j) = 0.d0
      enddo
c
c
      lower = -0.15d0
      upper = 0.1d0
      do k1 = 1,maxit
         muguess = (upper+lower)/2.d0
         CALL fcn1(xguess,muguess,x,mudiff)
         if (mudiff.gt.0.d0) lower = muguess
         if (mudiff.le.0.d0) upper = muguess
         do j = 1,n1
            xguess(j) = x(j)
         enddo
         write(*,40) k1,muguess,mudiff,lower,upper
 40      format(i6,5f9.5)
         if(dabs(mudiff).le.mutol) go to 120
      enddo
 120  continue
      do j = 1,n
         delta(j) = x(j)
      enddo
c
      do j = 1,n
         write(30,*) delta(j)
c         write(*,*) j,delta(j)
      enddo
c
      write(*,*) 'delta(n1) =',x(n1)
      write(*,*) 'pstay     =',pstay
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE fcn1(x1,muguess,xf,mudiff)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      USE numerical_libraries
c
      PARAMETER (tolc=0.000001d0)
      PARAMETER (n=1921,n1=n+1)
      INTEGER i,j
      DOUBLE PRECISION totpop,pop1(n1),pop2(n1),
     &   denom(n1),x1(n1),rhs(n1),numer,medval(n),
     &   share(n1),pshare(n1),xavg,xf(n1),test(n1),
     &   dum(n1,n1),mu,util(n1,n1),temp,stay,mudiff,
     &   muguess,stayp,dist(n1,n1),ddum,pstay,ddum1,
     &   rperc(n),mc
c
      common /data1/ pop1,pop2,medval,rperc
      common /data2/ stayp,pstay
      common /distance/ dist
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
      mu = muguess
      totpop = 0.d0
      DO i = 1,n1
         totpop = totpop+pop2(i)
      ENDDO
c
 300  continue
      do j = 1,n1
         do k = 1,n1
            mc = 0.d0
            if ((j.ne.k).and.(j.ne.n1).and.
     &         (k.ne.n1)) mc = (2910.d0+(0.03d0*(1.d0-
     &         rperc(j))*medval(j))+(0.03d0*(1.d0-
     &         rperc(k))*medval(k)))/24.556d0
            if ((j.ne.k).and.((j.eq.n1).or.
     &         (k.eq.n1))) mc = 528.71d0
            util(k,j) = x1(k)-x1(j)-mu*mc
         enddo
      enddo
      do j = 1,n1
         denom(j) = 0.d0
         do k = 1,n1
            denom(j) = denom(j)+dexp(util(k,j))
         enddo
      enddo
      do i = 1,n1
         rhs(i) = 0.d0
         do j = 1,n1
            numer = dexp(util(i,j))
            rhs(i) = rhs(i)+(numer/denom(j))*
     &               pop1(j)
         enddo
      enddo
      stay = 0.d0
      do i = 1,n
         stay = stay + (dexp(util(i,i))/denom(i))*pop1(i)
      enddo
      stay = stay/(totpop-pop2(n1))
c
      pstay = 0.d0
      do i = 1,n
         pstay = pstay + (dexp(util(i,i))/denom(i))/(1.d0*n)
      enddo
c
      do i = 1,n1
         share(i) = pop2(i)/totpop
         pshare(i) = rhs(i)/totpop
      enddo
c
      xavg = 0.d0
      DO j = 1,n1
         xf(j) = x1(j) + (dlog(share(j)) - dlog(pshare(j)))
c         write(*,*) j,share(j),pshare(j)
         xavg = xavg+(xf(j)/(1.d0*n1))
      ENDDO
 225  format(i10,f15.6,f15.6)
c
      temp = xavg
      DO j = 1,n1
         xf(j) = xf(j)-temp
      enddo
c      
      DO j = 1,n1
         test(j) = dabs(xf(j)-x1(j))
         IF (test(j).gt.tolc) then
            DO k = 1,n1
               x1(k) = xf(k)
            ENDDO
            GO TO 300
         ENDIF
      ENDDO
      mudiff = stayp - pstay
c
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
