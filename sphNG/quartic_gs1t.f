      SUBROUTINE QUARTIC_GS1T(u1term,u0term,uold,soln,moresweep)
      
      REAL tiny
      PARAMETER (tiny=1.0E-30)

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'

      REAL  E,U,cv,dt,planck,rho,s5,s6,s7,y1,ua,ub1,uc1,ub2,uc2,kappa
      INTEGER  sooty,im,pid,inpt,I,rtst
      REAL  a,b,d,f,ue0,lp,p,q,r,s,t,a0,a1,a2,a3,a4,ub,uc,t1,t2,yy
      REAL  z1,z2,z3,z4
      COMPLEX  ca,cb,cc,cd,ce,cf
      COMPLEX  p1,p2,p3,y2,y3,y,t0
!      COMPLEX,DIMENSION(4)  roots
!      COMPLEX,DIMENSION(5)  coeffs
!      REAL,DIMENSION(4)  real_roots
      CHARACTER*1  yn
      COMPLEX  r3,t3,q3
      REAL  tmp,ueo,uen,tb,banana,gamma,c1,c2,soln,tst1,tst2,tst
      REAL  tr,tm,tsoln1,tsoln2,tmin,tmax,y1a,fac,divv,c4
      REAL  tsoln3,tsoln4,c3,c5,tsoln
      REAL  d1,d2,d3,d4,e0,e1,e2,e3,e4,d5,quantity1,biggest_term
      REAL  f1,f2,f3,f4,f5,g0,g1,g2,g3,g4,laeg,la1,lg1,le1,vhgr,vlwr
      LOGICAL  swapg,swape,NR,moresweep
      DIMENSION z1(2),z2(2),z3(2),z4(2)

                                !QUARTIC4 - solver for trapezoidal T^4
      lightspeed = c /udist * utime
      uradconst = radconst / uergcc


      GOTO 300

      
c      PRINT *,"Entered Q4 with ",sooty," sooty"
      

      ue0=kappa/uradconst*E-planck/cv**4*U**4

c      PRINT *,"Original value of x:",ue0

      d1=kappa*rho/uradconst*E-1.0/2.0*dt/sooty*rho/uradconst*uradconst*
     $     lightspeed*(kappa*rho*E/uradconst-planck/cv**4*U**4)
      d2=planck/cv**4
      d3=U+0.5*dt/sooty*uradconst*lightspeed*(kappa/uradconst*E-planck/
     $     cv**4*U**4)
      d4=0.5*kappa/uradconst*dt/sooty*uradconst*lightspeed*rho
      d5=1.0/2.0*dt/sooty*uradconst*lightspeed

      
      c1=rho*E-1.0/2.0*dt/sooty*uradconst*lightspeed*rho*
     $     (kappa/uradconst*rho*E-planck/cv**4*U**4)
      c2=planck/cv**4*uradconst/kappa
      c3=U+0.5*dt/sooty*uradconst*lightspeed*(kappa/uradconst*rho*E-
     $     planck/cv**4*U**4)
      c4=0.5*kappa/uradconst*dt/sooty*uradconst*lightspeed*rho
      c5=0.5*dt/sooty*uradconst*lightspeed*kappa/uradconst

      f1=uradconst*rho*E/kappa-1.0/2.0*dt/sooty*rho*uradconst**2/kappa*
     $     lightspeed*(kappa/uradconst*rho*E-planck/cv**4*U**4)
      f2=planck/cv**4*uradconst**2/kappa*2
      f3=U+0.5*dt/sooty*uradconst*lightspeed*(kappa/uradconst*
     $     rho*E-planck/cv**4*U**4)
      f4=0.5*kappa/uradconst*dt/sooty*uradconst*lightspeed*rho
      f5=0.5*dt/sooty*uradconst*lightspeed*kappa**2/uradconst**2


 300  CONTINUE
      
                                ! a4=1.0
                                ! a3=4.0*c2*c3*c5**3/c2/c5**4
                                ! a2=6.0*c2*c3**2*c5**2/c2/c5**4
                                ! a1=(4.0*c2*c3**3*c5+1.0+c4)/c2/c5**4
                                ! a0=(c2*c3**4-c1)/c2/c5**4

                                !Simplified versions
      a4=1.0
      a3=0.0
      a2=0.0
      a1=u1term
      a0=u0term

      e4=1.0
      e3=0.0
      e2=0.0
      e1=u1term
      e0=u0term

      g4=1.0
      g3=0.0
      g2=0.0
      g1=u1term
      g0=u0term
      
      
      NR=.FALSE. 
      
      z1(2)=0.0
      z2(2)=0.0
      z3(2)=0.0
      z4(2)=0.0

      swapg=.FALSE.
      swape=.FALSE.

      la1=ABS(LOG10(ABS(a1))-10)
      le1=ABS(LOG10(ABS(e1))-10)
      lg1=ABS(LOG10(ABS(g1))-10)
      
      laeg=MIN(la1,le1,lg1)
      test1=a1
      test2=e1
      test3=g1
      
      IF(laeg.EQ.lg1) THEN
         swapg=.TRUE.
         
         a0=g0
         a1=g1
         a2=g2
         a3=g3
         a4=g4

c         vhgr=1.0d1*uradconst**2/kappa**2*(kappa/uradconst*rho*E-
c     $     planck/cv**4*U**4)
c         vlwr=1.0d-1*uradconst**2/kappa**2*(kappa/uradconst*rho*E-
c     $     planck/cv**4*U**4)
      ELSE IF(laeg.EQ.le1) THEN
         a0=e0
         a1=e1
         a2=e2
         a3=e3
         a4=e4
         swape=.TRUE.

c         vhgr=1.0d1*kappa/uradconst*E-planck/cv**4*U**4
c         vlwr=1.0d-1*kappa/uradconst*E-planck/cv**4*U**4
      ELSE
c      vhgr=1.0d1*uradconst/kappa*(kappa/uradconst*E-planck/cv**4*U**4)
c      vlwr=1.0d-1*uradconst/kappa*(kappa/uradconst*E-planck/cv**4*U**4)

      END IF

c      tm=u/cv
c      tr=(e*rho/uradconst)**0.25

c      tmin=MIN(tm,tr)
c      tmax=MAX(tm,tr)

c     !     Real root of equation:
c     !     y^3 -4 a0 y - a1^2 = 0
      
      quantity1=-54*a2*a1**3*a3+12*a3**2*a0*a2**3-3*a1**2*a3**2*a2**2-
     $     432*a2*a0*a1**2+18*a1**2*a3**2*a0+576*a1*a3*a0**2-432*a2*
     $     a0**2*a3**2+240*a1*a3*a2**2*a0-54*a2*a1*a3**3*a0+384*a2**2*
     $     a0**2-48*a0*a2**4+12*a1**2*a2**3+81*a3**4*a0**2-768*a0**3+
     $     81*a1**4+12*a1**3*a3**3

      biggest_term=MAX(ABS(-54*a2*a1**3*a3),ABS(12*a3**2*a0*a2**3),
     $     ABS(-3*a1**2*a3**2*a2**2),ABS(-432*a2*a0*a1**2),
     $     ABS(18*a1**2*a3**2*a0),ABS(576*a1*a3*a0**2),
     $     ABS(-432*a2*a0**2*a3**2),ABS(240*a1*a3*a2**2*a0),
     $     ABS(-54*a2*a1*a3**3*a0),ABS(384*a2**2*a0**2),
     $     ABS(-48*a0*a2**4),ABS(12*a1**2*a2**3),ABS(81*a3**4*a0**2),
     $     ABS(-768*a0**3),ABS(81*a1**4),ABS(12*a1**3*a3**3))

      !PRINT *,quantity1,biggest_term
      IF(quantity1.LT.0.0.AND.ABS(quantity1)/biggest_term.LT.1d-12) THEN
!        PRINT *,"q1 ",quantity1,biggest_term
         quantity1=0.0
         PRINT *,"q1 ",quantity1
!        moresweep=.TRUE.
      ELSE IF(quantity1.LT.0.0) THEN
         PRINT *,"QUARTIC4: Quantity1 is negative. "
         PRINT *,"Quantity1:",quantity1,biggest_term
         PRINT *,"Returning to TRAP with moresweep2=.TRUE."
         moresweep=.TRUE.
         GOTO 541
                                !    RETURN
         STOP
      END IF
      
!     PRINT *,"Getting y1",quantity1
      y1 = (-36*a2*a1*a3-288*a2*a0+108*a1**2+108*a3**2*a0+8*a2**3+12*
     $     SQRT(quantity1))**(1.D0/3.D0)/6.0-6*(a1*a3/3-4.D0/3.D0*a0-
     $     a2**2/9.0)/(-36*a2*a1*a3-288*a2*a0+108*a1**2+108*a3**2*a0+
     $     8*a2**3+12*SQRT(quantity1))**(1.D0/3.D0)+a2/3.0
           
      z1(2)=0.0
      z2(2)=0.0
      z3(2)=0.0
      z4(2)=0.0
      

      z1(1)=0.0
      z2(1)=0.0
      z3(1)=0.0
      z4(1)=0.0
      
                                !Solution to quartic
                                !This is solution to two quadratics
                                !v^2+(a2/2... etc)

      ub=(a3**2/4.0+y1-a2)

!     PRINT *,"Ub:",ub
      
      IF(ub.LT.0.0) THEN
         PRINT *,"QUARTIC4: Error, imaginary co-eff b to quadratic"
         PRINT *,test1,test2,test3
         PRINT *,u1term,u0term
         PRINT *,"Returning to TRAPIMPL with moresweep2=.TRUE."
         moresweep=.TRUE.
         RETURN
      END IF

      uc=((y1/2.0)**2-a0)
!     PRINT *,"Uc:",uc
      
      IF(uc.LT.0.0) THEN
         PRINT *,"QUARTIC4: Error, imaginary co-eff c to quadratic"
         STOP
      END IF

      ub1=a3/2.0+SQRT(ub)
      ub2=a3/2.0-SQRT(ub)
      IF (a1.LT.0.0 .AND. a0.GT.0.0) THEN
         uc1=y1/2.0+SQRT(uc)
         uc2=y1/2.0-SQRT(uc)
      ELSE
         uc1=y1/2.0-SQRT(uc)
         uc2=y1/2.0+SQRT(uc)
      ENDIF

c      if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3q',
c     &     y1,uc,ub,a2,a3,a1,a0

                                ! PRINT *,"ubc",ub1,ub2,uc1,uc2

c      if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3'

      IF(ub2**2-4.0*uc2.GT.0.0.AND.ub1**2-4.0*uc1.LT.0.0) THEN
c         if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3a',
c     &        ub2**2-4.0*uc2,ub1**2-4.0*uc1,ub2,uc2,ub1,uc1

         IF (ABS(a3).GE.tiny) THEN
            IF(ABS((a3/2.0)-SQRT(ub))/ABS(a3/2.0).LT.1d-6) THEN
         PRINT *,"QUARTIC4: Error, big - big / big too big for co-eff b"
               STOP
            END IF
         ENDIF

            z2(1)=0.5*((-ub2)-SQRT(ub2**2-4.0*uc2))
            z1(1)=0.5*((-ub2)+SQRT(ub2**2-4.0*uc2))

            z2(2)=1
            z1(2)=1
      ELSE IF(ub2**2-4.0*uc2.LT.0.0.AND.ub1**2-4.0*uc1.GT.0.0) THEN
c         if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3b',
c     &        ub2**2-4.0*uc2,ub1**2-4.0*uc1,ub2,uc2,ub1,uc1

         IF(ABS(2.0*a0/(y1**2)).LT.1d-6) THEN

            IF (ABS(a3).LT.tiny) THEN
               PRINT *,"a3 < tiny at position 1 ",a3
               STOP
            ELSE
            IF(ABS((2.0*a0/y1)/(((a3/2.0)+SQRT((a3**2/4.0)+y1-a2))**2)).
     $           GT.1d-6) THEN
            PRINT *,"QUARTIC4: Second Taylor expansion no longer valid"
               PRINT *,"QUARTIC4: Value is ",ABS((a0/y1)/((a3/2.0)+
     $              SQRT((a3**2/4.0)+y1-a2)))
               STOP
            END IF


            z4(1)=-1.0*((a0/y1)/((a3/2.0)+SQRT((a3**2/4.0)+y1-a2)))
            ENDIF

         ELSE
            z4(1)=0.5*((-ub1)+SQRT(ub1**2-4.0*uc1))

         END IF

                                !   PRINT *,"z34"
         z3(1)=0.5*((-ub1)-SQRT(ub1**2-4.0*uc1))
         z3(2)=1
         z4(2)=1
         

      ELSE  IF(ub2**2-4.0*uc2.LT.0.0.AND.ub1**2-4.0*uc1.LT.0.0) THEN
c         if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3c',
c     &        ub2**2-4.0*uc2,ub1**2-4.0*uc1,ub2,uc2,ub1,uc1


     !!  NUMERICAL SOLUTION IF ONLY IMAGINARY ARE RETURNED ANALYTICALLY !!
     !! 
         
         PRINT *,"QUARTIC4: All imaginary roots for quartic"
                                !    PRINT *,vhgr
c     !     PRINT *,vlwr
         NR=.TRUE.
                                !   PRINT *,a0,a1,a2,a3,a4
c      CALL RTSAFE(soln,vlwr,vhgr,1d-8,a0,a1,a2,a3,a4)
c      PRINT *,"soln:",soln
         STOP
      ELSE
c         if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 3c',
c     &        ub2**2-4.0*uc2,ub1**2-4.0*uc1,ub2,uc2,ub1,uc1

         IF (ABS(a3).GE.tiny) THEN
            IF(ABS((a3/2.0)-SQRT(ub))/ABS(a3/2.0).LT.1d-6) THEN
         PRINT *,"QUARTIC4: Error, big - big / big too big for co-eff b"
               STOP
            END IF
         ENDIF

         IF(ABS(2.0*a0/(y1**2)).LT.1d-6) THEN

            IF (ABS(a3).GE.tiny) THEN
            IF((ABS(2.0*a0/y1)/(((a3/2.0)+SQRT((a3**2/4.0)+y1-a2))**2)).
     $           GT.1d-6) THEN
            PRINT *,"QUARTIC4: Second Taylor expansion no longer valid"
               PRINT *,"QUARTIC4: Value is ",ABS((a0/y1)/((a3/2.0)+
     $           SQRT((a3**2/4.0)+y1-a2)))
               STOP
            END IF

            z4(1)=-1.0*((a0/y1)/((a3/2.0)+SQRT((a3**2/4.0)+y1-a2)))
            ELSE
               PRINT *,"QUARTIC4: a3<tiny at position 2 ",a3
            ENDIF

         ELSE
            z4(1)=0.5*((-ub1)+SQRT(ub1**2-4.0*uc1))
         END IF

         z3(1)=0.5*((-ub1)-SQRT(ub1**2-4.0*uc1))
         
         z3(1)=0.5*((-ub1)+SQRT(ub1**2-4.0*uc1))
         z2(1)=0.5*((-ub2)-SQRT(ub2**2-4.0*uc2))
         z1(1)=0.5*((-ub2)+SQRT(ub2**2-4.0*uc2))

         z4(2)=1
         z3(2)=1
         z2(2)=1
         z1(2)=1
      
      END IF


!     PRINT *,z1(1),z2(1),z3(1),z4(1)


!     END IF
      
 150  GOTO 555
c            if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 4'
      IF(swapg) THEN
      z1(1)=kappa**2/uradconst**2*z1(1)
      z2(1)=kappa**2/uradconst**2*z2(1)
      z3(1)=kappa**2/uradconst**2*z3(1)
      z4(1)=kappa**2/uradconst**2*z4(1)
      IF(NR) soln=kappa**2/uradconst**2*soln
      ELSE IF(swape) THEN
      z1(1)=z1(1)
      z2(1)=z2(1)
      z3(1)=z3(1)
      z4(1)=z4(1)
      IF(NR) soln=soln
      ELSE
      z1(1)=kappa/uradconst*z1(1)
      z2(1)=kappa/uradconst*z2(1)
      z3(1)=kappa/uradconst*z3(1)
      z4(1)=kappa/uradconst*z4(1)
      IF(NR) soln=kappa/uradconst*soln
      END IF
 555  CONTINUE

cC$OMP CRITICAL(quart)
c      PRINT *,"Quart ",z1(1),z2(1),z3(1),z4(1)
c      PRINT *,u1term,u0term
c      PRINT *,"     ",quantity1,biggest_term
c      PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
cC$OMP END CRITICAL(quart)


!     GOTO 160  
c            if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 5',
c     &     z1(2),z2(2),z3(2),z4(2)
      IF(z1(2).EQ.1.AND.z2(2).EQ.1.AND.z3(2).EQ.0.AND.z4(2).EQ.0) THEN
        
         IF(z1(1).GT.0.0 .AND. z2(1).LE.0.0) THEN
            soln = z1(1)
         ELSEIF (z1(1).LE.0.0 .AND. z2(1).GT.0.0) THEN
            soln = z2(1)
         ELSEIF (z1(1).LE.0.0 .AND. z2(1).LE.0.0) THEN
            PRINT *,"Failed 1 ",z1(1),z2(1),z3(1),z4(1),uold
            PRINT *,u1term,u0term
            PRINT *,"     ",quantity1,biggest_term
            PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
            STOP
         ELSEIF (ABS(z1(1)-uold).GT.ABS(z2(1)-uold)) THEN
            soln = z2(1)
            IF (ABS((soln-uold)/uold).GT.1.0) THEN
               PRINT *,"Change big 1a ",z1(1),z2(1),z3(1),z4(1),uold
               PRINT *,u1term,u0term
               PRINT *,"     ",quantity1,biggest_term
               PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
               STOP
            ENDIF
         ELSE
            soln = z1(1)
            IF (ABS((soln-uold)/uold).GT.1.0) THEN
               PRINT *,"Change big 1b ",z1(1),z2(1),z3(1),z4(1),uold
               PRINT *,u1term,u0term
               PRINT *,"     ",quantity1,biggest_term
               PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
               STOP
            ENDIF
         ENDIF

      ELSE IF(z3(2).EQ.1.AND.z4(2).EQ.1.AND.z1(2).EQ.0.AND.z2(2).EQ.0)
     $     THEN
         
         IF(z3(1).GT.0.0 .AND. z4(1).LE.0.0) THEN
            soln = z3(1)
         ELSEIF (z3(1).LE.0.0 .AND. z4(1).GT.0.0) THEN
            soln = z4(1)
         ELSEIF (z3(1).LE.0.0 .AND. z4(1).LE.0.0) THEN
            PRINT *,"Failed 2 ",z1(1),z2(1),z3(1),z4(1),uold
            PRINT *,u1term,u0term
            PRINT *,"     ",quantity1,biggest_term
            PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
            STOP
         ELSEIF (ABS(z3(1)-uold).GT.ABS(z4(1)-uold)) THEN
            soln = z4(1)
            IF (ABS((soln-uold)/uold).GT.1.0) THEN
               PRINT *,"Change big 2a ",z1(1),z2(1),z3(1),z4(1),uold
               PRINT *,u1term,u0term
               PRINT *,"     ",quantity1,biggest_term
               PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
               STOP
            ENDIF
         ELSE
            soln = z3(1)
            IF (ABS((soln-uold)/uold).GT.1.0) THEN
               PRINT *,"Change big 2b ",z1(1),z2(1),z3(1),z4(1),uold
               PRINT *,u1term,u0term
               PRINT *,"     ",quantity1,biggest_term
               PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
               STOP
            ENDIF
         ENDIF


      ELSE IF(z3(2).EQ.0.AND.z4(2).EQ.0.AND.z1(2).EQ.0.AND.z2(2).EQ.0)
     $     THEN

            PRINT *,"QUARTIC4: All imaginary in "
            PRINT *,"     ",u1term,u0term
            PRINT *,"     ",quantity1,biggest_term
            PRINT *,"     ",y1,ub,uc,ub1,ub2,uc1,uc2
            STOP

      ELSE                      !four solutions

      rtst=0
      write (*,*) 'Four solutions ',tsoln1,tsoln2,z1(1),z2(1),z3(1),
     &     z4(1)
      STOP
      IF(tsoln1.GT.0.0.AND.tsoln1.GE.tmin.AND.tsoln1.LE.tmax) THEN
         rtst=rtst+1
         soln=z1(1)
      END IF
         
      IF(tsoln2.GT.0.0.AND.tsoln2.GE.tmin.AND.tsoln2.LE.tmax) THEN
         rtst=rtst+1
         soln=z2(1)
      END IF

      IF(tsoln1.GT.0.0.AND.tsoln1.GE.tmin.AND.tsoln1.LE.tmax) THEN
         rtst=rtst+1
         soln=z3(1)
      END IF


      IF(tsoln1.GT.0.0.AND.tsoln1.GE.tmin.AND.tsoln1.LE.tmax) THEN
         rtst=rtst+1
         soln=z4(1)
      END IF

      IF(rtst.NE.1) THEN
         PRINT *,"QUARTIC4: There are four solutions and I'm
     $        incapable of"
         PRINT *,"picking one. Solns are: "
         PRINT *,z1(1),z2(1),z3(1),z4(1)
         PRINT *,tsoln1,tsoln2,tsoln3,tsoln4
         PRINT *,"Min..Max T :",tmin,tmax
         PRINT *,"I'm going back to trapimpl with moresweep2=.TRUE."
         moresweep=.TRUE.
!         STOP
      END IF
      END IF
c            if (HREAL.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'q 6'

!      PRINT *,"Solution chosen is: ",soln!,tmin,tmax
 541   CONTINUE

      END! SUBROUTINE QUARTIC4


