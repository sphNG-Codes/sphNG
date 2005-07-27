      PROGRAM temperature

      t1e20 = 10.

      d0 = 1.0e-20
      gam0 = 1.
      d1 = 1.0e-13
      d1 = 4.0e-13
      gam1 = 1.4
      gam1 = 5./3.
      d2 = d1 * (2000./t1e20)**(1.0/(gam1 - 1.0))
      gam2 = 1.1
      d3 = 1.0e-03
      gam3 = 5./3.

      dchange1 = 10**((LOG10(d0) + LOG10(d1))/2.)
      dchange2 = 10**((LOG10(d1) + LOG10(d2))/2.)
      dchange3 = 10**((LOG10(d2) + LOG10(d3))/2.)

c      WRITE (*,*) d0
c      WRITE (*,*) dchange1
c      WRITE (*,*) d1
c      WRITE (*,*) dchange2
c      WRITE (*,*) d2
c      WRITE (*,*) dchange3
c      WRITE (*,*) d3

      DO i = 0, 200
         density = 1.0e-20 * 10**(i/10.0)

c         IF (density.LT.dchange1) THEN
c            gamma = gam0
c            temperature = t1e20
c         ELSEIF (density.LT.dchange2) THEN
c         IF (density.LT.dchange2) THEN
            x = LOG10(density/d1)
            expx = EXP(x)
            fraction = expx/(expx + 1.0/expx)
            gamma = gam0 + (gam1 - gam0)*fraction
            t1 = t1e20
            t2 = t1e20*(density/d1)**(gam1-1.0)
            IF (density.LT.d1) t2 = (t2/t1e20)**0.7*t1e20
            temperature = t1 + t2
            temp2 = SQRT(t1**2 + t2**2)
            temp3 = (t1**0.9 + t2**0.9)**(1.0/0.9)
c            temp2 = t1e20*(density/d1)**(gam1-1.0)
c         ELSEIF (density.LT.dchange3) THEN
c            x = LOG10(density/d2)
c            expx = EXP(3.0*x)
c            fraction = expx/(expx + 1.0/expx)
c            gamma = gam1 + (gam2 - gam1)*fraction
c            temperature = t1e20*(dchange2/d1)**(gam1-1.0)
c     &           /(dchange2/d2)**(gamma-1.0)
c     &           *(density/d2)**(gamma-1.0)
c         ELSE
c            x = LOG10(density/d3)
c            expx = EXP(3.0*x)
c            fraction = expx/(expx + 1.0/expx)
c            gamma = gam2 + (gam3 - gam2)*fraction
c            temperature = t1e20*(gamma-1.0)
c         ENDIF

         WRITE (*,*) LOG10(density), gamma, LOG10(temperature), 
     &           LOG10(temp2), LOG10(temp3)

      END DO

      END

