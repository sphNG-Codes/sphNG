      SUBROUTINE planetesimals(np,ntot)
c************************************************************
c                                                           *
c  This subroutine positions planetesimals in orbits        *
c  defined using orbital elements: semi-major axis,         *
c  eccentricity, inclination, argument of periastron,       *
c  longitude of ascending node, and initial longitude of    *
c  planetesimal (equivalent to mean anomaly).               *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/xforce'

      REAL mina, maxa

      npart = np + ntot

      WRITE (*,1101)
 1101 FORMAT ('Please enter the maximum eccentricity: ')
      READ (*,*) eccentricity_p

      WRITE (*,1102)
 1102 FORMAT ('Please enter the maximum inclination: ')
      READ (*,*) inclination_p

      IF (ibound.EQ.102) THEN
         WRITE (*,1103)
 1103    FORMAT ('Enter the minimum and maximum semi-major axis:')
         READ (*,*) mina, maxa
         min_rplan = mina
         max_rplan = maxa
      ELSEIF (ibound.EQ.100) THEN
         mina = 1.0-variation
         maxa = 1.0+variation
      ENDIF

      DO i = ntot + 1, npart
         CALL orbital_elements (i, eccentricity_p, inclination_p, mina,
     &        maxa)
      END DO

      DO i = ntot + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END
      
      SUBROUTINE orbital_elements (i, eccentricity, inclination, mina,
     &     maxa)

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/gtime'

      REAL mina, maxa, inclination
      REAL e, in, p, asc, l, a, vdenom
      REAL e_anom, coeff, com, coeff_sum, nu
      INTEGER i, j, k, factorial

 100  e = ran1(1)*eccentricity  !eccentricity
      in = ran1(1)*inclination  !inclination
      p = ran1(1)*2.*pi         !argument of periastron
      asc = ran1(1)*2.*pi       !longitude of ascending node
      l = ran1(1)*2.*pi         !Called longitude in Wyatt 2003, but
                                !is equivalent to the mean anomaly
c         a = ran1(1)*(maxa-mina) + mina !semi-major axis
      IF (ibound.EQ.102) THEN
         IF (abs(sdprof+2.0).LT.tiny) THEN
            a = exp((log(maxa)-log(mina))*ran1(1) + log(mina))
         ELSE
            a = ((maxa**(2.+sdprof) - mina**(2.+sdprof))*
     &           ran1(1) + mina**(2.+sdprof))**(1./(2.+sdprof))
         ENDIF
      ELSEIF (ibound.EQ.100) THEN
         a = ((maxa**1.5 - mina**1.5)*ran1(1) + mina**1.5)
     &        **(2.0/3.0)
      ENDIF
c
c--Series expansion used to obtain e_anom from l = e-anom - e*sin(e_anom)
c
      e_anom = l
      DO j = 1, 20
         coeff = (1./(2.**(j-1)*factorial(j)))
         coeff_sum = 0.0
         DO k = 0, INT(j/2)
            com = factorial(j)/(factorial(k)*factorial(j-k))
            coeff_sum = coeff_sum + ((-1.)**k*com*(j-2.*k)**(j-1)*
     &           sin(l*REAL(j-1.*k)))
         ENDDO
         coeff = coeff*coeff_sum
         e_anom = e_anom + coeff*e**j
      ENDDO
      
c--Alternative method of calculating e_anom.
c         E0 = l
c 30      CONTINUE
c         Eorig = E0
c         E1 = E0 + (l + e*sin(E0) - E0)/(1. - e*cos(E0))
c         IF (abs(1. - abs(E1/Eorig)) .GT. 1.E-8) THEN
c            GOTO 30
c         ENDIF
c         e_anom = E1
      
      nu = 2.*ATAN(sqrt((1.+e)/(1.-e))*TAN(e_anom/2.))
      
 40   IF (nu.lt.0.0) THEN
         nu = nu + (2.*pi)            
         IF (nu.lt.0.0) THEN
            print *, 'nu < 0.0'
            goto 40
         ENDIF
      ENDIF
c
c--Shouldn't be used, but ensures the true anomaly (nu) is less than
c  pi if e_anom is less than pi, as required.
c
 50   IF (e_anom.LT.pi .AND. nu.GT.pi) THEN
         nu = nu - pi
         print *, i, nu
         goto 50
      ENDIF
      
c         re = ATAN((sqrt(1-e**2)*sin(nu))/(e+cos(nu)))
      
c
c--Convert orbital elements into cartesian coordinates
c
      radius = a*(1.-e*cos(e_anom))
      
c         write(26,1234) i, radius,e,in,p,asc,a,e_anom,nu
c 1234    FORMAT (I6,8(2X,1PE12.5))
      
      
      xyzmh(1,i) = radius*(cos(asc)*cos(p+nu) -
     &     sin(asc)*cos(in)*sin(p+nu))
      xyzmh(2,i) = radius*(sin(asc)*cos(p+nu) +
     &     cos(asc)*cos(in)*sin(p+nu))
      xyzmh(3,i) = radius*sin(in)*sin(p+nu)
c
c--Ensure no planetesimal is initially placed too near a planet,
c  or beyond relevant boundaries.
c
      IF (ibound.EQ.102 .AND. nptmass.GE.1) THEN
         DO j = 1, nptmass
            rtemp = sqrt((xyzmh(1,i)-xyzmh(1,listpm(j)))**2
     &           + (xyzmh(2,i)-xyzmh(2,listpm(j)))**2
     &           + (xyzmh(3,i)-xyzmh(3,listpm(j)))**2)
            IF (rtemp.LE.(xyzmh(5,listpm(j))*pradfac(j)*2.1))
     &           GOTO 100
         ENDDO
      ELSEIF (ibound.EQ.102 .AND. irotpot.EQ.1) THEN
         rtemp = sqrt((xyzmh(1,i)-rorbit_orig)**2 +
     &        xyzmh(2,i)**2 + xyzmh(3,i)**2)
         IF (rtemp.LE.(rplanet*pradfac(1)*2.1)) GOTO 100
      ELSEIF (ibound.EQ.100 .AND. gt.EQ.0.0) THEN
         IF ((radius.GT.(1.0+variation)) .OR. 
     &        ((radius.LT.(1.0-variation)))) GOTO 100
         IF (abs(atan2(xyzmh(2,i),xyzmh(1,i))).GT.phibound) goto 100
         rtemp = sqrt((xyzmh(1,i) - 1.0)**2 + xyzmh(2,i)**2 +
     &        xyzmh(3,i)**2)
         IF (rtemp.LE.(rplanet*pradfac(1)*2.1)) GOTO 100
      ELSEIF (ibound.EQ.100) THEN
         IF ((radius.GT.(1.0+variation)) .OR.
     &        ((radius.LT.(1.0-variation)))) GOTO 100         
      ENDIF
         
      xyzmh(5,i) = 0.1

      vdenom = a*(1.-e**2)
      vxyzu(1,i) = (xyzmh(1,i)*e/(sqrt(vdenom)*radius))*sin(nu) -
     &     (sqrt(vdenom)/radius)*(cos(asc)*sin(p+nu) +
     &     sin(asc)*cos(p+nu)*cos(in))
      vxyzu(2,i) = (xyzmh(2,i)*e/(sqrt(vdenom)*radius))*sin(nu) -
     &     (sqrt(vdenom)/radius)*(sin(asc)*sin(p+nu) -
     &     cos(asc)*cos(p+nu)*cos(in))
      vxyzu(3,i) = (xyzmh(3,i)*e/(sqrt(vdenom)*radius))*sin(nu) +
     &     (sqrt(vdenom)/radius)*sin(in)*cos(p+nu)

      RETURN
      END SUBROUTINE orbital_elements



      FUNCTION factorial(n)
      INTEGER i, n, factorial

      IF (n.EQ.0) THEN
         factorial = 1
      ELSE
         factorial = n
         DO i = n-1, 1, -1
            factorial = factorial*i
         ENDDO
      ENDIF
      
      END
