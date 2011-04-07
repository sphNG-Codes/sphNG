      SUBROUTINE planetesimals(np)
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

      REAL eccentricity, inclination, periastron, ascending, longitude
      REAL mina, maxa
      REAL e, in, p, asc, l, a, vdenom
      REAL e_anom, coeff, com, coeff_sum, nu
      INTEGER i, j, k, factorial

      npart = np + nptmass

      WRITE (*,1101)
 1101 FORMAT ('Please enter the maximum eccentricity: ')
      READ (*,*) eccentricity

      WRITE (*,1102)
 1102 FORMAT ('Please enter the maximum inclination: ')
      READ (*,*) inclination

      WRITE (*,1103)
 1103 FORMAT ('Please enter the minimum and maximum semi-major axis:')
      READ (*,*) mina, maxa

      DO i = nptmass + 1, npart
 100     e = ran1(1)*eccentricity  !eccentricity
         in = ran1(1)*inclination  !inclination
         p = ran1(1)*2.*pi         !argument of periastron
         asc = ran1(1)*2.*pi       !longitude of ascending node
         l = ran1(1)*2.*pi         !Called longitude in Wyatt 2003, but
                                   !is equivalent to the mean anomaly
         a = ran1(1)*(maxa-mina) + mina !semi-major axis
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
     &              sin(l*REAL(j-1.*k)))
            ENDDO
            coeff = coeff*coeff_sum
            e_anom = e_anom + coeff*e**j
         ENDDO
         nu = 2.*ATAN(sqrt((1.+e)/(1.-e))*TAN(e_anom/2.))

 40      IF (nu.lt.0.0) THEN
            nu = nu + (2.*pi)            
            IF (nu.lt.0.0) goto 40
         ENDIF
c
c--Shouldn't be used, but ensures the true anomaly (nu) is less than
c  pi if e_anom is less than pi, as required.
c
 50      IF (e_anom.LT.pi .AND. nu.GT.pi) THEN
            nu = nu - pi
            print *, i, nu
            goto 50
         ENDIF

c         re = ATAN((sqrt(1-e**2)*sin(nu))/(e+cos(nu)))

c
c--Convert orbital elements into cartesian coordinates
c
         radius = a*(1-e*cos(e_anom))

c         write(26,1234) i, radius,e,in,p,asc,a,e_anom,nu
c 1234    FORMAT (I6,8(2X,1PE12.5))


         xyzmh(1,i) = radius*(cos(asc)*cos(p+nu) -
     &        sin(asc)*cos(in)*sin(p+nu))
         xyzmh(2,i) = radius*(sin(asc)*cos(p+nu) +
     &        cos(asc)*cos(in)*sin(p+nu))
         xyzmh(3,i) = radius*sin(in)*sin(p+nu)

c
c--Ensure no planetesimal is initially placed too near a planet.
c
         IF (nptmass.GE.1) THEN
            DO j = 1, nptmass
               IF (sqrt((xyzmh(1,i)-xyzmh(1,listpm(j)))**2 +
     &              (xyzmh(2,i)-xyzmh(2,listpm(j)))**2
     &              + (xyzmh(3,i)-xyzmh(3,listpm(j)))**2) .LE.
     &              (xyzmh(5,listpm(j))*pradfac(j)*2.1))
     &              GOTO 100
            ENDDO
         ENDIF
         
         xyzmh(5,i) = 0.1

         vdenom = a*(1.-e**2)
         vxyzu(1,i) = (xyzmh(1,i)*e/(sqrt(vdenom)*radius))*sin(nu) -
     &        (sqrt(vdenom)/radius)*(cos(asc)*sin(p+nu) +
     &        sin(asc)*cos(p+nu)*cos(in))
         vxyzu(2,i) = (xyzmh(2,i)*e/(sqrt(vdenom)*radius))*sin(nu) -
     &        (sqrt(vdenom)/radius)*(sin(asc)*sin(p+nu) -
     &        cos(asc)*cos(p+nu)*cos(in))
         vxyzu(3,i) = (xyzmh(3,i)*e/(sqrt(vdenom)*radius))*sin(nu) +
     &        (sqrt(vdenom)/radius)*sin(in)*cos(p+nu)

            amx = xyzmh(2,i)*vxyzu(3,i) - xyzmh(3,i)*vxyzu(2,i)
            amy = xyzmh(3,i)*vxyzu(1,i) - xyzmh(1,i)*vxyzu(3,i)
            amz = xyzmh(1,i)*vxyzu(2,i) - xyzmh(2,i)*vxyzu(1,i)

            angmom = amx + amy + amz

c         write(27,1234) i, angmom,e,in,p,asc,a,e_anom,nu
c 1234    FORMAT (I6,8(2X,1PE12.5))

      END DO

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END

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
