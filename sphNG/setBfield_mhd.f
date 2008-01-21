       SUBROUTINE setBfield
c************************************************************
c                                                           *
c  This subroutine sets up the initial magnetic field       *
c                                                           *
c  The distribution can be   uniform    (iBfield=1)         *
c                            uniform toroidal (iBfield=2)   *
c                            non-zero div B   (iBfield=3)   *      
c                            uniform poloidal (not yet)     *
c                                                           *
c Note that we always setup B even if evolving B/rho        *
c B/rho is then constructed once we know the density by sum *
c                                                           *
c************************************************************

      INCLUDE 'idim'
            
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/setBfield'
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/cylinder'
      
      INTEGER iBfield      
      CHARACTER*1 isetB
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine setBfield')
99002 FORMAT (A1)

      Bzero = 0.
      Bxzero = 0.
      Byzero = 0.
      Bzzero = 0.

50    WRITE(*, 99102)
99102 FORMAT (/,' Which MHD variable do you want to evolve?', /,
     &   ' B     (flux density per unit volume)         : Bvol b)', /,
     &   ' B/rho (flux density per unit mass)           : Brho r)', /,
     &   ' Euler potentials (B = grad alpha x grad beta): eulr e)', /,
     &   ' Vector potential (B = curl A)                : vecp v)')
      WRITE(*,*) ' (NB: a simulation can be restarted in B or B/rho'
      WRITE(*,*) '  from any run, but an Euler potential run can only '
      WRITE(*,*) '  restarted from an Euler potential run.)'
      WRITE(*,*)
      READ (*, 99024) varmhd
99024 FORMAT (A4)
      IF (varmhd.EQ.'b' .or. varmhd.EQ.'B') varmhd = 'Bvol'
      IF (varmhd.EQ.'r') varmhd = 'Brho'
      IF (varmhd.EQ.'e') varmhd = 'eulr'
      IF (varmhd.EQ.'v') varmhd = 'vecp'
      IF (varmhd.NE.'Bvol' .and. varmhd.NE.'Brho'
     &    .and. varmhd.NE.'eulr' .and. varmhd.NE.'vecp') GOTO 50
     
      WRITE(*,99025) ' Magnetic field variable = ',varmhd
99025 FORMAT(A,A4)

100   WRITE(*, 99003)
99003 FORMAT(//, ' INITIAL MAGNETIC FIELD GEOMETRY ', //,
     &        '       uniform cartesian field : 1 ',/,
     &        '        uniform toroidal field : 2 ',/,
     &        '          non-zero div B field : 3 ',/,
     &        '          tokamak Btheta field : 4 ',/,
     &        'force free z-periodic cylinder : 5')
      READ (*, *) iBfield
      IF (iBfield.LT.1.OR.iBfield.GT.5) GOTO 100
c
c--set initial mean pressure for use in beta calculations
c
      przero = 2./3.*thermal*rhozero
c
c--set critical mass-to-flux ratio in code units
c  (code units are G=1, mu_0=1 -- see Price & Monaghan 2004a for details)
c  c1 is a normalisation factor taken from Mouschovias & Spitzer 1976
c
      c1 = 0.53
      rmasstoflux_crit = 2./3.*c1*sqrt(5./pi)
c     set area for mass-to-flux calculation (assumed spherical at the moment)
      area = pi*rmax**2

c     field strength is fixed for torus and force free cylinder
      IF (iBfield.EQ.4.OR.iBfield.EQ.5) GOTO 300 

c
c--set magnitude of B field a variety of ways
c
200   WRITE(*,99004)
99004 FORMAT(' Do you want to enter: ',/,
     &       '         magnetic field strength in Gauss (m)',/,
     &       ' or magnetic field strength in code units (c)',/,
     &       '            or Alfven speed in code units (a)',/,
     &       '                      or mean plasma beta (b)',/, 
     &       '                or the mass-to-flux ratio (f)?')
      READ (*,99002) isetB
        
      IF (isetB.EQ.'m') THEN
c
c--magnetic field strength in Gauss
c
         WRITE (*,99104) umagfd
99104    FORMAT(' Code units of mag flux density = ',1pe14.6)

         WRITE (*,99005)
99005    FORMAT(/,' Enter Mag flux density in Gauss')
         READ (*,*) Bzero
         Bzero = Bzero/umagfd         
      
      ELSEIF (isetB.EQ.'c') THEN
c
c--magnetic field strength in code units
c
         WRITE (*,99006) umagfd
99006    FORMAT(/,' Enter Mag flux density in units of ',1pe14.6,' G')
         READ (*,*) Bzero

      ELSEIF (isetB.eq.'a') THEN
c
c--Alfven speed in code units
c
         WRITE (*,99007)
99007    FORMAT (' Enter Alfven speed in code units')
         READ (*,*) valfven
         Bzero = SQRT(rhozero)*valfven

      ELSEIF (isetB.eq.'b') THEN
c
c--mean plasma beta
c
         WRITE (*,98007)
98007    FORMAT (' Enter mean plasma beta (gas/mag pressure)',/,
     &           ' (assuming uniform density and temperature, 0=hydro)')
         READ (*,*) betazero
         IF (betazero.LT.0.) THEN
            STOP 'beta must be > 0'
         ELSEIF (betazero.EQ.0) THEN
            Bzero = 0.
         ELSE
            Bzero = SQRT(2.*przero/betazero)
         ENDIF

      ELSEIF (isetB.eq.'f') THEN
c
c--mass-to-flux ratio for a spherical cloud
c
         WRITE(*,*) 'critical mass-to-flux ratio (code units) = ',
     &            rmasstoflux_crit
         WRITE (*,98107)
98107    FORMAT (' Enter spherical mass-to-flux ratio',/,
     &           ' (where < 1 = subcritical,',/,
     &           '          1 = critical,',/,
     &           '        > 1 = supercritical, ',/,
     &           '          0=inf=hydro )')
         READ (*,*) rmasstoflux
         IF (rmasstoflux.LT.0.) THEN
            STOP 'ratio must be >= 0'
         ELSEIF (rmasstoflux.LT.tiny) THEN
            Bzero = 0.
         ELSE
            Bzero = totmas/(area*rmasstoflux*rmasstoflux_crit)
         ENDIF

      ELSE
         GOTO 200
      ENDIF

300   CONTINUE

      IF (iBfield.EQ.1) THEN
c
c--uniform cartesian field
c
         WRITE (*,99008)
99008    FORMAT (' Enter Bx:By:Bz ratio')
         READ (*,*) fracx,fracy,fracz
         fractot = SQRT(fracx**2. + fracy**2. + fracz**2.)
         Bxzero = fracx*Bzero/fractot
         Byzero = fracy*Bzero/fractot
         Bzzero = fracz*Bzero/fractot
c
c--spit out actual settings
c
         WRITE (*,99009) Bxzero,Byzero,Bzzero
99009    FORMAT (' Bx_0 = ',1PE14.5,/,
     &           ' By_0 = ',1PE14.5,/,
     &           ' Bz_0 = ',1PE14.5)

         IF (varmhd.EQ.'eulr') THEN
            IF (abs(fracz-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = -Bzero*xyzmh(2,i)
                  Bevolxyz(2,i) = xyzmh(1,i)
                  Bevolxyz(3,i) = 0.
               ENDDO
            ELSEIF (abs(fracy-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = -Bzero*xyzmh(1,i)
                  Bevolxyz(2,i) = xyzmh(3,i)
                  Bevolxyz(3,i) = 0.
               ENDDO
            ELSEIF (abs(fracx-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = -Bzero*xyzmh(3,i)
                  Bevolxyz(2,i) = xyzmh(2,i)
                  Bevolxyz(3,i) = 0.
               ENDDO
            ELSE
             STOP 'mixed cartesian field NOT IMPLEMENTED for EULER POTS'
            ENDIF
         ELSEIF (varmhd.EQ.'vecp') THEN
c
c--for vector potential with cartesian fields
c  these should be entirely set as EXTERNAL fields (done below anyway)
            IF (.FALSE.) THEN
            
            WRITE(*,*) 'Setting vector potential to zero '
            WRITE(*,*) '=> using external B fields'
            DO i=1,npart
               Bevolxyz(1,i) = 0.
               Bevolxyz(2,i) = 0.
               Bevolxyz(3,i) = 0.
            ENDDO
            
            ELSE
c
c  the below is an alternative which should work but is a bit weird 
c  and is only implemented for the purposes of testing
c  NB: there is a Gauge choice doing it this way
c
            IF (abs(fracz-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = -Bzero*xyzmh(2,i)
                  Bevolxyz(2,i) = 0.
                  Bevolxyz(3,i) = 0.
                  ! alternative
c                  Bevolxyz(1,i) = 0.
c                  Bevolxyz(2,i) = Bzero*xyzmh(1,i)
c                  Bevolxyz(3,i) = 0.
               ENDDO
            ELSEIF (abs(fracy-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = Bzero*xyzmh(3,i)
                  Bevolxyz(2,i) = 0.
                  Bevolxyz(3,i) = 0.
c                  ! alternative
c                  Bevolxyz(1,i) = 0.
c                  Bevolxyz(2,i) = 0.
c                  Bevolxyz(3,i) = -Bzero*xyzmh(1,i)
               ENDDO
            ELSEIF (abs(fracx-1.0).LT.tiny) THEN
	       DO i=1,npart
                  Bevolxyz(1,i) = 0.
                  Bevolxyz(2,i) = -Bzero*xyzmh(3,i)
                  Bevolxyz(3,i) = 0.
                  ! alternative
c                  Bevolxyz(1,i) = 0.
c                  Bevolxyz(2,i) = 0.
c                  Bevolxyz(3,i) = Bzero*xyzmh(2,i)
               ENDDO
            ELSE
             STOP 'mixed cartesian field NOT IMPLEMENTED for VECTOR POT'
            ENDIF
            
            ENDIF       
         ELSE
            DO i=1,npart
               Bevolxyz(1,i) = Bxzero
               Bevolxyz(2,i) = Byzero
               Bevolxyz(3,i) = Bzzero
            ENDDO        
         ENDIF
      
      ELSEIF (iBfield.EQ.2) THEN
c
c--uniform toroidal field
c
         IF (varmhd.EQ.'eulr') THEN
            DO i=1,npart
               rsph2 = xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2
               rsph = sqrt(rsph2)
               theta = acos(xyzmh(3,i)/rsph)
               Bevolxyz(1,i)= -Bzero*theta
               Bevolxyz(2,i)= 0.5*(rsph2)
               Bevolxyz(3,i)= 0.
            ENDDO
         ELSE
            STOP 'not implemented for non-Euler potentials'
         ENDIF
         
      ELSEIF (iBfield.EQ.3) THEN ! this is NOT divergence free
c
c  Setup for the Bx peak advection problem in Dedner et al JCP 175, 645
c  Bx = r(x^2 + y^2)/sqrt(4pi) (ie div B .ne. 0) 
c  Basically to see how an initially non-zero div B propagates
c
         IF (varmhd(1:1).NE.'B') THEN
           STOP 'not implemented for eulr or vecp'
         ENDIF
         WRITE(*,99104) umagfd
         WRITE(*,99106) 'Enter Binit'
99106    FORMAT(A)
         READ(*,*) Bzero
         rbump = 4./sqrt(2.)        ! radius of the initial bump
         rbump2 = rbump*rbump

         Bxzero = 0.
         Byzero = 0.
         Bzzero = Bzero/sqrt(4.*pi)
         DO i=1,npart
            rr = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i) +
     &          xyzmh(3,i)*xyzmh(3,i)
            IF (rr.le.rbump2) THEN
               Bevolxyz(1,i) = ((rr/rbump2)**4 - 2.*(rr/rbump2)**2 
     &                            + 1.)/sqrt(4.*pi)*Bzero
            ELSE
               Bevolxyz(1,i) = 0.
            ENDIF
            Bevolxyz(2,i) = 0.
            Bevolxyz(3,i) = Bzzero
         ENDDO
      ELSEIF (iBfield.EQ.4) THEN
c
c  Setup for the tokamak confining field (Btheta)
c  corresponding to the external force iexf=9
c  taken from Wesson, J., Tokamaks, 3rd Ed, Oxford 2004, page 123
c
         WRITE(*,*) 'SETUP FOR TOKAMAK BTHETA FIELD'
         WRITE(*,*) ' J0 = ',currj0,' a = ',atorus,'nu = ',nutorus

         DO i=1,npart
c
c        get position in torus r coordinate
c
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
            rcyl = sqrt(xi**2 + yi**2)
            rintorus = sqrt((rcyl - Rtorus)**2 + zi**2)

            IF (rintorus.GT.tiny) THEN
c
c        set Btheta in torus co-ordinates
c
               ra2 = rintorus**2*da2
               Btheta = (0.5*currj0*atorus**2)/(nutorus+1)*
     &                  (1.- (1.-ra2)**(nutorus+1))/rintorus

               sintheta = zi/rintorus
               theta = ATAN2(zi,rcyl-Rtorus)
               phi = ATAN2(yi,xi)

               IF (varmhd.EQ.'eulr') THEN
c
c        set Euler potentials directly in torus co-ordinates
c
                  IF (nutorus.EQ.2) THEN
                     Bevolxyz(1,i) = currj0*atorus**2*rintorus*
     &               ((ra2**3)/6. - 0.75*ra2**2 + 1.5*ra2)*phi/6.
                  ELSE 
                     STOP 'not implemented for nutorus.ne.2'
                  ENDIF
                  Bevolxyz(2,i) = 0.5*rintorus**2
                  Bevolxyz(3,i) = 0.
               ELSEIF (varmhd(1:1).EQ.'B') THEN
c
c        otherwise transform to get Bx, By and Bz
c            
                  Bevolxyz(1,i) = -Btheta*sintheta*COS(phi)
                  Bevolxyz(2,i) = -Btheta*sintheta*SIN(phi)
                  Bevolxyz(3,i) = Btheta*COS(theta)
               ELSE
                  STOP 'not implemented for vecp'
               ENDIF
               Bzero = MAX(Bzero,Btheta)
            ELSE
               Bevolxyz(1,i) = 0.
               Bevolxyz(2,i) = 0.
               Bevolxyz(3,i) = 0.
            ENDIF
         ENDDO
      ELSEIF (iBfield.EQ.5) THEN
c
c  Setup for the Force Free Cylinder (Btheta)
c  corresponding to the external force iexf=10
c  taken from 
c
         WRITE(*,*) 'SETUP FOR THE FORCE FREE Z-PERIODIC CYLINDER'
         WRITE(*,*) ' Radius ',radius,' Length = ',length
         WRITE (*,*)' Muff = ',Muff, 'Amplitude =', ampl
         IF (varmhd(1:1).NE.'B') THEN
            WRITE (*,*) 'Not implemented for Euler potentials'
            STOP
         ElSE   
            DO i=1,npart
c
c        get position in cylindrical r coordinate
c
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               rr= sqrt(xi**2 + yi**2)
               IF (rr.GT.tiny) THEN
                  drr = 1./rr
               ELSE 
                  drr = 0.
               ENDIF
            STOP 'bessel functions commented out in repository version'
               !Bevolxyz(1,i) = -ampl*yi*drr*besj1(Muff*rr)
               !Bevolxyz(2,i) = ampl*xi*drr*besj1(Muff*rr)
               !Bevolxyz(3,i) = ampl*besj0(Muff*rr)
            ENDDO   
         ENDIF  
         Bzero = ampl
      ENDIF
c
c--spit out various information about the magnetic field we have set up
c
      Bzero2 = Bzero**2
      valfven = SQRT(Bzero2/rhozero)
      betazero = przero/(0.5*Bzero2)
      WRITE(*,98009) valfven,betazero
98009 FORMAT (' Alfven speed = ',1pe10.4,/,' Plasma beta  = ',1pe10.4)

c
c--spit out flux to mass ratio (assumes spherical geometry at the moment)
c
        IF (Bxzero.GT.tiny) WRITE(*,98010) 'x',
     &      totmas/(area*Bxzero),
     &      totmas/(area*Bxzero)/rmasstoflux_crit
        IF (Byzero.GT.tiny) WRITE(*,98010) 'y',
     &      totmas/(area*Byzero),
     &      totmas/(area*Byzero)/rmasstoflux_crit
        IF (Bzzero.GT.tiny) WRITE(*,98010) 'z',
     &      totmas/(area*Bzzero),
     &      totmas/(area*Bzzero)/rmasstoflux_crit

98010 FORMAT (' Mass to flux ratio (',a1,') = ',es10.4,
     &        ' r/rcrit = ',f9.4)
c
c--calculate angle of uniform field to x axis
c        
      angle = ACOS(Bxzero/Bzero)*180.0/pi
      WRITE (*,99010) angle
99010 FORMAT (' Angle between field and x axis = ',f7.3,' degrees')
c
c--spit out field strength in CGS units
c        
      WRITE(*,99011) Bzero*umagfd,Bzero
99011 FORMAT(' Initial field strength (Bzero) = ',1pe10.4,' G',
     &       ' (',1pe10.4,' in code units)')

c
c--set external field components
c (required for const. stress boundaries and current projection)
c
      Bextx = Bxzero
      Bexty = Byzero
      Bextz = Bzzero
      WRITE(*,99055) Bextx,Bexty,Bextz
99055 FORMAT(' Setting external B field = ',3(1PE12.4,2X),//) 

      RETURN
      END
