       SUBROUTINE setBfield
c************************************************************
c                                                           *
c  This subroutine sets up the initial magnetic field       *
c                                                           *
c  The distribution can be   uniform    (iBfield=1)         *
c			     random     (iBfield=2)	    *      
c                            uniform toroidal (not yet)     *
c                            uniform poloidal (not yet)     *
c							    *
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
      INCLUDE 'COMMONS/Bzero'
      
      INTEGER iBfield      
      CHARACTER*1 isetB
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine setBfield')
99002 FORMAT (A1)

100   WRITE(*, 99003)
99003 FORMAT(//, ' INITIAL MAGNETIC FIELD GEOMETRY ', //,
     &        '                 uniform field : 1 ',/,
     &        '                        random : 2 ')
      READ (*, *) iBfield
      IF (iBfield.LT.1.OR.iBfield.GT.2) GOTO 100


      IF(iBfield.EQ.1) THEN

200     WRITE(*,99004)
99004   FORMAT(' Do you want to enter the magnetic field strength (m)',
     &         ' or the Alfven speed (a) ?')	
        READ (*,99002) isetB
	IF (isetB.NE.'m'.AND.isetB.NE.'a') GOTO 200
	
        IF (isetB.EQ.'m') THEN
           WRITE (*,99104) umagfd
99104      FORMAT(' Code units of mag flux density = ',1pe14.6)

           WRITE (*,99005)
99005      FORMAT(/,' Enter Mag flux density in Gauss')     
	   WRITE (*,99006) 'Bx'
	   READ (*,*) Bxzero
	   WRITE (*,99006) 'By'
	   READ (*,*) Byzero
	   WRITE (*,99006) 'Bz'
	   READ (*,*) Bzzero
99006      FORMAT(' Enter ',A2,':')	
           Bxzero = Bxzero/umagfd
           Byzero = Byzero/umagfd
           Bzzero = Bzzero/umagfd
	ELSEIF (isetB.eq.'a') THEN
	   WRITE (*,99007)
99007      FORMAT (' Enter Alfven speed in code units')
	   READ (*,*) valfven
	   WRITE (*,99008)
99008      FORMAT (' Enter Bx:By:Bz ratio')
	   READ (*,*) fracx,fracy,fracz
	   fractot = SQRT(fracx**2. + fracy**2. + fracz**2.)
	   Bxzero = fracx*SQRT(rhozero)*valfven/fractot
	   Byzero = fracy*SQRT(rhozero)*valfven/fractot
	   Bzzero = fracz*SQRT(rhozero)*valfven/fractot
	   
	ENDIF
        

	WRITE (*,99009) Bxzero,Byzero,Bzzero
99009   FORMAT (' Bx_0 = ',1PE14.5,/,
     &          ' By_0 = ',1PE14.5,/,
     &          ' Bz_0 = ',1PE14.5)	   

	Bzero2 = Bxzero**2 + Byzero**2 + Bzzero**2
        Bzero = SQRT(Bzero2)
c
c--calculate angle of uniform field to x axis
c	
        angle = ACOS(Bxzero/Bzero)*180.0/pi
        WRITE (*,99010) angle
99010   FORMAT (' Angle between field and x axis = ',f7.3,' degrees')
		
	DO i=1,npart
	   Bevolxyz(1,i) = Bxzero
           Bevolxyz(2,i) = Byzero
	   Bevolxyz(3,i) = Bzzero
	ENDDO

      ELSEIF(iBfield.EQ.2) THEN		! this is NOT divergence free
        WRITE(*,99005) umagfd
        WRITE(*,99006) 'maximum B'
        READ(*,*) Binit	
	Binit5 = 0.5*Binit
	DO i=1,npart
	   Bevolxyz(1,i) = 2.0*(Binit*ran1(1) - Binit5)
	   Bevolxyz(2,i) = 2.0*(Binit*ran1(1) - Binit5)
	   Bevolxyz(3,i) = 2.0*(Binit*ran1(1) - Binit5)
	ENDDO
      ENDIF          
      
      RETURN
      END
