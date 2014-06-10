      PROGRAM testcv

      INCLUDE 'COMMONS/optbl'
      INCLUDE 'COMMONS/eostbl'
      INCLUDE 'COMMONS/units'

      REAL*4 rho
      
      CHARACTER*10 filename
      CHARACTER*7 filebase

      OPEN(UNIT=8,FILE='/home/mbate/tables/opacitytbl')
      DO i=1, opmxtg
         READ(8,*) (optable(i,j), j=1, opmxrh)
      ENDDO
      CLOSE(8)
      
      OPEN(UNIT=8,FILE='/home/mbate/tables/gasttbl',
     &     FORM='unformatted')
      DO i=1, tgmxu
         READ(8) (tgtable(i,j), j=1, tgmxrh)
      ENDDO
      CLOSE(8)
 
      OPEN(UNIT=8,FILE='/home/mbate/tables/molmasstbl',
     &     FORM='unformatted')
      DO i=1, mumxu
         READ(8) (mutable(i,j), j=1, mumxrh)
      ENDDO
      CLOSE(8)

      gg = 6.67E-8
      umass = 2.0E+33
      udist = 1.0E+16
      utime = DSQRT(DBLE(udist)**3/(DBLE(gg)*DBLE(umass)))
      udens = umass/udist**3
      uergg = DBLE(udist)**2/DBLE(utime)**2
      filebase = 'radtest'
      
      DO irho = -16,9

         rho = 10.0**irho/udens

         IF (irho.GE.0) THEN
            WRITE (filename,88001) filebase,irho            
         ELSEIF (ABS(irho).LT.10) THEN
            WRITE (filename,88002) filebase,irho
         ELSE
            WRITE (filename,88003) filebase,irho            
         ENDIF
88001    FORMAT(A6,'00',I1)
88002    FORMAT(A6,I2)
88003    FORMAT(A6,I3)

         print *,'opening ',filename,rho*udens
         OPEN (11,FILE=filename)

         u = 0.025
         WRITE (11,*) '# ',rho*umass/udist/udist/udist

         DO i = 1, 100
            u = u*1.258
            cv = getcv(rho,u)
            rkappa = getkappa(u,cv,rho)
         WRITE(11,9901) u*uergg,cv*uergg,
     &           1.0/get1overmu(rho,u),
     &           rkappa*udist*udist/umass,u/cv
 9901       FORMAT(5(1PE12.5,1X))
         END DO

         print *,'closing ',filename
         CLOSE(11)

      END DO
      END
