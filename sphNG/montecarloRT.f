      SUBROUTINE montecarloRT(nlst_in,nlst_end,nlstall,list,npart,
     &         xyzmh,vxyzu,trho)
      
      INCLUDE 'idim'
      
      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      REAL*4 trho(idim)
      DIMENSION list(idim)

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'

      DIMENSION temperature(idim)
      
      CHARACTER*7 where

      DATA where/'monte'/
c     
c--Only perform radiative transfer if ptmass exist, otherwise don't modify u()
c     
c     
c--Find mass of gas particle (assume all same mass at moment)
c   
      IF (nptmass.GT.0) THEN
        WRITE (*,*) 'Monte Carlo, icall ',icall,nlstall,nlst_end
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               gaspartmass = xyzmh(4,i)
               GOTO 10
            ENDIF
         END DO
 10      CONTINUE
c     
c--Call torus for Monte-Carlo radiative transfer     
c
         CALL torus(idim,xyzmh,trho,iphase,nptmass,listpm,
     &       udist,umass,utime,time,gaspartmass,temperature)
c     
c--Set thermal energy of gas particles using torus's temperatures     
c     
         igas = 0
         iterate = 0
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               igas = igas + 1
               cvold = getcv(trho(i),vxyzu(4,i))
 20            vxyzu(4,i) = cvold*MAX(temperature(igas),10.0)
               cvnew = getcv(trho(i),vxyzu(4,i))
               IF (ABS((cvnew-cvold)/cvold).GT.0.001) THEN
                  cvold = cvnew
                  iterate = iterate + 1
                  IF (iterate.LT.100) THEN
                     GOTO 20
                  ELSE
                     WRITE (iprint,*) 'ERROR - cv iteration ',
     &                  cvold,cvnew,temperature(igas),igas
                  ENDIF
               ENDIF
            ENDIF
         END DO
      ENDIF 

      RETURN

      END
