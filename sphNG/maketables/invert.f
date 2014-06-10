c
c Program to calculate T(u,rho) by inverting the table containing u(rho,T)
c
      PROGRAM INVERT

      IMPLICIT NONE

      INCLUDE '../COMMONS/eostbl'

      REAL tg1,rho,lrho,U,ltg1,ltg2,w,K,y1,y2,y3
      INTEGER I,J,nkrho1,N
      REAL cvtable(umxt,umxrh)

      OPEN(UNIT=8,FILE='specheattbl',FORM='unformatted')

      PRINT *,"Reading in CV table and making U table"
      DO I=1,umxt
        READ(8) (cvtable(i,j), j=1, umxrh)
        PRINT *,cvtable(i,2096),I
         tg1=10.0**((I-1)*0.005)
         DO J=1,umxrh
            utable(i,j)=(10.0**cvtable(i,j))*tg1
            IF(cvtable(i,j).EQ.0.0) THEN
               PRINT *,"Detected zero cv ",i,j,cvtable(i,j),tg1
               STOP
            ENDIF
         ENDDO
      END DO
      CLOSE(8)
      
      OPEN(UNIT=8, FILE='utbl', FORM='unformatted')
      DO i=1,umxrh
        WRITE(8) ((log10(utable(j,i))), j=1, umxt)
      END DO
      CLOSE(8)

      PRINT *,log10(utable(umxt,umxrh))
      PRINT *,log10(utable(umxt,1)),tg1
      PRINT *,"Done making U table"
      
      DO J=1,umxrh
         rho=10.0**((J-eostbl_rho1)*0.005)
         lrho=(J-eostbl_rho1)*0.005
          nkrho1=INT(lrho/0.005)+eostbl_rho1
          IF(lrho.LT.0.0) nkrho1=nkrho1-1
          nkrho1=J

         PRINT *,"Making density column ",nkrho1,lrho
         K=7.725
         DO I=1,tgmxu
            K=K+0.005
            u=10.0**REAL(K)

            DO N=1,tgmxu-1

               IF(utable(N,nkrho1).LE.u.AND.utable(N+1,nkrho1).GT.u)
     $              THEN

                  ltg1=((N)*0.005)
                  ltg2=((N+1)*0.005)
                  
                  w=(u-utable(N,nkrho1))/
     $                 (utable(N+1,nkrho1)-utable(N,nkrho1))
  
                  tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
c                  write(7,*) rho,u,tgtable(I,nkrho1)
                  GOTO 100
               ENDIF
            ENDDO ! N-loop Brute force loop to search for temperatures

!            PRINT *,"Can't seem to get u",I,u,utable(N+1,nkrho1),K
!            PRINT *,u,log10(u),801*0.005
!            tgtable(I,nkrho1)=99.00
!            STOP
            IF(u.GE.utable(tgmxu,nkrho1)) THEN
            !If failed, interpolate from last two points
            ltg1=((1000)*0.005)
            ltg2=1001*0.005
            w=(u-utable(tgmxu-1,nkrho1))/
     $           (utable(tgmxu,nkrho1)-utable(tgmxu-1,nkrho1))
            
            tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
            ELSE
               ltg1=((1)*0.005)
            ltg2=2*0.005
            w=(u-utable(1,nkrho1))/
     $           (utable(2,nkrho1)-utable(1,nkrho1))
            
            tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
               ENDIF
 100     ENDDO ! I loop ( U to tgmxu)
      ENDDO ! J loop (density to tgmxrh)


      OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
  
 66    FORMAT(7(1PE12.5,1x))
      DO I=1,tgmxu
        WRITE(8) ((tgtable(i,j)), j=1, tgmxrh)
  
      END DO
      CLOSE(8)

      PRINT *,(tgtable(50,150))
      PRINT *,(tgtable(50,50))

      OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
 
      DO I=1,tgmxu
         READ(8) (tgtable(i,j), j=1, tgmxrh)
      ENDDO
      
      CLOSE(8)
      PRINT *,tgtable(60,53)
      PRINT *,tgtable(94,356)

c      K=7.725
c      DO I=1,tgmxu
c         K=K+0.005
c         u=10.0**(REAL(K))
c         y1=tgtable(I,2000)
c         WRITE(57,*)  u,10**y1
c        y2=tgtable(I,3999)
c         WRITE(58,*)  u,10**y2
c        y3=tgtable(I,1)
c         WRITE(59,*)  u,10**y3
c      ENDDO
      END
