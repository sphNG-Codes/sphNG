
      PROGRAM INVERT

      IMPLICIT NONE

      REAL tg1,rho,lrho,U,ltg1,ltg2,w,K,y1,y2,y3
      INTEGER I,J,nkrho1,N
      REAL utable(1801,4601),tgtable(1799,4601),cvtable(1801,4601)
      
!      INCLUDE 'COMMONS/cvtbl'

      OPEN(UNIT=8,FILE='specheattbl',FORM='unformatted')
  
      PRINT *,"Reading in CV table and making U table"
      DO I=1,1801
        READ(8) (cvtable(i,j), j=1, 4601)
        PRINT *,cvtable(i,2096),I
         tg1=10.0**((I-1)*0.005)
         DO J=1,4601
            utable(i,j)=(10.0**cvtable(i,j))*tg1
				IF(cvtable(i,j).EQ.0.0) THEN
					PRINT *,"Detected zero cv"
					STOP
				ENDIF
         ENDDO
      END DO
      CLOSE(8)
      DO I=1,1801
        write(3,*) utable(i,1),I
                IF(MOD(I,100).EQ.0) THEN
          WRITE(94,66) utable(i,1),utable(i,100),utable(i,500),
     $     utable(i,1000),utable(i,2000),
     $     utable(i,3000),utable(i,4600)
        ENDIF
      ENDDO
      PRINT *,log10(utable(1801,4601))
      PRINT *,log10(utable(1801,1)),tg1
      PRINT *,"Done making U table"
      
      DO J=1,4601
         rho=10.0**((J-4001)*0.005)
         lrho=(J-4001)*0.005
          nkrho1=INT(lrho/0.005)+4001
          IF(lrho.LT.0.0) nkrho1=nkrho1-1
          nkrho1=J
!          PRINT *,lrho,nkrho1,J
!         IF(nkrho1.NE.J) THEN
!           PRINT *,"Failed"
!           STOP
!           ENDIF
         PRINT *,"Making density column ",nkrho1,lrho
         K=7.725
         DO I=1,1799
            K=K+0.005
            u=10.0**REAL(K)
 !           IF(nkrho1.EQ.1) PRINT *,u,utable(1,1),nkrho1,lrho,rho
 !           IF(nkrho1.EQ.2) STOP
            DO N=1,1798

               IF(utable(N,nkrho1).LE.u.AND.utable(N+1,nkrho1).GT.u)
     $              THEN
                  ltg1=((N)*0.005)
                  ltg2=((N+1)*0.005)
                  
                  w=(u-utable(N,nkrho1))/
     $                 (utable(N+1,nkrho1)-utable(N,nkrho1))
  
                  tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
                  write(7,*) rho,u,tgtable(I,nkrho1)

                  GOTO 100
               ENDIF
            ENDDO ! N-loop Brute force loop to search for temperatures

!            PRINT *,"Can't seem to get u",I,u,utable(N+1,nkrho1),K
!            PRINT *,u,log10(u),801*0.005
!            tgtable(I,nkrho1)=99.00
!            STOP
            IF(u.GE.utable(1799,nkrho1)) THEN
            !If failed, interpolate from last two points
            ltg1=((1000)*0.005)
            ltg2=1001*0.005
            w=(u-utable(1798,nkrho1))/
     $           (utable(1799,nkrho1)-utable(1798,nkrho1))
            
            tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
            ELSE
               ltg1=((1)*0.005)
            ltg2=2*0.005
            w=(u-utable(1,nkrho1))/
     $           (utable(2,nkrho1)-utable(1,nkrho1))
            
            tgtable(I,nkrho1)=((1.0-w)*ltg1+w*ltg2)
               ENDIF
 100     ENDDO ! I loop ( U to 1800)
      ENDDO ! J loop (density to 4601)


      OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
  
 66    FORMAT(7(1PE12.5,1x))
      DO I=1,1799
        WRITE(8) ((tgtable(i,j)), j=1, 4601)
  
      END DO
      CLOSE(8)

      PRINT *,(tgtable(50,150))
      PRINT *,(tgtable(50,50))


      OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
  

      DO I=1,1799
         READ(8) (tgtable(i,j), j=1, 4601)
      ENDDO
      
      CLOSE(8)
      PRINT *,tgtable(60,53)
      PRINT *,tgtable(94,356)
      K=7.725
      DO I=1,1799
         K=K+0.005
         u=10.0**(REAL(K))
         y1=tgtable(I,2000)
         WRITE(57,*)  u,10**y1
        y2=tgtable(I,3999)
         WRITE(58,*)  u,10**y2
        y3=tgtable(I,1)
         WRITE(59,*)  u,10**y3
      ENDDO
      END
