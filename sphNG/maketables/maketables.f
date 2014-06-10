c-------------------------------------------------------------------------
c This is the driver program for creating the specific heat and mean
c molecular weight tables. Written by Stuart Whitehouse
c
c Was previously part of specificheat.f, program part now in separate 
c file so subroutines can be called from other programs (DJP 23/5/11)
c
c-------------------------------------------------------------------------
      PROGRAM MAKETABLE

      INCLUDE '../COMMONS/eostbl'

      REAL*8 ltm,lrho,tm,rho,specific,uoverT,cvv,muu,mu,Rg
      DIMENSION cvv(umxrh+1,umxt+1),muu(umxrh+1,umxt+1)
      INTEGER I,J

!		PRINT *,"Rho, Tg"
!		READ *,rho,tm
!		CALL GENERATECV(rho,tm,cv,mu)
!		PRINT *,"CV: ",cv,cv/8.3144d7
!		PRINT *,"Mu: ",mu
!		STOP

      Rg=8.3145d7

!     OPEN(UNIT=55,FILE='specheattbl')

      OPEN(UNIT=8,FILE='specheattbl',FORM='unformatted')
!      OPEN(UNIT=10,FILE='molmasstbl',FORM='unformatted')

      DO nrho=-30000,3000,5! log10(density) -20,3,0.005

         IF(MOD(nrho,500).EQ.0) PRINT *,nrho
         lrho=nrho/1000.0
         i=nrho/5+eostbl_rho1
         rho=10.0**lrho

         DO ntm=0,9000,5!0.05,4,0.005
            ltm=ntm/1000.0
            j=ntm/5+1
            tm=10.0**ltm
            CALL GENERATEU(rho,tm,specific,uoverT,mu)
c            write (95,99001) tm,specific/Rg,uoverT/Rg,mu
99001       FORMAT(4(1PE12.5,1X))
            cvv(I,J)=log10(uoverT)
            muu(I,J)=log10(mu)
         END DO

c         STOP

         IF(MOD(nrho,500).EQ.0) PRINT *,cvv(I,900)

      ENDDO
      PRINT *,"Doing making logcv table. Writing to disk..."
!J is temperature J = 1 to umxt rows
!I is density I = 1 to umxrh columns

      DO j=1,umxt
         WRITE(8) (cvv(i,j), i=1, umxrh)
!         WRITE(10) (muu(i,j), i=1, umxrh)
         PRINT *,cvv(2096,j),muu(2096,j)
c         WRITE(66,*) cvv(4600,j)
      END DO
!      PRINT *,cvv(20,50)
      PRINT *,"Complete"
      CLOSE(55)

      END 
            
c=========================================================================
