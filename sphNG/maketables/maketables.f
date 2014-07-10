c-------------------------------------------------------------------------
c
c--This is the driver program for creating the specific heat and mean
c     molecular weight tables. Originally written by Stuart Whitehouse
c
c--Was previously part of specificheat.f, program part now in separate 
c     file so subroutines can be called from other programs (DJP 23/5/11)
c
c--This program makes a table called "specheattbl", which gives the ratio
c     of the internal energy to the temperature of the gas, which would
c     be c_v = U/T if c_v was a constant.  However, the table is a function
c     of density and temperature.  The code sphNG does not evolve temperature,
c     it evolves specific internal energy, u.  So the code actually wants 
c     quantities as functions of  density and u.  These tables are made
c     by running "invert.f" after the "specheattbl" has been created.
c     This program inverts the table to give U/T as a function of density
c     and u, and produces a table of u  as a function of density and T.
c     It also creates a table of 1/mu, where mu is the mean particle mass 
c     which is also indexed using density and u.
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

      OPEN(UNIT=8,FILE='specheattbl',FORM='unformatted')

      DO nrho=-30000,3000,5! log10(density) -30,3,0.005

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
            muu(I,J)=log10(1./mu)
         END DO

         IF(MOD(nrho,500).EQ.0) PRINT *,cvv(I,900)

      ENDDO
      PRINT *,"Doing making logcv table. Writing to disk..."
!J is temperature J = 1 to umxt rows
!I is density I = 1 to umxrh columns

      DO j=1,umxt
         WRITE(8) (cvv(i,j), i=1, umxrh)
         PRINT *,cvv(2096,j),muu(2096,j)
      END DO
      CLOSE(8)

      PRINT *,"Complete"

      END 
            
c=========================================================================
