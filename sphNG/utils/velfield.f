      PROGRAM VELFIELD

C==============================================================
C  Generates a random turbulent velocity field which is
C  divergence-free.  Follows Dubinski, Narayan & Phillips 1985.
C  Based on original code for Zeldovich shift from Volker Bromm.
C  Velocity field code written by M. R. Bate (21/11/2000).
C
C==============================================================

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER(PI=3.1415927D0)
      PARAMETER(NGRID=16)

      REAL*8 MASSB,MASSDM,nindex,kmod,kdotq,LF,MASSTOT
      DIMENSION phix(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION phiy(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION phiz(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION pow(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION ampx(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION ampy(2*NGRID,2*NGRID,2*NGRID)
      DIMENSION ampz(2*NGRID,2*NGRID,2*NGRID)

      REAL*4 vel(2*NGRID,2*NGRID,2*NGRID,4)


      RTOT=1.0d0
      LF=2.0d0*RTOT
      AA=0.5d-1/32.52d0
      nindex=-6.0d0
c
c--Standard case:
c
      iseed=17
      iseed2=11
      print*,'Enter 2 seeds (integers, e.g. 17,11)'
      read (*,*) iseed
      read (*,*) iseed2
      print*,'Enter index (e.g. n=17/3 for Kolmogorov, 6 for Burgers)'
      read (*,*) nindex
c
c--Gives P(k)=k^{nindex) power spectrum
c     Note: rayldev returns the square root of the -Log of a random number
c
      powsum=0.d0
      do k=1,NGRID
        do j=1,2*NGRID
          do i=1,2*NGRID
            phix(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            phiy(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            phiz(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            kmod=sqrt(REAL((i-NGRID)**2+(j-NGRID)**2+k*k))
            pow(i,j,k)=AA*(kmod)**(nindex)
            powsum=powsum+pow(i,j,k)
            ampx(i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
            ampy(i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
            ampz(i,j,k)=rayldev(iseed2)*sqrt(pow(i,j,k))
          enddo
        enddo
      enddo

      sigma2=0.d0
      DO kk = 1, 2*NGRID
         print*,'kk ',kk
         DO jj = 1, 2*NGRID
         print*,'jj ',jj
            DO ii = 1, 2*NGRID
               xdmjk = (ii-NGRID-0.5)/REAL(NGRID)
               ydmjk = (jj-NGRID-0.5)/REAL(NGRID)
               zdmjk = (kk-NGRID-0.5)/REAL(NGRID)
               velx = 0.
               vely = 0.
               velz = 0.

C$OMP PARALLEL DO default(none)
C$OMP& shared(LF,xdmjk,ydmjk,zdmjk,ampx,ampy,ampz)
C$OMP& shared(phix,phiy,phiz)
C$OMP& private(i,j,k,kx,ky,kz,kdotq,contrib)
C$OMP& reduction(+:velx,vely,velz)
               do k=1,NGRID
                  do j=1,2*NGRID
                     do i=1,2*NGRID
                        kx = i-NGRID
                        ky = j-NGRID
                        kz = k
cc                        kmod=(2.d0*PI/LF)*
cc     &                       sqrt(REAL(kx**2+ky**2+kz**2))
                        kdotq=(2.d0*PI/LF)*(REAL(kx)*XDMJK+
     &                        REAL(ky)*YDMJK+REAL(kz)*ZDMJK)
                        contrib = ampz(i,j,k)*(2.d0*PI/LF)*REAL(ky)*
     &                       2.d0*SIN(kdotq+phiz(i,j,k))
     &                       - ampy(i,j,k)*(2.d0*PI/LF)*REAL(kz)*
     &                       2.d0*SIN(kdotq+phiy(i,j,k))
                        velx = velx + contrib
                        contrib = ampx(i,j,k)*(2.d0*PI/LF)*REAL(kz)*
     &                       2.d0*SIN(kdotq+phix(i,j,k))
     &                       - ampz(i,j,k)*(2.d0*PI/LF)*REAL(kx)*
     &                       2.d0*SIN(kdotq+phiz(i,j,k))
                        vely = vely + contrib
                        contrib = ampy(i,j,k)*(2.d0*PI/LF)*REAL(kx)*
     &                       2.d0*SIN(kdotq+phiy(i,j,k))
     &                       - ampx(i,j,k)*(2.d0*PI/LF)*REAL(ky)*
     &                       2.d0*SIN(kdotq+phix(i,j,k))
                        velz = velz + contrib
                     enddo
                  enddo
               enddo
C$OMP END PARALLEL DO

               vel(ii,jj,kk,1) = velx
               vel(ii,jj,kk,2) = vely
               vel(ii,jj,kk,3) = velz
            ENDDO
         ENDDO
      ENDDO

      DO kk = 1, 2*NGRID
         DO jj = 1, 2*NGRID
            DO ii = 1, 2*NGRID
               vel(ii,jj,kk,4) = SQRT(vel(ii,jj,kk,1)**2 + 
     &              vel(ii,jj,kk,2)**2 + vel(ii,jj,kk,3)**2)
            ENDDO
         ENDDO
      ENDDO

      OPEN (22, FILE='cube_v1.dat', FORM='unformatted')
      WRITE (22)(((vel(i,j,k,1),i=1,2*NGRID),
     &     j=1,2*NGRID),
     &     k=1,2*NGRID)
      CLOSE(22)
      OPEN (22, FILE='cube_v2.dat', FORM='unformatted')
      WRITE (22)(((vel(i,j,k,2),i=1,2*NGRID),
     &     j=1,2*NGRID),
     &     k=1,2*NGRID)
      CLOSE(22)
      OPEN (22, FILE='cube_v3.dat', FORM='unformatted')
      WRITE (22)(((vel(i,j,k,3),i=1,2*NGRID),
     &     j=1,2*NGRID),
     &     k=1,2*NGRID)
      CLOSE(22)
      OPEN (22, FILE='cube_v4.dat', FORM='unformatted')
      WRITE (22)(((vel(i,j,k,4),i=1,2*NGRID),
     &     j=1,2*NGRID),
     &     k=1,2*NGRID)
      CLOSE(22)
      STOP

      END

C*** The hopefully BETTER random generator ****
C** NOTE THAT this function is DOUBLE PRECISION!! 

C                                                                               
C
C     THE FUNCTION SUBROUTINE RAN3.
C                                                                               
C      THIS SUROUTINE GENERATES A UNIFORM RANDOM DEVIATE ON (0,1).              
C                                                                               
        DOUBLE PRECISION FUNCTION RAN3(IDUM)                                    
C
        IMPLICIT DOUBLE PRECISION (M)
          DOUBLE PRECISION FAC
         PARAMETER (MBIG=4000000., MSEED=1618033., MZ=0., FAC=1./MBIG)
          DIMENSION MA(55)
              SAVE
            DATA IFF/0/
          IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
             IFF=1
             MJ=MSEED-IABS(IDUM)
             MJ=MOD(MJ,MBIG)
             MA(55)=MJ
             MK=1
            DO 100 I=1,54
              II=MOD(21*I,55)
              MA(II)=MK
              MK=MJ-MK
              IF(MK.LT.MZ)MK=MK+MBIG
              MJ=MA(II)
100         CONTINUE
            DO 200 K=1,4
              DO 300 I=1,55
                MA(I)=MA(I)-MA(1+MOD(I+30,55))
                IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
300           CONTINUE
200         CONTINUE
             INEXT=0
             INEXTP=31
             IDUM=1
         ENDIF
            INEXT=INEXT+1
            IF(INEXT.EQ.56)INEXT=1
            INEXTP=INEXTP+1
            IF(INEXTP.EQ.56)INEXTP=1
            MJ=MA(INEXT)-MA(INEXTP)
            IF(MJ.LT.MZ)MJ=MJ+MBIG
            MA(INEXT)=MJ
            RAN3=MJ*FAC
             RETURN
             END

      REAL*8 FUNCTION rayldev(idum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CU    USES ran4
      REAL*8 ran4
      rayldev=dsqrt(-LOG(ran4(idum)))
      return
      END

C*** The hopefully BETTER random generator ****
C** NOTE THAT this function is DOUBLE PRECISION!! 

C                                                                               
C
C     THE FUNCTION SUBROUTINE RAN4.
C                                                                               
C      THIS SUROUTINE GENERATES A UNIFORM RANDOM DEVIATE ON (0,1).              
C                                                                               
        DOUBLE PRECISION FUNCTION RAN4(IDUM)                                    
C
        IMPLICIT DOUBLE PRECISION (M)
          DOUBLE PRECISION FAC
         PARAMETER (MBIG=4000000., MSEED=1618033., MZ=0., FAC=1./MBIG)
          DIMENSION MA(55)
              SAVE
            DATA IFF/0/
          IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
             IFF=1
             MJ=MSEED-IABS(IDUM)
             MJ=MOD(MJ,MBIG)
             MA(55)=MJ
             MK=1
            DO 100 I=1,54
              II=MOD(21*I,55)
              MA(II)=MK
              MK=MJ-MK
              IF(MK.LT.MZ)MK=MK+MBIG
              MJ=MA(II)
100         CONTINUE
            DO 200 K=1,4
              DO 300 I=1,55
                MA(I)=MA(I)-MA(1+MOD(I+30,55))
                IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
300           CONTINUE
200         CONTINUE
             INEXT=0
             INEXTP=31
             IDUM=1
         ENDIF
            INEXT=INEXT+1
            IF(INEXT.EQ.56)INEXT=1
            INEXTP=INEXTP+1
            IF(INEXTP.EQ.56)INEXTP=1
            MJ=MA(INEXT)-MA(INEXTP)
            IF(MJ.LT.MZ)MJ=MJ+MBIG
            MA(INEXT)=MJ
            RAN4=MJ*FAC
             RETURN
             END

