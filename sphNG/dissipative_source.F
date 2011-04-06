      SUBROUTINE dissipative_source(nlst_in,nlst_end,nlst_tot,listp,
     &     itime,npart,ntot,xyzmh,dha,trho,vsound,alphaMMpass,Bxyz)
c************************************************************
c                                                           *
c  This subroutine computes the dissipative source terms    *
c     for artificial viscosity and resistivity.             *
c     These used to be calculated in forcei before MPI.     *
c     MRB 4/4/2011                                          *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2)
      REAL*4 trho(idim2),vsound(idim2)
      REAL*4 dha(1+isizealphaMM,idim2),alphaMMpass(isizealphaMM,idim2)
      DIMENSION listp(idim2)
      DIMENSION Bxyz(3,imhd2)

      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/compact'

c      DATA epsil2/1.E-4/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine dissipative_source')
c
c--Calculate dissipative source terms
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,npart,listp)
C$OMP& shared(xyzmh,dha,trho,vsound,divv,curlv)
C$OMP& shared(iphase,iprint,ifsvi)
C$OMP& shared(alphaMMpass,alphamin)
C$OMP& shared(Bxyz,divcurlB)
C$OMP& shared(ivar)
C$OMP& private(ipart,k)
C$OMP& private(hi,hi1,rhoi,rho1i,vsoundi)
C$OMP& private(fi,tdecay1,source)
C$OMP& private(Bxi,Byi,Bzi,B2i)
C$OMP& private(vsigi,vs2i,dB2)
      DO n = nlst_in, nlst_end
         ipart = ivar(3,n)
         
         IF (ipart.EQ.0) THEN
            print*,' ERROR : ipart= ',ipart,'n = ',n,nlst_in
            CALL quit
         ENDIF
c
c--Check for MPI code
c
         IF (ipart.GT.npart) THEN
            WRITE (iprint,*) 'Error: dissipative_source called ',
     &           'for non-local particle ',ipart,npart
            CALL quit
         ENDIF
      
         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ENDIF
c
c--Zero dissipation for particle ipart [Note: dha(1,) is historical]
c
         DO k = 1, 1+isizealphaMM
            dha(k,ipart) = 0.0
         END DO
c
c--Dissipation below only applies to gas particles
c
         IF (iphase(ipart).EQ.0) THEN
            IF (ifsvi.EQ.6) THEN
               hi = xyzmh(5,ipart)
               hi1 = 1./hi
               rhoi = trho(ipart)
               rho1i = 1./rhoi

               IF (imhd.EQ.idim .AND. iphase(ipart).EQ.0) THEN
                  Bxi = Bxyz(1,ipart)
                  Byi = Bxyz(2,ipart)
                  Bzi = Bxyz(3,ipart)
                  B2i = Bxi**2 + Byi**2 + Bzi**2
               ELSE
                  B2i = 0.
               ENDIF

               vsoundi = vsound(ipart)
               vs2i = vsoundi**2 + B2i*rho1i
c
c--Artificial viscosity and resistivity source terms
c
c--Morris & Monaghan switch source and decay terms for both artificial 
c     viscosity and artificial resistivity (see Price & Monaghan 2005)
c
               vsigi = SQRT(vs2i)
               tdecay1 = 0.2*vsigi*hi1
c
c--Balsara factor
c               adivi = ABS(divv(ipart)*rho1i)
c               acurlvi = ABS(curlv(ipart)*rho1i)
c               fi = adivi/(adivi+acurlvi+epsil2*vsigi*hi1)

cc         dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
cc     &           MIN(divv(ipart)*rho1i+0.5*vsoundi*hi1,0.0)
               dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*
     &              tdecay1 - MIN(divv(ipart)*rho1i,0.0) !!*fi
               IF (imhd.EQ.idim) THEN
                  dB2 = divcurlB(1,ipart)**2 + divcurlB(2,ipart)**2
     &                 + divcurlB(3,ipart)**2 + divcurlB(4,ipart)**2
                  IF (B2i.GT.tiny) THEN
                     source = SQRT(MAX(dB2*rho1i,
     &                    vs2i*divcurlB(1,ipart)**2/B2i))
                  ELSE
                     source = SQRT(dB2*rho1i)
                  ENDIF
                  dha(3,ipart) = (alphamin(2)-alphaMMpass(2,ipart))*
     &                 tdecay1 + source
               ENDIF
            ENDIF
         ENDIF

      END DO
C$OMP END PARALLEL DO
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine dissipative_source')

      RETURN
      END
