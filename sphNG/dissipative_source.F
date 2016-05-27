      SUBROUTINE dissipative_source(nlst_in,nlst_end,nlst_tot,listp,
     &     itime,npart,ntot,xyzmh,vxyzu,dha,trho,vsound,alphaMMpass,
     &     Bxyz,Bevolxyz,dBxyz,ekcle)
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

      DIMENSION xyzmh(5,mmax2),vxyzu(4,idim2)
      REAL*4 trho(idim2),vsound(idim2)
      REAL*4 dha(1+isizealphaMM,idim2),alphaMMpass(isizealphaMM,idim2)
      DIMENSION listp(idim2)
      DIMENSION Bxyz(3,imhd2)
      DIMENSION Bevolxyz(imhdevol,imhd3),dBxyz(imhdevol,imhd3)
      DIMENSION ekcle(5,iradtrans2)

      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/raddust'

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
C$OMP& shared(alphaMMpass,alphamin,alphamax)
C$OMP& shared(Bxyz,divcurlB,gradB,dBxyz,Bevolxyz,overclean_fac)
C$OMP& shared(ivar)
C$OMP& shared(dh2dt,dust_tk,vxyzu,ekcle,utime,udens,icall)
C$OMP& private(ipart,k)
C$OMP& private(hi,hi1,rhoi,rho1i,vsoundi)
C$OMP& private(fi,tdecay1,source)
C$OMP& private(Bxi,Byi,Bzi,B2i)
C$OMP& private(vsigi,vs2i,dB2,psi,gradBmag)
C$OMP& private(dust_temp,gas_temp,xnH2,ih2pos)
      DO n = nlst_in, nlst_end
         ipart = ivar(3,n)
         
         IF (ipart.EQ.0) THEN
            print*,' ERROR : ipart= ',ipart,'n = ',n,nlst_in
            CALL quit(1)
         ENDIF
c
c--Check for MPI code
c
         IF (ipart.GT.npart) THEN
            WRITE (iprint,*) 'Error: dissipative_source called ',
     &           'for non-local particle ',ipart,npart
            CALL quit(1)
         ENDIF
      
         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit(1)
         ENDIF
c
c--Zero dissipation for particle ipart [Note: dha(1,) is historical]
c
         DO k = 1, 1+isizealphaMM
            dha(k,ipart) = 0.0
         END DO
c
c--Dissipation and H_2 formation below only applies to gas particles
c
         IF (iphase(ipart).EQ.0) THEN

#ifdef RT
c
c--Evolution of H_2
c
            IF (idustRT.GT.0) THEN
c
c--Position value to allow full timestep and half timestep derivatives
c     to be used (as in normal R-K integrator)
c
               ih2pos = MIN(2,icall)
               dust_temp = dust_tk(1,ipart)
               gas_temp = vxyzu(4,ipart)/ekcle(3,ipart)
               xnH2 = trho(ipart)*udens/(gmw*mH)
               dh2dt(ih2pos,ipart) = h2_formation(ipart,gas_temp,
     &              dust_temp,xnH2) - h2_destruction(ipart,xnH2,.TRUE.)
c
c--Convert to code units
c
               dh2dt(ih2pos,ipart) = dh2dt(ih2pos,ipart)*
     &              utime*gmw*mH/(2.0*trho(ipart)*udens)
            ENDIF
#endif

c
c--Dissipative terms
c
            IF (ifsvi.GE.6) THEN
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
c                  dB2 = divcurlB(1,ipart)**2 + divcurlB(2,ipart)**2
c     &                 + divcurlB(3,ipart)**2 + divcurlB(4,ipart)**2
c                  IF (B2i.GT.tiny) THEN
c                     source = SQRT(MAX(dB2*rho1i,
c     &                    vs2i*divcurlB(1,ipart)**2/B2i))
c                  ELSE
c                     source = SQRT(dB2*rho1i)
c                  ENDIF
c                  dha(3,ipart) = (alphamin(2)-alphaMMpass(2,ipart))*
c     &                 tdecay1 + source

                  gradBmag = 0.
                  
c                  write (*,*) ipart,SQRT(xyzmh(1,ipart)**2 + 
c     &                 xyzmh(2,ipart)**2 + xyzmh(3,ipart)**2),
c     &                 gradB(9,ipart),gradB(7,ipart),gradB(8,ipart),
c     &                 gradB(6,ipart)


                  DO k = 1, 9
                     gradBmag = gradBmag + gradB(k,ipart)**2
                  END DO
                  gradBmag = SQRT(gradBmag)

                  IF (B2i.GT.tiny) THEN
                     dha(3,ipart) = MAX(MIN(hi * gradBmag /
     &                    SQRT(B2i),alphamax(2) ),alphamin(2) )
                  ELSE
                     dha(3,ipart) = 0.0
                  ENDIF

                  IF (imhdevol.EQ.4) THEN
                     psi = Bevolxyz(4,ipart)

                     vs2i = vs2i * overclean_fac**2
                     vsigi = vsigi * overclean_fac

c Use vsigi since we are _actually_ evolving psi/c_h NOT psi.
c Talk to Terry for why :-)
                     dBxyz(4,ipart) = - vsigi*divcurlB(1,ipart)
     &                    - 0.8*vsigi*hi1*psi
     &                    - 0.5*divv(ipart)*rho1i*psi
                  ENDIF
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
