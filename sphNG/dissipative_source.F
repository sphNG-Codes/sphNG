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
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/HY09accel'
      INCLUDE 'COMMONS/tstopvar'

c      DATA epsil2/1.E-4/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine dissipative_source')
c
c--Position value to allow full timestep and half timestep derivatives
c     to be used (as in normal R-K integrator)
c
      ih2pos = MIN(2,icall)
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
C$OMP& shared(dHIIdt,umass,nptmasstot,stellarrad,ih2pos)
C$OMP& shared(stellar_xyzmhrti,HIIion)
C$OMP& shared(dtmax,imaxstep,isteps)
C$OMP& shared(grav_accel,omega_disc,iexf,xmass,udist)
C$OMP& private(ipart,k,dti)
C$OMP& private(hi,hi1,rhoi,rho1i,vsoundi)
C$OMP& private(fi,tdecay1,source)
C$OMP& private(Bxi,Byi,Bzi,B2i)
C$OMP& private(vsigi,vs2i,dB2,psi,gradBmag)
C$OMP& private(dust_temp,gas_temp,xnH2,isink,h2_over_dist2,recomb)
C$OMP& private(rgrav_accel_code_units,radius,vel_dot_grav_accel)
C$OMP& private(vel2,vel_tangential,ratio_vrad_to_vtan,ratio_vrad_to_cs)
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
c--Dissipation, H_2 formation, and stellar feedback below only applies 
c     to gas particles
c
         IF (iphase(ipart).EQ.0) THEN
c
c--Compute rate of ionisation
c
            IF (istellar_ionisation) THEN
               dHIIdt(ih2pos,ipart) = 0.

               dti = (dtmax/imaxstep)*isteps(ipart)
               DO isink = 1, nptmasstot
                  h2_over_dist2 = xyzmh(5,ipart)**2/(
     &                 (xyzmh(1,ipart)-stellar_xyzmhrti(1,isink))**2 + 
     &                 (xyzmh(2,ipart)-stellar_xyzmhrti(2,isink))**2 + 
     &                 (xyzmh(3,ipart)-stellar_xyzmhrti(3,isink))**2)

                  dHIIdt(ih2pos,ipart) = dHIIdt(ih2pos,ipart) + 
     &                 h2_over_dist2 * MAX(0.0,
     &                 stellar_xyzmhrti(8,isink)/(4.0*pi)*(gmw/gmw*mH) - 
     &                 stellarrad(3,isink,ipart)*2.7E-13/(gmw/gmw*mH) 
     &                 - stellarrad(4,isink,ipart)/(dti*utime)
     &                 )
               END DO
c
c--Account for recombinations, and convert to code units
c
               dHIIdt(ih2pos,ipart) = (dHIIdt(ih2pos,ipart)
c     &              - (trho(ipart)*udens/mH)*2.7E-13*HIIion(ipart)*
c     &              xyzmh(4,ipart)*umass
     &              ) /
     &              (xyzmh(4,ipart)*umass/1.0)*utime
c
c--Local recombination rate (g/s) for a fully ionised particle,
c     in code units
c
               recomb = (trho(ipart)*udens/mH)*2.7E-13*
     &              xyzmh(4,ipart)*umass
     &              / (xyzmh(4,ipart)*umass/1.0)*utime

               IF (dHIIdt(ih2pos,ipart).EQ.0.) THEN
                  dHIIdt(ih2pos,ipart) = dHIIdt(ih2pos,ipart)
     &                 - HIIion(ipart)**2  * recomb
               ENDIF

            ENDIF
c
c--For dust diffusion in disc, need omega_disc values
c
            IF (idustFluid.NE.0) THEN
               IF (iexf.EQ.11) THEN ! For dust settling test only
                  omega_disc(ipart) = SQRT(xmass/((5.0)**3))
               ELSE
c
c--The rgrav_accel is equal to Omega*v_Keplerian
c     Note: this is not properly defined until forcei has been called
c     but this is okay because dust evolution is not called for icall=1
c
                  rgrav_accel_code_units = SQRT(grav_accel(1,ipart)**2 +
     &                 grav_accel(2,ipart)**2 +
     &                 grav_accel(3,ipart)**2)
                  IF (rgrav_accel_code_units.EQ.0.) THEN
                     WRITE (*,*) 'ERROR - DissSour - rgrav_accel.EQ.0.'
                     CALL quit(0)
                  ENDIF
                  radius = SQRT(xyzmh(1,ipart)**2 + xyzmh(2,ipart)**2
     &                 + xyzmh(3,ipart)**2)
                  omega_disc(ipart) = SQRT(rgrav_accel_code_units/
     &                 (radius + tiny))
c
c--Determine whether or not particle is in a disc (determines whether
c     or not the relative velocities for dust grain growth include:
c       pressure gradient driven terminal velocity (envelope only)
c       turbulent relative velocities (disc only)
c       radial drift (disc only)
c     Based on ratio of radial velocity to tangential velocity (with
c       the direction being defined by the local gravitational 
c       acceleration), and the ratio of the tangential velocity to
c       the circular orbital speed (defined by local gravitational 
c       acceleration and radius).
c
                  vel_dot_grav_accel = (
     &                 vxyzu(1,ipart)*grav_accel(1,ipart) +
     &                 vxyzu(2,ipart)*grav_accel(2,ipart) +
     &                 vxyzu(3,ipart)*grav_accel(3,ipart)) /
     &                 rgrav_accel_code_units
                  vel2 = vxyzu(1,ipart)**2 + vxyzu(2,ipart)**2 +
     &                 vxyzu(3,ipart)**2
                  vel_tangential = SQRT( vel2 - vel_dot_grav_accel**2 )
                  ratio_vrad_to_vtan = vel_dot_grav_accel/vel_tangential
                  ratio_vrad_to_cs = vel_dot_grav_accel/vsound(ipart)

c         ratio_vtan_to_circular_orbit = vel_tangential/
c     &        SQRT(rgrav_accel_code_units*(radius+solarr/udist))
c         HY09discflag(ipart) = ratio_vtan_to_circular_orbit
c
c--Disc flag is that ABS(ratio_vrad_to_vtan).LT.0.5
c     I also tried: ratio_vtan_to_circular_orbit.GT.0.6 and
c     ABS(ratio_vrad_to_vtan).LT.0.5
c     but the additional criterion only helps at the outer boundary
c     of the cloud where dust coagulation is negligible anyway.
c     This seems to provide good disc extraction, for example, 
c     excluding particles in the high "atmosphere" of a disc.
c
c         IF (ratio_vtan_to_circular_orbit.GT.0.6) THEN
                  IF (ABS(ratio_vrad_to_vtan).GT.0.5) THEN
                     omega_disc(ipart) = 1.0E-10
c
c     As pointed out in Bate (2022), when including envelope turbulence
c     in non-rotating cloud calculations it is unreasonable to have the 
c     first hydrostatic core have envelope-type turbulence rather than
c     more quiescent disc-type turbulence, so the additional criteria
c     below was added to use disc-type turbulence if the radial
c     velocity is smaller than 10% of the sound speed.  This selects
c     the first hydrostatic core in non-rotating calculations as being
c     disc-type turbulence.
c
                  ELSEIF (ABS(ratio_vrad_to_cs).GT.0.1) THEN
                     omega_disc(ipart) = 2.0E-10
c
c     Disc scale height must be less than 100 au for any reasonable disc
c
                  ELSEIF (vsound(ipart)/omega_disc(ipart).GT.
     &                    100.0*au/udist) THEN
                     omega_disc(ipart) = 3.0E-10
                  ENDIF
               ENDIF
            ENDIF
#ifdef RT
c
c--Evolution of H_2
c
            IF (idustRT.GT.0) THEN
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
     &                    - 0.5*vsigi*hi1*psi
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
