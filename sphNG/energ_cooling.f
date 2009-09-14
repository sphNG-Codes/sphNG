      SUBROUTINE energ_cooling(ipart, npart, ntot,ti, vxyzu, dvxyzu,
     & xyzmh,h2ratio,delt,abhpq,abeq,abHIq,abco,trho)

c************************************************************
c                                                           *
c  This routine computes the change in internal energy or   *
c     in entropy.                                           *
c                                                           *
c    varsta = entropy : compute change in specific entropy  *
c    varsta = intener : compute change in specific internal *
c                          energy                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/tcooling'
      INCLUDE 'COMMONS/h2'

      DIMENSION vxyzu(4,idim), dvxyzu(4,idim), xyzmh(5,mmax)
      DIMENSION h2ratio(idim_h2),abund(9),abco(idim_h2)
      DIMENSION abeq(idim_h2),abhpq(idim_h2),abHIq(idim_h2),ratesq(11) 
      REAL omega1,gamma1,chi,dlnrdt,rscale,drdt
      REAL np1,nh1,nh21,dtclare,nh1t,cdensH2,cdens,tstep
      REAL rate,nH2,rate_diss,cdens2,dttest,ylamq,dlq
      REAL dphot,dchem,xfac,tstep2,tsteptest
      REAL delt,tempiso,pdv,coeff1,coeff2,th2,h2mol,nmol
      REAL abundc,abundo,abundsi,abc,gmwvar,dqexp
      REAL abcp,phrates,phrates2,k0,k1,gamma_chx,gamma_co,beta
      REAL*4 trho(idim)
      INTEGER nstep, ipart, nstep2
c
c--Reads is h2ratio,length of current timestep and abundances of protons, electrons,HI, CO
c--Also reads in density (trho) 
c
c
c--abundances of C, O and Si and ionization needed for Simon's chemistry
c
      abundc = 2d-4
      abundo = 4.5d-4
      abundsi = 3d-5
c
c--some variables needed for CO chemistry
c
      k0=5d-16
      k1=5d-10
c
c--For MPI code do not need to add energy changes that are added only added
c     by the process that actually holds the particle
c
      IF (ipart.GT.ntot) THEN
         xfac = 0.0
         ipartntot = ipart + ntot + 2
      ELSE
         xfac = 1.0
         ipartntot = ipart
      ENDIF
c
c--Initialisation
c
      gamma1 = gamma - 1.
c
c--pdv term is stored in dvxyzu(4,ipart) in forcei to use memory efficiently
c
      pdv = dvxyzu(4,ipart)
      dvxyzu(4,ipart) = 0.0

ccccccccccccccccc Beginning of chemistry cccccccccccccccccc
      IF (iener.EQ.4) THEN
      IF (xfac.EQ.1) THEN
c
ccc iener=4 means Simon's chemistry is being used
c
ccc timestep in seconds
c
         dtclare=delt*utime
c                 
ccc mean molecular calculated taking into account the molecular gas. This is 1.27 when there is no
ccc molecular hydrogen
c
         gmwvar=(2.*h2ratio(ipart)+(1.-2.*h2ratio(ipart))+0.4)/
     &      (0.1+h2ratio(ipart)+(1.-2.*h2ratio(ipart)))
c
ccc this first section is to update the H2 fraction. But this requires using a very small timestep and 
ccc subcycling. First an estimate of the timestep to update the H2 is determined. 
c
ccc np1=total number density inclusive of protons - nh1=number density of HI inclusive of protons
ccc - nh21=number density of H2
c 
         np1=(trho(ipart)*udens/mp)*5./7.
         nh1=np1*(1.-2.*h2ratio(ipart))
         nh21=np1*h2ratio(ipart)
c
ccc Temperature
c
         tempiso = 2./3.*vxyzu(4,ipart)/(Rg/gmwvar/uergg)
         IF (ipart.EQ.1) PRINT*,'densities',np1,nh1,nh21,tempiso,gmwvar
c
ccc Calculate photodissociation rate
c
ccc column density of H2
c
         cdensH2=nh21*dphot
c
ccc total column density 
c
         cdens=np1*dphot
         rate=-1.0*cdens*sigma
c
ccc Draine & Bertoldi photodissociation rate
c
         cdens2=cdensH2/5.e14
         nH2=0.965/(1.+cdens2/3.)**2.+0.035/(1.+cdens2)**0.5
         nH2=nH2*exp(-0.00085*(1.+cdens2)**0.5)
         rate_diss=nH2*exp(rate)*rate_diss0
c
         th2=10000.
c		 
ccc From Bergin et. al. 2004, d(n(h2))/dt=(formation-destruction)
ccc Here we are looking to dt such that n(H2) goes to zero, i.e dt=0-n(H2)/(formation-destruction)
ccc First just check that dt is not 0. - but actually this should not occur anyway.
c
         IF ((np1*nh1*tempiso**0.5*Rconst-
     &       (rate_diss+rate_cr)*nh21).EQ.0) THEN
			nstep=5000
         ELSE
c		 
ccc - calculate dt according to dt=0-n(H2)/(formation-destruction)
c 		 
			dttest=-h2ratio(ipart)*np1/(np1*nh1*tempiso**0.5*Rconst-
     &         (rate_diss+rate_cr)*nh21)
c
ccc - if dttest>0 then h2 is being photodissociated. Set timestep to correspond to 1/10th of the time to 
ccc completely photodissociate to 0 (this step is basically ensuring that the amount of H2 does not
ccc become negative. 
c
         IF (dttest.GT.0.) th2=0.1*dttest/dtclare 
c
ccc Set number of steps for subcycling. In the event that the h2 is 0 (i.e. for the first timestep) or h2
ccc being created, the number of timesteps is 201 (chosen by running the simulation initially and finding
ccc a sensible value). 
c 		 
         IF (th2.GT.0) nstep=MAX(INT(1./th2)+1,201)
         END IF
         tstep=dtclare/nstep

         IF (ipart.EQ.1)  PRINT*,'nh2',rate,dttest,nstep,tstep,th2
c
ccc Now actually update the H2 fraction - yes the H2 bit should be called by a separate subroutine as
ccc it was also in the last section, but I've not done that yet
c
         DO i=1,nstep
c
	      np1=(trho(ipart)*udens/mp)*5./7.
		  nh1=np1*(1.-2.*h2ratio(ipart))
	      nh21=np1*h2ratio(ipart)
c         
ccc Calculate photodissociation rate
c
ccc column density of H2
c
          cdensH2=nh21*dphot
c
ccc total column density 
c
          cdens=np1*dphot
          rate=-1.0*cdens*sigma
c
ccc Draine & Bertoldi photodissociation rate
c
          cdens2=cdensH2/5.e14
          nH2=0.965/(1.+cdens2/3.)**2.+0.035/(1.+cdens2)**0.5
          nH2=nH2*exp(-0.00085*(1.+cdens2)**0.5)
          rate_diss=nH2*exp(rate)*rate_diss0
c
ccc Determine number density of molecular hydrogen
c
		  nmol=max(np1*nh1*tempiso**0.5*Rconst*tstep-(rate_diss
     &      +rate_cr)*tstep*nh21+h2ratio(ipart)*np1,0.)
c         
ccc There should also be a timestep criterion so that H2 abundances do not exceed 1.
ccc But in my calculations so far, the ratios don't reach 1 (there is no sharp cut off
ccc in H2 fraction) indicating that the minimum 
ccc number of timesteps of 201 prevents this. But for higher density calculations,
ccc a criterion should be included. This is the number density of h2
c
		  h2mol=min(np1*0.5,nmol)
c
c and ratio of H2 to HI+H2 (i.e. n(H2)/n(HI+H2)
c
          h2ratio(ipart)=h2mol/np1
c
ccc Update abundance of CO - this is wrong at present, but CO not included in any other part of 
ccc calculations at present, so does not matter.
c
c          abcp= max(0d0, abundc - abco(ipart))
c          IF (h2ratio(ipart).EQ.0.) THEN
c             beta=0.
c              abco(ipart)=0.
c          ELSE IF (h2ratio(ipart).GT.0.) THEN
c             beta = k1 * (abundo-abco(ipart))/(k1 *(abundo-abco(ipart))+ 
c     &      gamma_chx / (h2ratio(ipart) * np1))
c             abco(ipart)= max(abco(ipart)+(k0*abcp*beta*np1**2 
c     &      -gamma_co*abco(ipart)*np1)*tstep,0.)
c          END IF
c
          abcp= max(0d0, abundc - abco(ipart))
          phrates2=phrates*exp(-0.267*sigma*cdens)
          gamma_chx=5d0*phrates2
          gamma_co=phrates2
c
          beta = k1 * (abundo-abco(ipart))/(k1 *(abundo-abco(ipart))+ 
     &       gamma_chx / (h2ratio(ipart) * np1))
          IF (abco(ipart).EQ.0.) THEN
            nstep2=INT(trho(ipart)*1000)
c 
            tstep2=tstep/nstep2
          ELSE
           tsteptest=-abco(ipart)/(k0*abcp*beta*np1**2 
     &      -gamma_co*abco(ipart)*np1)
c           nstep2=max(INT(tstep*10./tsteptest),INT(trho(ipart)*10.),1)
           nstep2=max(INT(tstep*10./tsteptest),1)
           tstep2=tstep/nstep2
          END IF
          DO j=1,nstep2
            abcp= max(0d0, abundc - abco(ipart))
            phrates2=phrates*exp(-0.267*sigma*cdens)
            gamma_chx=5d0*phrates2
            gamma_co=phrates2
           IF (h2ratio(ipart).EQ.0.) THEN
             beta=0.
             abco(ipart)=0.
           ELSE IF (h2ratio(ipart).GT.0.) THEN
             beta = k1 * (abundo-abco(ipart))/(k1 *(abundo-abco(ipart))+ 
     &       gamma_chx / (h2ratio(ipart) * np1))
             abco(ipart)= max(abco(ipart)+(k0*abcp*beta*np1**2 
     &      -gamma_co*abco(ipart)*np1)*tstep2,0.)
           END IF
          END DO
c
         END DO
c
ccc End of updating H2 ratio. Now actually need to do heating and cooling.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--If were not including H2, could set h2ratio to a small value (e.g. 1.e-7) and just
c--have this part to calculate heating and cooloing (need nh1 and np1 though).
c
ccc column density of HI excluding protons
c
         nh1t=nh1*abHIq(ipart)
c
ccc update abundances of HI, protons and electrons requires total number density
ccc and column number density of HI 
c
         CALL hchem(tempiso,np1,nh1t*dchem,abeq(ipart),
     &     abhpq(ipart),coeff1,coeff2)
         IF (ipart.EQ.1) PRINT*,'1',delt,tempiso,np1,coeff1,abhpq(1)
c
ccc New number density of HI
c
         nh1t=nh1t+(coeff1-coeff2*nh1t)*dtclare
c
ccc Calculate new abundance of HI
c
         abHIq(ipart)=nh1t/nh1
c
ccc Make sure abHI not >1 - should be taken care of by Simon's chemistry.
c
         IF (abHIq(ipart).GE.1) abHIq(ipart)=1.
c
ccc calculate abundance of protons and electrons
c
         abhpq(ipart)=1.0-abHIq(ipart)
         abeq(ipart)=abhpq(ipart)+abc
c
ccc Set abunances of H2, HI, electrons and protons for cooling routine
c
         abund(1)=h2ratio(ipart)
         abund(2)=(1-2.*h2ratio(ipart))*abHIq(ipart)
         abund(3)=(1-2.*h2ratio(ipart))*abhpq(ipart)+abc
         abund(4)=(1-2.*h2ratio(ipart))*abhpq(ipart)
         abund(5) = abundo
         abund(6) = 0d0
         abund(7) = abundc
         abund(8) = 0d0
         abund(9) = abundsi 
c
c Call cooling routine, requiring total density, some distance measure and 
c abundances
         IF (ipart.EQ.1) PRINT*,'cool',tempiso,np1,dlq,ylamq,ratesq
         CALL cool_func(tempiso,np1,dlq,abund,ylamq,ratesq)
c     
ccc compute change in u
c
         dvxyzu(4,ipart)=(-1.*ylamq/(trho(ipart)*udens))
     &     *utime**3./udist**2.
         IF (ipart.EQ.1) PRINT*,'dvxyzu',dvxyzu(4,1)
         END IF
         END IF

ccccccccccccccc End of Chemistry (oh no, must be Exeter!)cccccccccc
c--rest of subroutine essentially the same apart from adding iener=4 to shock heating and pdv.
c
c--dq from expansion
c
      chi = 4.0 - 3.0*gamma
      CALL scaling(ti, rscale, drdt, dlnrdt)
      dqexp = chi*dlnrdt

      IF (varsta.NE.'entropy') THEN
         IF (iener.EQ.2) THEN
c
c--Cooling instability
c
c            critlambda = 0.1
c            rhoratio = 100.0
c            rho_zero = 0.001
c            rhocrit = 3.0*rho_zero

            critlambda = 0.1
            rhoratio = 10.0
            rho_zero = 0.001
            rhocrit = 3.0*rho_zero

            dvxyzu(4,ipart) = 1.5/critlambda*(1.0 - trho(ipart)/rhocrit)
            IF (vxyzu(4,ipart).GE.1.5 .AND. dvxyzu(4,ipart).GT.0.0) 
     &           dvxyzu(4,ipart) = 0.0
            IF (vxyzu(4,ipart).LE.1.5/rhoratio .AND. 
     &           dvxyzu(4,ipart).LT.0.) dvxyzu(4,ipart) = 0.0
         ELSEIF ((iener.EQ.1).OR.(iener.EQ.3).OR.(iener.EQ.4)) THEN

c
c--Compute change in specific internal energy
c
c  a) pdv term first
c
            IF (iexpan.EQ.0) THEN
               dvxyzu(4,ipart) = dvxyzu(4,ipart)+pdv*pr(ipart)/
     &              trho(ipart)**2

            ELSE
               dvxyzu(4,ipart) = dvxyzu(4,ipart)+pdv*pr(ipart)/
     &              trho(ipart)**2 - vxyzu(4,ipart)*dqexp*xfac
            ENDIF
c
c  b) Shock dissipation if appropriate
c
            IF (ichoc.NE.0) THEN
               dvxyzu(4,ipart) = dvxyzu(4,ipart) + 0.5*dq(ipart)
            ENDIF

c
c  c) Cooling for adiabatic-with-cooling eos
c
            IF(encal.EQ.'c' .AND. xfac.EQ.1.0) THEN
               r = SQRT(xyzmh(1,ipartntot)**2 + xyzmh(2,ipartntot)**2 + 
     &               xyzmh(3,ipartntot)**2 )

c Omega1 is 1/omega
c !! I think listpm(1) should be the array index of the 
c !! central object. If not, this is wrong.
               omega = sqrt(xyzmh(4,listpm(1))/ r**3)
               omega1 = 1./omega

c Cooling times used by Rice et al 2003 are 5/omega and 3/omega
               tcool = coolingrate * omega1

c               PRINT *, dvxyzu(4,ipart),vxyzu(4,ipart)/tcool,omega
               dvxyzu(4,ipart) = dvxyzu(4,ipart) - vxyzu(4,ipart) 
     &              / tcool
            ENDIF
         ENDIF
c
c--Compute change in specific entropy
c
c  a) Shock dissipation
c
      ELSEIF (ichoc.NE.0) THEN
          dvxyzu(4,ipart) = vxyzu(4,ipart)*gamma1*dq(ipart)*
     &        trho(ipart)*0.5/pr(ipart) - pr(ipart)*dqexp*xfac
      ENDIF

         IF (ipart.EQ.1) PRINT*,'dvxyzu2',dvxyzu(4,1)
      RETURN
      END
