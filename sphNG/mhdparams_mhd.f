      SUBROUTINE mhdparams
c************************************************************
c                                                           *
c  This routine computes MHD quantities for printout        * 
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/outmhd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/numpa'
      REAL*8 etai
      REAL rhoi
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine mhdparams')

      betamhdav = 0.
      betamhdmax = 0.
      betamhdmin = 1.E30
      divBmax = 0.
      divBav = 0.
      curlBav = 0.
      curlBmax = 0.
      div2curlBav = 0.
      div2curlBmax = 0.
      omegamhdav = 0.
      omegamhdmax = 0.
      omegtol = 1.E-2
      fracdivBok = 0.
      fluxtot = 0.
      crosshel = 0.
      Bmin = 1.E30
      Bmean = 0.
      Bmax = 0.
      valphaBmax = -1.
      valphaBmin = 10.0
      etamax = 0.
      etamin = 1.E30
      etaav = 0.
      etarealoverartav = 0.
      etarealoverartmin= 1.E30
      etarealoverartmax = 0.
      ietacount = 0
      
c
c--Calculate quantities
c
      DO i=1,npart
         rhoi = rho(i)
         rho1i = 0.
         IF (rhoi.NE.0) rho1i = 1./rhoi
c         IF (varmhd.EQ.'Bvol') THEN
c            Bxi = Bevolxyz(1,i)
c            Byi = Bevolxyz(2,i)
c            Bzi = Bevolxyz(3,i)         
c         ELSEIF (varmhd.EQ.'Brho') THEN
c            Bxi = Bevolxyz(1,i)*rho(i)
c            Byi = Bevolxyz(2,i)*rho(i)
c            Bzi = Bevolxyz(3,i)*rho(i)
c         ELSE
            Bxi = Bxyz(1,i)
            Byi = Bxyz(2,i)
            Bzi = Bxyz(3,i)
c         ENDIF

         B2i = Bxi**2 + Byi**2 + Bzi**2
         Bi = SQRT(B2i)
         IF (Bi.GT.0.) THEN
            Bi1 = 1./Bi
            betamhdi = pr(i)/(0.5*B2i)
         ELSE
            Bi1 = 0.
            betamhdi = 0.
         ENDIF
c
c--Max/min/ave B
c         
         IF (Bi.LT.Bmin) Bmin = Bi
         IF (Bi.GT.Bmax) Bmax = Bi
         Bmean = Bmean + Bi
c
c--Plasma beta minimum/maximum/average
c         
         betamhdav = betamhdav + betamhdi
         IF (betamhdi.GT.betamhdmax) betamhdmax = betamhdi
         IF (betamhdi.LT.betamhdmin) betamhdmin = betamhdi
c
c--Maximum/average divergence of B
c         
         divBi = abs(divcurlB(1,i))
         IF (divBi.GT.divBmax) divBmax = divBi
         divBav = divBav + divBi
c
c--max/av magnitude of current
c
         curlBxi= divcurlB(2,i)
         curlByi= divcurlB(3,i)
         curlBzi= divcurlB(4,i)
         curlBmag2i= curlBxi*curlBxi + curlByi*curlByi + curlBzi*curlBzi
         curlBmagi= SQRT(curlBmag2i)
         curlBav= curlBav + curlBmagi
         IF (curlBmagi.GT.curlBmax) curlBmax= curlBmagi
c
c--ratio of divB to curlB
c
         div2curlBi= divBi/(curlBmagi + tiny)
         div2curlBav= div2curlBav + div2curlBi
         IF (div2curlBi.GT.div2curlBmax) div2curlBmax= div2curlBi
c
c--|div B| x smoothing length / |B| (see e.g. Cerqueira and Gouveia del Pino 1999) 
c  this quantity should be less than ~0.01.
c
         omegamhdi = divBi*xyzmh(5,i)*Bi1
         IF (omegamhdi.LT.omegtol) fracdivBok = fracdivBok + 1.
         IF (omegamhdi.GT.omegamhdmax) omegamhdmax = omegamhdi
         omegamhdav = omegamhdav + omegamhdi
c
c--Conserved magnetic flux (int B dV)
c
         pmassi = xyzmh(4,i)
         fluxtot = fluxtot + pmassi*divcurlB(1,i)*rho1i
c
c--Conserved Cross Helicity (int v.B dV)
c
         crosshel = crosshel + pmassi*(vxyzu(1,i)*Bxi 
     &                         + vxyzu(2,i)*Byi + vxyzu(3,i)*Bzi)*rho1i
c
c--magnetic resistivity parameter
c
         valphaBmax = MAX(valphaBmax, alphaMM(2,i))
         valphaBmin = MIN(valphaBmin, alphaMM(2,i))
         
c         IF (iresist.GE.2) THEN
            etai = etafunc(rhoi,vxyzu(4,i),etamhd,iresist)
            etamax = MAX(etamax, etai)
            etamin = MIN(etamin, etai)
            etaav = etaav + etai
            
            vsigi = sqrt(vsound(i)**2 + B2i*rho1i)
            IF (ifsvi.EQ.6) THEN
               etaart = alphaMM(2,i)*vsigi*xyzmh(5,i)
            ELSE
               etaart = alphamin(2)*vsigi*xyzmh(5,i)
            ENDIF

            etaartificial(i) = etaart
            etareal(i) = etai

            IF (etaart.GT.0.) THEN
               ratio = etai/etaart
               etarealoverartmin = min(etarealoverartmin,ratio)
               etarealoverartmax = max(etarealoverartmax,ratio)
               etarealoverartav = etarealoverartav + ratio
               ietacount = ietacount + 1
            ENDIF
c         ENDIF
      ENDDO
      
      denom = 1./FLOAT(npart)
      betamhdav = betamhdav*denom
      fracdivBok = 100.*fracdivBok*denom
      omegamhdav = omegamhdav*denom
      divBav = divBav*denom
      curlBav = curlBav*denom
      div2curlBav = div2curlBav*denom
      Bmean = Bmean*denom
      etaav = etaav*denom
      ietacount = max(ietacount,1)
      etarealoverartav = etarealoverartav/FLOAT(ietacount)

      RETURN
      END
