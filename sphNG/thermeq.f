      SUBROUTINE thermeq
c************************************************************
c                                                           *
c     This subroutine calculates the equilibrium            *
c              temperature for different densities          *
c and initialises the proton mass and heatin rate           *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/tcooling'

c Parameters for ccoling curve:
c Mass of proton

      mpenerg=1.67e-24
      
c Heating rate from Vazquez-Semadeni 2006, divided by hydrogen mass
c as required for energy equation (Koyama & Inutsuka 2002).
c Koyama & Inutsuka 2002 use d(rho*u)/dt rather than du/dt which
c is used here for SPH.
      
      heatingr=2.0e-26/mpenerg

      do i=1,35100
        tguess(i)=(i-1)/10000.+1.0
        tguess(i)=10**(tguess(i))
	tg=tguess(i)
	cooling(i)=1.0e7*exp(-1.184*1e5/(tg+1000.))
	cooling(i)=1./(cooling(i)+1.4e-2*sqrt(tg)*exp(-92./tg))
       end do
	   
      RETURN
      END
