      subroutine potential(x,y,z,ti,p)

      real r, phi, gamma,Kn, Bn,Dn,sum
      real x,y,p,ti,t0


      integer n

      INCLUDE 'COMMONS/potent'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'

c     Convert to cylindrical polars:

      if (x.NE.0 .and. y.NE.0) then
         phi=atan(y/x)
         if (x.LT.0) phi=pi+phi
      else
         if (x.EQ.0 .and. y.GT.0) phi=pi/2.
         if (x.EQ.0 .and. y.LT.0) phi=-pi/2.
         if (y.EQ.0 .and. x.GT.0) phi=0
         if (y.EQ.0 .and. x.LT.0) phi=pi
      end if

      r=sqrt(x**2.+y**2.)

c     Time at initial conditions (100Myr)

      t0=3.153e+15/utime

c     Calculate the potential

      gamma=NN*(phi+phir*(t0+ti)-log(r/r0)/tan(alpha))
c      gamma=NN*(phi+pi/2.-log(r/r0)/tan(alpha))

      sum=0

      do n=1,3

         Kn=n*NN/(r*sin(alpha))
         Bn=Kn*Hz*(1.+0.4*Kn*Hz)
         Dn=(1.+Kn*Hz+0.3*(Kn*Hz)**2.)/(1.+0.3*Kn*Hz)

         sum=sum+(Cz(n)/(Dn*Kn))*cos(n*gamma)
     &        *(1/cosh((Kn*z)/Bn))**Bn

      end do
         
      spiral=-4.*pi*Hz*p0*exp(-(r-r0)/rS)*sum
c      logpart=Co*log(Rc**2.+r**2.+(z/0.7)**2.)
     
      p=spiral


c      OPEN (19, FILE = 'consts', STATUS = 'UNKNOWN')
c      WRITE(19,*)phir,ti,phir*ti,spiral,logpart

      RETURN
c     end subroutine potentialsub
      END
