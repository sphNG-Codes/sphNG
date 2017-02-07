      subroutine dragtable

c************************************************************
c                                                           *
c  Routine to setup a table to record energy lost by solids *
c  being dragged by gas, so as to record energy release as  *
c  they are accreted by planets.                            *
c                                                           *
c  Output stored in a 3-dimensional grid using spherical    *
c  polar coordinates as the indices.                        *
c                                                           *
c  Ben Ayliffe, 17th November 2011.                         *
c                                                           *
c************************************************************

#ifdef MPIALL
#define IMPLICITNONE
#include "mpi_sup.h"
#else
      IMPLICIT NONE
#endif

      INCLUDE 'idim'

      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/logun'

c
c--Setup r coordinate
c      
      dragtrmax = log10(5.*(coremass_orig/(3.*xmass))**(1./3.))
      dragtrmin = log10(rplanet)
      dragtdr = (dragtrmax-dragtrmin)/(1.*idragresr)
c
c--Setup phi coordinate
c      
      dragtpmax = pi
      dragtpmin = -pi
      dragtdp = (dragtpmax-dragtpmin)/(1.*idragresp)
c
c--Setup theta coordinate
c      
      dragttmax = pi
      dragttmin = 0.0
      dragtdt = (dragttmax-dragttmin)/(1.*idragrest)

      print *, dragtrmin, dragtrmax, dragtdr
      print *, dragtpmin, dragtpmax, dragtdp
      print *, dragttmin, dragttmax, dragtdt

#ifdef MPICOPY
      IF ((2.*idragresr*idragresp*idragrest).GT.(idim*6+imhd*9)) THEN
         WRITE (iprint,*) 'realtransfer5to15 array used for MPICOPY ',
     &     'communication is insufficiently large for the dragtable',
     &     ' resolution. Enlarge idim to accomodate it, or rewrite',
     &     ' the code to use an alternative array.'
         STOP
      ENDIF
#endif

      end subroutine dragtable


      subroutine dragcoordinates(radius, phi, theta, ir, ip, it)
      
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/radtran3'

      REAL radius, phi, theta, logradius
      INTEGER ir, ip, it

      IF (phi.GT.dragtpmax .OR. theta.GT.dragttmax) THEN
         print *, 'Angles beyond drag table limits'
         print *, 'Check angle setup: +-1 or something'
         print *, phi, dragtpmax
         print *, theta, dragttmax
         STOP
      ENDIF

      logradius = log10(radius)

      IF (logradius.GE.dragtrmax) THEN
         ir = idragresr+1
         ip = 1
         it = 1
         RETURN
      ENDIF

      ir = INT((logradius - dragtrmin)/dragtdr + 1)
      ip = INT((phi - dragtpmin)/dragtdp + 1)
      it = INT((theta - dragttmin)/dragtdt + 1)

c--Bodies may slightly penetrate surface, thus ir < 1
      ir = MAX(1,ir)

      IF (ir.GT.idragresr .OR. ip.GT.idragresp .OR. it.GT.idragrest. OR.
     &     ir.LT.1 .OR. ip.LT.1 .OR. it.LT.1) THEN
         print *, 'Index beyond allowed range'
         print *, 'ir, ip, it : ', ir, ip, it
         print *, 'res r, p, t: ', idragresr, idragresp, idragrest
         STOP
      ENDIF

      end subroutine dragcoordinates
