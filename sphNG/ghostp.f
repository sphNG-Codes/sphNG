      SUBROUTINE ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz, 
     &     dustvar)
c************************************************************
c                                                           *
c     This subroutine acts as an interface to all the       *
c     routines which compute the list of ghost particles    *
c       for treating the boundaries.                        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(imhdevol,imhd)
      DIMENSION dustvar(idim_dustFluid)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/typef'

      nghost = 0
      IF (ibound.EQ.1) CALL ghostp1(npart,xyzmh,vxyzu,ekcle,Bevolxyz,
     &     dustvar)
      IF (ibound.EQ.2) CALL ghostp2(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
      IF (ibound.EQ.3 .OR. ibound.EQ.8 .OR. ibound/10.EQ.9) 
     &        CALL ghostp3(npart,xyzmh,vxyzu,ekcle,Bevolxyz,dustvar)
      IF (ibound.EQ.100) 
     &        CALL ghostp100(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
      IF (ibound.EQ.11) 
     &        CALL ghostp11(npart,xyzmh,vxyzu,ekcle,Bevolxyz,dustvar)
      IF (ibound.EQ.102 .OR. ibound.EQ.103) 
     &        CALL ghostp102(npart,xyzmh,vxyzu,ekcle,Bevolxyz,dustvar)
      IF (ibound.EQ.104) 
     &        CALL ghostp104(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
      ntot = npart + nghost
      
      RETURN
      END
