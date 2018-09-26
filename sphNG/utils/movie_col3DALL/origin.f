      SUBROUTINE origin (n1, n2, frac)
c************************************************************
c                                                           *
c  this subroutine computes the origin of the coordinate    *
c  system.                                                  *
c                                                           *
c************************************************************
c
      INCLUDE 'idim'
      PARAMETER (iline=301)
c
      INCLUDE 'COMMONS/part'

      COMMON /origi/ cmx, cmy, cmz
c
c--center of mass
c
      cmx = 0.
      cmy = 0.
      cmz = 0.
      tmass = 0.
c
      DO 10 i=1,npart
         tmass = tmass + xyzmh(4,i)
         cmx = cmx + xyzmh(4,i)*xyzmh(1,i)
         cmy = cmy + xyzmh(4,i)*xyzmh(2,i)
   10    cmz = cmz + xyzmh(4,i)*xyzmh(3,i)

      cmx = cmx/tmass
      cmy = cmy/tmass
      cmz = cmz/tmass
c
      DO 20 i=1,npart
         xyzmh(1,i) = xyzmh(1,i) - cmx*frac
         xyzmh(2,i) = xyzmh(2,i) - cmy*frac
         xyzmh(3,i) = xyzmh(3,i) - cmz*frac
   20 CONTINUE

      WRITE(*,*) cmx,cmy,cmz
c
      RETURN
      END
