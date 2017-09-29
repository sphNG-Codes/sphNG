      SUBROUTINE addump
c************************************************************
c                                                           *
c  This routine reads two existing dumps and create a third *
c  one out of it.                                           *
c                                                           *
c************************************************************

      INCLUDE 'idim' 

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/active'

      CHARACTER*7 where
      CHARACTER*1 iok
      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums(8)

      DATA where/'addump'/
c
c--allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine addump')

      nstep = 1

      fmas1 = 0.
      fmas2 = 0.
      WRITE (*, 99002)
99002 FORMAT (' PARTICLE SET UP', //)
      WRITE (*, 99003)
99003 FORMAT (' self-gravity included? (y/n)')
      READ (*, 99004) iok
99004 FORMAT (A1)
      igphi = 0
      IF ( iok.EQ.'y' ) igphi = 1
      WRITE (*, 99005)
99005 FORMAT (' give the position and the velocity of the cm of body 1')
      READ (*, *) xx1, yy1, zz1m, vvx1, vvy1, vvz1
      WRITE (*, 99006)
99006 FORMAT (' give the position and the velocity of the cm of body 2')
      READ (*, *) xx2, yy2, zz2, vvx2, vvy2, vvz2
      WRITE (*, 99007)
99007 FORMAT (' do you want to rotate body 2? (y/n)')
      READ (*, 99004) iok
      rotang = 0.
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99008)
99008    FORMAT (' give rotation angle (rad.)')
         READ (*, *) rotang
      ENDIF
c
c--read first object
c
      WRITE (*,*) 'Reading first file and storing data...'
      CALL rdump (idisk1, ichkl, 0)
      np1 = npart
c
c--put first object into place
c
      DO 100 j = 1, np1
         fmas1 = fmas1 + xyzmh(4,j)
         xyzmh(1,j) = xyzmh(1,j) + xx1
         xyzmh(2,j) = xyzmh(2,j) + yy1
         xyzmh(3,j) = xyzmh(3,j) + zz1
         vxyzu(1,j) = vxyzu(1,j) + vvx1
         vxyzu(2,j) = vxyzu(2,j) + vvy1
         vxyzu(3,j) = vxyzu(3,j) + vvz1
 100  CONTINUE
c
c--read second dump
c
      WRITE (*,*) 'Reading second file and storing data...'
      CALL rdump (idisk2, ichkl, np1)
      np2 = npart
c
c--rotate second object if needed
c
      crotang = cos(rotang)
      srotang = sin(rotang)
      DO 200 j = np1 + 1, np1 + np2
         fmas2 = fmas2 + xyzmh(4,j)
         xyzmh(1,j) = crotang*xyzmh(1,j) + srotang*xyzmh(2,j) + xx2
         xyzmh(2,j) = -srotang*xyzmh(1,j) + crotang*xyzmh(2,j) + yy2
         xyzmh(3,j) = xyzmh(3,j) + zz2
         vxyzu(1,j) = vxyzu(1,j) + vvx2
         vxyzu(2,j) = vxyzu(2,j) + vvy2
         vxyzu(3,j) = vxyzu(3,j) + vvz2
 200  CONTINUE

      WRITE (*, 99011)
99011 FORMAT (' what is the equation of state variable:', /,
     &        ' specific internal energy :  intener', /,
     &        ' specific entropy         :  entropy')
      READ (*, 99012) varsta
99012 FORMAT (A7)
      WRITE (*, 99013)
99013 FORMAT (' is the equation of state a gamma-law? (y/n)')
      READ (*, 99004) iok
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99014) gamma1, gamma2
99014    FORMAT (' what is the value of the adiabatic index gamma?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) gamma
         WRITE (*, 98014) RK21, RK22
98014    FORMAT (' what is the value of the adiabatic constant RK2?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) RK2
         WRITE (*, 98015) rhozero1, rhozero2
98015    FORMAT (' what is the value of the  constant rhozero?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) rhozero
      ENDIF
      npart = np1 + np2
      nactive = npart
      n1 = np1
      n2 = np2
      WRITE (*, 99015) n1, n2, npart, nptmass
99015 FORMAT (' set up completed', /,
     &        ' particles in object 1          :', I6, /,
     &        ' particles in object 2          :', I6, /,
     &        ' total number of particles used :', I6, I4)
      WRITE (*, 99016)
99016 FORMAT (//, ' END SET UP')

      RETURN

 300  CALL error(where, 1)

      RETURN
      END
