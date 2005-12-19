      SUBROUTINE evol
c************************************************************
c                                                           *
c  This routine drives the time evolution of the system.    *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/debpt'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/accurpt'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/debugit'
      INCLUDE 'COMMONS/secret'
      INCLUDE 'COMMONS/accrem'
      INCLUDE 'COMMONS/avail'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/current'

      CHARACTER*7 where

      DATA where/'evol'/
c
c--Define options of this run
c
      imax = 1073741824
      imaxstep = imax/2

      nlistinactive = 0

      CALL options
c
c--Set all quantities needed for the run
c
      CALL preset
c
c--Initialise random number generator
c
      iseed = -4357
      rnd = ran1(iseed)
      WRITE(iprint,*)'Random seed: ',iseed
c
c--Set number of active particles
c
      nactive = 0
      DO i = 1, npart
         IF (iphase(i).NE.-1) nactive = nactive + 1
      END DO

      DO i = 1, nptmass
         j = listpm(i)
         listrealpm(j) = i
         pmasspt = xyzmh(4,j)
         ptmsyn(i) = pmasspt
         ptmadd(i) = 0.0
         xmomsyn(i) = pmasspt*vxyzu(1,j)
         ymomsyn(i) = pmasspt*vxyzu(2,j)
         zmomsyn(i) = pmasspt*vxyzu(3,j)
         xmomadd(i) = 0.0
         ymomadd(i) = 0.0
         zmomadd(i) = 0.0
      END DO
c
c--Write all quantities on listings
c
      CALL header(where)
c
c---------------------------
c---- E V O L U T I O N ----
c---------------------------
c
      ncount = 0
      icount = 0
      nbuild = 0
      iaccr = 0
      nlstacc = 0
      iptout = 0
      iptdump = 0
      isheld = .FALSE.
      tkeep = 4.0
ccc      tkeep = 0.0

      itest = 0

      DO i = 1, idim
c         torqt(i) = 0.0
c         torqg(i) = 0.0
c         torqp(i) = 0.0
c         torqv(i) = 0.0
c         torqc(i) = 0.0
         iremove(i) = -1
         iavail(i) = 0
         iscurrent(i) = .FALSE.
      END DO

      DO i = 1, 1000000
c
c--Evolve one timestep
c
         CALL integs
c
c--Check if saving time has come
c
         CALL save

      END DO

      RETURN
      END
