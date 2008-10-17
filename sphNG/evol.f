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
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/sync'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/gtime'

      CHARACTER*7 where

      DATA where/'evol'/
c
c--Define options of this run
c
      imax = 1073741824
      imaxstep = imax/2

      CALL options
c
c--Set all quantities needed for the run
c
      CALL preset(0)
c
c--Initialise random number generator
c
      IF (gt.EQ.0.0) THEN
         iseed = -4357
         rnd = ran1(iseed)
         WRITE(iprint,*)'Random seed: ',iseed
      ENDIF
c
c--Set number of active particles
c
      nactive = 0
      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            nactive = nactive + 1
            rhomaxsync = MAX(rhomaxsync,rho(i))
         ENDIF
      END DO      

      DO i = 1, idim
c         torqt(i) = 0.0
c         torqg(i) = 0.0
c         torqp(i) = 0.0
c         torqv(i) = 0.0
c         torqc(i) = 0.0
         iremove(i) = -1
         iscurrent(i) = .FALSE.
         it1(i) = 0
         listrealpm(i) = 0
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

c Initialise quantities for cooling curve if required

      IF (iener.EQ.3) THEN
 
       CALL thermeq

      END IF

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
      iresort = 0
      isheld = .FALSE.
      tkeep = 4.0
ccc      tkeep = 0.0

      itest = 0

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
