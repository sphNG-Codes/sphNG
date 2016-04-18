      SUBROUTINE makedust
c************************************************************
c                                                           *
c  This subroutine reads in an existing dump file and       *
c     brings dead particles to life as dust particles.      *
c     Originally used to add dust to a settled disc.        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/logun'

      CHARACTER*11 ifile, ofile
      CHARACTER*1 igeom, iok
      
1000  FORMAT (A7)
1001  FORMAT (A1)
1002  FORMAT (A11)

#ifdef MPI
c
c--IF MPI setup, then read set up input from file
c
      iread = 10
      OPEN (iread, FILE='makedust.txt', STATUS='unknown')
#endif
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF
      PRINT *, 'Name of binary dump file to add dust to ?'
      READ (iread, 1000) ifile

      print *,iread,ifile

      PRINT *, 'Name of output binary dump file ?'
      READ (iread, 1000) ofile

      print *,iread,ofile

c      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'Name of corresponding ASCII ifile ?'
      READ (iread, 1002) inname

      PRINT *, 'Do you want to reset time to (almost) zero ?'
      READ (iread, 1001) iok

      print *,iread,inname



      imax = 1073741824
      imaxstep = imax/2

      nactive = 0
      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO

      PRINT *, 'About to open ',ifile

      OPEN (idisk1, FILE = ifile, FORM = 'unformatted')
      PRINT *, 'reading file ', ifile
c
c--Read dump file
c
      CALL options

      file1 = ofile
      print *, 'Name0=', file1

      CALL rdump_wrapper(idisk1,ichkl,1)
      CLOSE (idisk1)
      IF (ichkl.EQ.1) THEN
         PRINT*, 'ERROR READING DUMP FILE'
         CALL quit(0)
      ENDIF
c
c--End reading of dump file
c--------------------------
c
      DO i = 1, npart
         IF (iphase(i).GE.0) nactive = nactive + 1
      END DO

      min_rplan = 1.0
      max_rplan = 3.0

      PRINT *,'Enter the minimum and maximum radii for the dust',
     &     ' (e.g. ',min_rplan,max_rplan,' )'
      READ (iread,*) min_rplan, max_rplan

      PRINT *,'Do you want cylindrical or spherical radius?'
      READ (iread,1001) igeom

      gas_mass_inside = 0.
      ngas_inside = 0.
      DO i = 1, npart
         IF (iphase(i).EQ.0) THEN
            IF (igeom.EQ.'s') THEN
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            ELSE
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
            ENDIF
            IF (SQRT(r2).GT.min_rplan.AND.SQRT(r2).LT.max_rplan) THEN
               gas_mass_inside = gas_mass_inside + xyzmh(4,i)
               ngas_inside = ngas_inside + 1
            ENDIF
         ENDIF
      END DO
      PRINT *,''
      PRINT *,'Total gas mass inside dust radii ',gas_mass_inside
      PRINT *,'Total gas particles inside dust radii ',ngas_inside

      PRINT *,''
      PRINT *,'There are ',nactive,' active particles ',
     &     'out of ',npart,' total'
 50   PRINT *,'How many dust particles do you want?'
      READ (iread,*) ndustpart
#ifdef MPI
      ndustpart = ndustpart/numproc
#endif
      IF (nactive+ndustpart.GT.idim) THEN
         PRINT *,'ERROR - number must be <= to ',idim-nactive
         GOTO 50
      ENDIF

      PRINT *,'Enter desired dust to gas ratio (normally < 1)'
      READ (iread,*) dust_gas_ratio
      dust_particle_mass = dust_gas_ratio*gas_mass_inside/ndustpart
      PRINT *,'Gas particle mass, Dust particle mass ',
     &     gas_mass_inside/ngas_inside,dust_particle_mass

      ndust = 0
      DO j = 1, idim
         IF (iphase(j).EQ.-1 .OR. j.GT.npart) THEN
c
c--Disk distribution
c
 1237       i = MIN(INT(ran1(1)*npart)+1,npart)
            IF (iphase(i).NE.0) GOTO 1237

            IF (igeom.EQ.'s') THEN
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            ELSE
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
            ENDIF
            IF (SQRT(r2).LT.min_rplan.OR.SQRT(r2).GT.max_rplan)
     &           GOTO 1237

            ndust = ndust + 1

            xi = xyzmh(1,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)
            yi = xyzmh(2,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)
            zi = xyzmh(3,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)

            xyzmh(1,j) = xi
            xyzmh(2,j) = yi
            xyzmh(3,j) = zi
            xyzmh(4,j) = dust_particle_mass
            xyzmh(5,j) = xyzmh(5,i)*
     &           (REAL(ngas_inside)/ndustpart)**(1./3.)

            vxyzu(1,j) = vxyzu(1,i)
            vxyzu(2,j) = vxyzu(2,i)
            vxyzu(3,j) = vxyzu(3,i)
            vxyzu(4,j) = vxyzu(4,i)

            rho(j)     = dust_gas_ratio*rho(i)
            iphase(j)  = 11
            isteps(j)  = isteps(i)

            IF (ndust.EQ.ndustpart) THEN
               npart = j
               GOTO 100
            ENDIF
         ENDIF
      ENDDO

      PRINT *,'FAILED -- only managed to create ',ndust,' particles ',
     &     'before running out of space'
      STOP

 100  nlstacc = 0
      nlistinactive = 0
      nkill = 0
      n1 = 0
      nactive = 0
      DO i=1, npart
         IF (iphase(i).EQ.-1) THEN
            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

         ELSE
            n1 = i
            nactive = nactive + 1
         ENDIF
      ENDDO
      PRINT *, 'Final Numbers ',npart, nactive, n1
c
c--Dump new file
c
c
c--Zero time (but not quite, so that the timesteps are still used when
c     starting the code).
c
      IF (gt.NE.0. .AND. iok.EQ.'y') gt = 1.0E-20

      CALL wrinsph

      ifulldump = 0
      nfullstep = 1
      PRINT *,'writing dump file'

c
c--Transform into original frame of reference
c
      CALL chanref(1)
c
c--Write dump file
c
#ifdef MPI
      IF (iproc.EQ.0) THEN
         CALL file
      ENDIF
#endif

      CALL wdump_wrapper(idisk1)
c
c--End writing of full dump file
c-------------------------------
c

#ifdef MPI
      CLOSE(iread)
#endif

      PRINT 2000, file1
 2000 FORMAT ('file ', A10, 'has been created')
 
      RETURN
      END
