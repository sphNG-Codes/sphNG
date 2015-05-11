      SUBROUTINE makedust
c************************************************************
c                                                           *
c  This subroutine reads in an existing dump file and       *
c     brings dead particles to life as dust particles.      *
c     Originally used to add dust to a settled disc.        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

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

      CHARACTER*40 ifile, ofile
      
1000  FORMAT (A40)
1001  FORMAT (A1)
1002  FORMAT (A11)
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF
      PRINT *, 'Name of binary dump file to add dust to ?'
      READ (*, 1000) ifile

      PRINT *, 'Name of output binary dump file ?'
      READ (*, 1000) ofile

      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'Name of corresponding ASCII ifile ?'
      READ (*, 1002) inname

      imax = 1073741824
      imaxstep = imax/2

      nactive = 0
      DO i = 1, npart
         isort(i) = i
         iorig(i) = i
      END DO

      OPEN (UNIT = 8, FILE = ifile, FORM = 'unformatted')
      PRINT *, 'reading file ', ifile
c
c--Read dump file
c
      CALL options

      file1 = ofile
      print *, 'Name0=', file1

      CALL rdump_wrapper(8,ichkl,1)
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
      READ (*,*) min_rplan, max_rplan

      gas_mass_inside = 0.
      ngas_inside = 0.
      DO i = 1, npart
         IF (iphase(i).EQ.0) THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
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
      READ (*,*) ndustpart
      IF (nactive+ndustpart.GT.npart) THEN
         PRINT *,'ERROR - number must be <= to ',npart-nactive
         GOTO 50
      ENDIF

      PRINT *,'Enter desired dust to gas ratio (normally < 1)'
      READ (*,*) dust_gas_ratio
      dust_particle_mass = dust_gas_ratio*gas_mass_inside/ndustpart
      PRINT *,'Gas particle mass, Dust particle mass ',
     &     gas_mass_inside/ngas_inside,dust_particle_mass

      ndust = 0
      DO j = 1, npart
         IF (iphase(j).EQ.-1) THEN
c
c--Disk distribution
c
 1237       i = MIN(INT(ran1(1)*npart)+1,npart)
            IF (iphase(i).NE.0) GOTO 1237

            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
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
            xyzmh(5,j) = xyzmh(5,i)

            vxyzu(1,j) = vxyzu(1,i)
            vxyzu(2,j) = vxyzu(2,i)
            vxyzu(3,j) = vxyzu(3,i)
            vxyzu(4,j) = vxyzu(4,i)

            rho(j)     = dust_gas_ratio*rho(i)
            iphase(j)  = 11
            isteps(j)  = isteps(i)

            IF (ndust.EQ.ndustpart) GOTO 100
         ENDIF
      ENDDO

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
      CALL wrinsph

      ifulldump = 0
      nfullstep = 1
      PRINT *,'writing dump file'
      CALL wdump_wrapper(7)
c
c--End writing of full dump file
c-------------------------------
c
      CLOSE (8)

      CLOSE (7)

      PRINT 2000, ofile
 2000 FORMAT ('file ', A10, 'has been created')
 
      RETURN
      END
