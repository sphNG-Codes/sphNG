      subroutine read_dragtable

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/logun'
      
      INTEGER i, j, k
      INTEGER ires_r, ires_p, ires_t
      REAL trmin, trmax, tdr
      REAL tpmin, tpmax, tdp
      REAL ttmin, ttmax, tdt
      REAL dragpmass, draghmass
      REAL tol

      tol = 1.0E-2
c
c--Read in the resolution and check it matches current code
c
      IF (numplanet.NE.1) THEN
         WRITE (*,*) 'ERROR - numplanet must be 1 for dragtable'
         CALL quit(0)
      ENDIF
      
      READ (idragfile, ERR=100) ires_r, ires_p, ires_t

      IF (ires_r.NE.idragresr .OR. ires_p.NE.idragresp .OR.
     &     ires_t.NE.idragrest) THEN
 111     FORMAT('Drag table resolution being read in does not ',
     &        'match resolution set in radtran3',/)
 112     FORMAT('(code, table) ', A2, ' = ', I6, 1X, I6)
         WRITE (iprint,111)
         WRITE (iprint,112) 'ir', idragresr, ires_r
         WRITE (iprint,112) 'ip', idragresp, ires_p
         WRITE (iprint,112) 'it', idragrest, ires_t
         CALL quit(0) 
      ENDIF
c
c--Read in spatial scales and check they match the current code
c
      READ (idragfile, ERR=100) trmin, trmax, tdr, 
     &     tpmin, tpmax, tdp, ttmin, ttmax, tdt

      IF (abs(trmin-dragtrmin).GT.tol .OR.
     &     abs(trmax-dragtrmax).GT.tol .OR.
     &     abs(tdr-dragtdr).GT.tol) THEN
 113     FORMAT ('Drag table: Radial scale changed')
 114     FORMAT ('(read in) min, max, delta ', 3(1PE12.5,1X))
 115     FORMAT ('(new    ) min, max, delta ', 3(1PE12.5,1X))
         WRITE (iprint,113)
         WRITE (iprint,114) trmin, trmax, tdr
         WRITE (iprint,115) dragtrmin, dragtrmax, dragtdr
         CALL quit(0)
      ENDIF
c
c--Read in planet mass, and check it matches current state
c
      READ (idragfile, ERR=100) dragpmass, draghmass

      IF (abs(dragpmass-planetmass(1)).GT.tiny) THEN
         IF (abs(dragpmass-planetmass(1)).LT.1.0E-6) THEN
 116        FORMAT ('Planet mass does not match in dragfile')
 117        FORMAT ('(read in here   ) mass ', 1PE15.8)
 118        FORMAT ('(read in options) mass ', 1PE15.8)
            WRITE (iprint, 116)
            WRITE (iprint,117) dragpmass
            WRITE (iprint,118) planetmass(1)
            WRITE (iprint, *) 'Drag file mass overwrites ifile'
            planetmass(1) = dragpmass
         ELSE
 119        FORMAT ('ERROR : Mass difference too large.')
            WRITE (iprint, 116)
            WRITE (iprint,117) dragpmass
            WRITE (iprint,118) planetmass(1)
            WRITE (iprint, 119)
            CALL quit(0) 
         ENDIF
      ENDIF

c
c--If all checks successful, repopulate the table from file
c
      do i = 1, idragresr
         do j = 1, idragresp
            do k = 1, idragrest
               READ (idragfile, ERR=100) dragenergy(i,j,k)
            enddo
         enddo
      enddo

      CLOSE (idragfile)

      RETURN
      
 100  print *, 'An error has occurred reading drag file'
      CALL quit(0) 

      end subroutine read_dragtable

