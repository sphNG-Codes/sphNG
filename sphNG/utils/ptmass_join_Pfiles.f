      PROGRAM join_ptmass
c
c--Joins together 2 ptmass P-files into a single one
c
      PARAMETER (nfilemax = 1000)
      
      CHARACTER*21 infile(nfilemax), outfile
      CHARACTER*1 ians

      INTEGER*8 j, imerged1, imerged2

      INTEGER*8 listsink_old(1000),listsink_new(1000)
      DIMENSION pmass_old(1000),pmass_new(1000)
      DIMENSION tform(100000000)

      WRITE (*,*) 'Enter number of P-files'
      READ (*,*) nfiles

      IF (nfiles.GT.nfilemax) THEN
         WRITE (*,*) 'ERROR - too many files'
         STOP
      ENDIF

      WRITE (*,*) 'Enter filenames'
      DO i = 1, nfiles
         READ (*,99001) infile(i)
99001    FORMAT(A21)
      END DO

      WRITE (*,*) 'Enter output filename:' 
      READ (*,99001) outfile

      WRITE (*,*) 'Is it an accretion run?'
      READ (*,99004) ians
99004 FORMAT(A1)

      ilow = 15
      ihi  = 16
      ifile = 3
      timelast = 0.
      nold = 0

      OPEN (ilow, FILE=infile(1), STATUS='unknown', FORM='unformatted')
      OPEN (ihi, FILE=infile(2), STATUS='unknown', FORM='unformatted')
      OPEN (17, FILE=outfile, STATUS='unknown', FORM='unformatted')

 50   READ (ihi, END=500,ERR=500) fftime2, time2, nptmass2
c      write (*,*) ihi,fftime2

 100  READ (ilow, END=200,ERR=200) fftime, time, nptmass
c      write (*,*) ilow,fftime

      IF (time.LT.time2 .AND. time.GT.timelast) THEN
         IF (nptmass.GT.0 .OR. (nptmass.LT.0 .AND. 
     &        ABS(nptmass).LT.nold+10)) 
     &        WRITE (17) fftime, time, nptmass
         timelast = time
         IF (nptmass.LT.0) THEN
            READ (ilow, END=200, ERR=200) imerged1, imerged2
            WRITE (*,*) 'MERGER ',imerged1, imerged2, nptmass, nold, 
     &           nnew
            IF (ABS(nptmass).GT.nold+10) THEN
               WRITE (*,*) 'ERROR in FILE -- ignoring merger'
               imerged1 = 0
               imerged2 = 0
            ELSE
               WRITE (17) imerged1, imerged2
            ENDIF
            GOTO 100
         ELSE
            nnew = nptmass
            DO i = 1, nptmass
               IF (ians.EQ.'y' .OR. ians.EQ.'Y') THEN
                  READ (ilow, END=200, ERR=200) j, x, y, z, 
     &                 vx, vy, vz, pmass, rho, 
     &                 nactotal, ptmassinner,spinx,spiny,spinz,angaddx,
     &                 angaddy, angaddz, naccrete, anglostx, anglosty,
     &                 anglostz, nkill
                  WRITE (17) j, x, y, z, vx, vy, vz,pmass,rho,nactotal,
     &                 ptmassinner, spinx, spiny, spinz, angaddx, 
     &                 angaddy, angaddz, naccrete, anglostx, anglosty,
     &                 anglostz, nkill
               ELSE
                  READ (ilow, END=200,  ERR=200) j, x, y, z, 
     &                 vx, vy, vz, pmass, rho, 
     &                 nactotal, ptmassinner,spinx,spiny,spinz,angaddx, 
     &                 angaddy, angaddz, naccrete

                  listsink_new(i) = j
                  IF (tform(listsink_new(i)).EQ.0.) THEN
                     tform(listsink_new(i)) = fftime
                     print *,'NEW ',nnew,nold,fftime,j
                  ENDIF
                  pmass_new(i) = pmass

                  IF (j.EQ.imerged1 .OR. j.EQ.imerged2) THEN
                     WRITE (*,*) 'mass ',j,pmass,nold,nnew
                     DO kk = 1, nold
                        IF (listsink_old(kk).EQ.imerged1) 
     &                       print *,'Old mass1 ',pmass_old(kk),
     &                       tform(listsink_old(kk)),fftime
                        IF (listsink_old(kk).EQ.imerged2) 
     &                       print *,'Old mass2 ',pmass_old(kk),
     &                       tform(listsink_old(kk)),fftime
                     END DO
                     IF (j.EQ.imerged1) imerged1 = 0
                     IF (j.EQ.imerged2) imerged2 = 0
                  ENDIF                     

                  WRITE (17) j, x, y, z, vx,vy,vz,pmass,rho,nactotal,
     &                 ptmassinner, spinx, spiny, spinz, angaddx, 
     &                 angaddy, angaddz, naccrete
               ENDIF
            END DO

            kk_keep = 0
            IF (nnew.LT.nold) THEN
               DO kk = 1, nold
                  DO jj = 1, nnew
                     IF (listsink_old(kk).EQ.listsink_new(jj)) GOTO 1000
                  END DO
                  WRITE (*,*) 'Couldnt find1 ',kk,pmass_old(kk),
     &                 tform(listsink_old(kk)),listsink_old(kk)
                  kk_keep = kk
 1000             CONTINUE
               END DO
            ENDIF

            DO jj = 1, nnew
               listsink_old(jj) = listsink_new(jj)
               pmass_old(jj) = pmass_new(jj)
            END DO
            nold = nnew

         ENDIF
         GOTO 100
      ELSE
         READ (ilow, END=200, ERR=200)
         GOTO 100
      ENDIF
         
 200  WRITE (17) fftime2, time2, nptmass2
      timelast = time2
      IF (nptmass2.LT.0) THEN
         READ (ihi, END=500, ERR=500) imerged1, imerged2
         WRITE (17) imerged1, imerged2
         WRITE (*,*) 'MERGER2 ',imerged1, imerged2
      ELSE
         nnew = nptmass2
         DO i = 1, nptmass2
            IF (ians.EQ.'y' .OR. ians.EQ.'Y') THEN
               READ (ihi, END=500, ERR=500) j, x, y, z, 
     &              vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx,
     &              angaddy, angaddz, naccrete, anglostx, anglosty,
     &              anglostz, nkill
               WRITE (17) j, x, y, z, vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete, anglostx, anglosty,
     &              anglostz, nkill
            ELSE
               READ (ihi, END=500, ERR=500) j, x, y, z, 
     &              vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete

               listsink_new(i) = j
               IF (tform(listsink_new(i)).EQ.0.) THEN
                  tform(listsink_new(i)) = fftime
                  print *,'NEW2 ',nnew,nold,fftime,j
               ENDIF

               pmass_new(i) = pmass

               IF (j.EQ.imerged1 .OR. j.EQ.imerged2) THEN
                  WRITE (*,*) 'mass2 ',j,pmass,nold,nnew
                  IF (j.EQ.imerged1) imerged1 = 0
                  IF (j.EQ.imerged2) imerged2 = 0
                  DO kk = 1, nold
                     IF (listsink_old(kk).EQ.j) 
     &                    print *,'Old mass2 ',pmass_old(kk)
                  END DO
               ENDIF

               WRITE (17) j, x, y, z, vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete
            ENDIF
         END DO

         IF (nnew.LT.nold) THEN
            DO kk = 1, nold
               DO jj = 1, nnew
                  IF (listsink_old(kk).EQ.listsink_new(jj)) GOTO 1001
               END DO
               WRITE (*,*) 'Couldnt find2 ',kk,pmass_old(kk)
 1001          CONTINUE
            END DO
         ENDIF

         DO jj = 1, nnew
            listsink_old(jj) = listsink_new(jj)
            pmass_old(jj) = pmass_new(jj)
         END DO
         nold = nnew

      ENDIF

      itemp = ilow
      ilow = ihi
      ihi = itemp
      CLOSE (ihi)
      IF (ifile.LE.nfiles) THEN
        WRITE (*,*) 'Open file ',ifile,infile(ifile)
        OPEN (ihi, FILE=infile(ifile), STATUS='unknown', 
     &                                     FORM='unformatted')
        ifile = ifile + 1
        GOTO 50
      ENDIF

 300  READ (ilow, END=400) fftime, time, nptmass
c      write (*,*) ilow,fftime

      IF (time.GT.timelast) THEN
         WRITE (17) fftime, time, nptmass
         timelast = time
         IF (nptmass.LT.0) THEN
            READ (ilow, END=400, ERR=400) imerged1, imerged2
            WRITE (17) imerged1, imerged2
            WRITE (*,*) 'MERGER3 ',imerged1,imerged2
         ELSE
            nnew = nptmass
            DO i = 1, nptmass
            IF (ians.EQ.'y' .OR. ians.EQ.'Y') THEN
               READ (ilow, END=400, ERR=400) j, x, y, z, 
     &              vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx,
     &              angaddy, angaddz, naccrete, anglostx, anglosty,
     &              anglostz, nkill
               WRITE (17) j, x, y, z, vx, vy, vz, pmass, rho, nactotal,
     &              ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete, anglostx, anglosty,
     &              anglostz, nkill
            ELSE
               READ (ilow, END=400, ERR=400) j, x, y, z, 
     &              vx, vy, vz, pmass, rho, 
     &              nactotal, ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete

               listsink_new(i) = j
               IF (tform(listsink_new(i)).EQ.0.) THEN
                  tform(listsink_new(i)) = fftime
                  print *,'NEW3 ',nnew,nold,fftime,j
               ENDIF
               pmass_new(i) = pmass

               IF (j.EQ.imerged1 .OR. j.EQ.imerged2) THEN
                  WRITE (*,*) 'mass3 ',j,pmass,nold,nnew
                  DO kk = 1, nold
                     IF (listsink_old(kk).EQ.imerged1) 
     &                    print *,'Old mass3-1 ',pmass_old(kk)
                     IF (listsink_old(kk).EQ.imerged2) 
     &                    print *,'Old mass3-2 ',pmass_old(kk)
                  END DO
                  IF (j.EQ.imerged1) imerged1 = 0
                  IF (j.EQ.imerged2) imerged2 = 0
               ENDIF

               WRITE (17) j, x, y, z, vx, vy, vz, pmass, rho, nactotal,
     &              ptmassinner, spinx, spiny, spinz, angaddx, 
     &              angaddy, angaddz, naccrete
            ENDIF
            END DO

            IF (nnew.LT.nold) THEN
               DO kk = 1, nold
                  DO jj = 1, nnew
                     IF (listsink_old(kk).EQ.listsink_new(jj)) GOTO 1002
                  END DO
                  WRITE (*,*) 'Couldnt find3 ',kk,pmass_old(kk),
     &                 tform(listsink_old(kk)),listsink_old(kk)
 1002             CONTINUE
               END DO
            ENDIF

            DO jj = 1, nnew
               listsink_old(jj) = listsink_new(jj)
               pmass_old(jj) = pmass_new(jj)
            END DO
            nold = nnew


         ENDIF
         GOTO 300
      ELSE
         READ (ilow, END=400, ERR=400)
         GOTO 300
      ENDIF

 400  CONTINUE

      CLOSE (ilow)
      CLOSE (17)

      STOP

 500  CONTINUE

      WRITE (*,*) 'Error - file ',ihi,' only has one record'
      WRITE (*,*) 'Leave it out of the list of filenames'
      STOP

      END
