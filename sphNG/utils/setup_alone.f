      PROGRAM setup_alone

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE '../COMMONS/savernd'
      INCLUDE '../COMMONS/rbnd'
      INCLUDE '../COMMONS/sort'

      nlistinactive = 0

c
c--Initialize physical, mathematical and astronomical constants
c
      CALL constan
c
c--Define code's units
c
      CALL unit
c
c--Initialise starting date and time
c  (some versions of the zz file require this)
c
      CALL getdat(iday, imon, iyear)
c
c--Initialise random number generator arrays
c
      iyr = 0
      DO i = 1, NTAB
         iv(i) = 0
      END DO
c
c--Set up new run
c
c      WRITE (*, 99006)
c99006 FORMAT (' name of binary file (7 char. max)')
c      READ (*, 99005) file1
c      WRITE (*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec,' ',file1
c      OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
c     &     RECL=imaxrec)

      CALL setpart
      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO
      nfullstep = 1

      file1 = 'XX00000'
      CALL file
      OPEN (11,FILE=file1,FORM='unformatted')

      CALL wdump(11)

      CALL quit
      END


      SUBROUTINE quit
      INCLUDE 'idim'

      STOP
      END

      SUBROUTINE endrun
      CALL quit
      END

      SUBROUTINE ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

      STOP
      END

      SUBROUTINE insulate(ival, ntot, npart, dumxyzmh, f1vxyzu)

      STOP
      END

      SUBROUTINE labrun
      RETURN
      END

      SUBROUTINE preset(ihopin)
      RETURN
      END

      SUBROUTINE hcalc
      RETURN
      END

      SUBROUTINE getcv(rho2,u2)
      RETURN
      END

      SUBROUTINE get1overmu(rho2,u2)
      RETURN
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
 1              v1=2.*ran1(idum)-1.
         v2=2.*ran1(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END
