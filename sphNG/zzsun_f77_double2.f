      SUBROUTINE zzsun2

c      USE f90_unix

c************************************************************
c                                                           *
c  This routine contains all machine dependant routines     *
c                                                           *
c************************************************************

c************************************************************
c
c  NOTE: Must be compiled separately NOT using -dbl
c
c************************************************************
      REAL*4 tarray(2)
      REAL*8 tused
c
c--Get time used since begining
c
      ENTRY getused(tused)
      tused = etime(tarray)
      RETURN

      END
