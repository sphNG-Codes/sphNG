      SUBROUTINE lunit
c************************************************************
c                                                           *
c  This routine attributes logical unit numbers.            *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/binfile'

c      iprint = 9
      iprint = 6
      iterm = 50 
      idisk1 = 11
      idisk2 = 12
      idisk3 = 13
      iptprint = 14
      iaccpr = 15
      ikillpr = 16
      ireasspr = 17
      inotify = 18

      imaxrec = 60*idim

      RETURN
      END
