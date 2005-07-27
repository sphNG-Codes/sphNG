      SUBROUTINE mainop
c************************************************************
c                                                           *
c  This subroutine attributes first logical units and then  *
c     determines the main task of current run               *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/actio'
c
c--Define logical units first
c
      CALL lunit
c
c--Get main option of job from command line argument
c
      CALL getcom(job, inname)

      IF (job(1:9).NE.'evolution') iprint = 6 

      RETURN
      END
