      SUBROUTINE chekopt
*************************************************************
*                                                           *
*   This subroutine checks for possible incompatibilities   * 
*      between the different options selected; it also      *
*      warns the user of possible problems related to the   *
*      present limitations of the code                      *
*                                                           *
*************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'

      CHARACTER*7 where

      DATA where/'chekopt'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine chekopt')

      IF (iener.EQ.1 .AND. damp.NE.0.) CALL error(where, 2)

      IF (encal.EQ.'i' .AND. (iener.EQ.1 .OR. ichoc.EQ.1)) 
     &                                   CALL error(where,3)

      IF ((encal.EQ.'a' .OR. encal.EQ.'r') .AND. 
     &     (iener.EQ.0 .OR. ichoc.EQ.0)) CALL error(where,4)

      IF ((encal.EQ.'p' .OR. encal.EQ.'v' .OR. encal.EQ.'x') .AND.  
     &     (iener.EQ.1 .OR. ichoc.EQ.1)) CALL error(where,5)

      RETURN
      END
