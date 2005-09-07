      SUBROUTINE wdump(idisk1)
c************************************************************
c                                                           *
c  This routine writes a dump on disk                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/crap'

      CHARACTER*7 where
      CHARACTER*100 fileident
      INTEGER*4 int1, int2
      INTEGER*8 number8
      DIMENSION nums(8)
      REAL dummy(3),hfact

      DATA where/'wdump'/
c
c--Write
c
      irec = irec + 1
      iresort = iresort + 1
      ifulldump = ifulldump + 1
c
c--Write dump file in SPMHD format
c-----------------------------------


      DO j=1,3
         dummy(j) = 0.
      ENDDO
      hfact = 1.2
      WRITE (idisk1, ERR=100) gt,npart,npart,gamma,hfact,3,3,10,1,0,
     &                        0,0,0,dummy(1:3),dummy(1:3)
c--Default real
      WRITE (idisk1, ERR=100) ((xyzmh(j,isort(i)),j=1,3), i=1,npart)
      WRITE (idisk1, ERR=100) ((vxyzu(j,isort(i)),j=1,3), i=1,npart)
      WRITE (idisk1, ERR=100) (xyzmh(5,isort(i)), i=1,npart)
      WRITE (idisk1, ERR=100) (real(rho(isort(i))), i=1, npart)
      WRITE (idisk1, ERR=100) (vxyzu(4,isort(i)), i=1,npart)
      WRITE (idisk1, ERR=100) (xyzmh(4,isort(i)), i=1,npart)

      ENDFILE idisk1
      BACKSPACE idisk1

      ipos = irec

      RETURN
c
c--An error as occured while writing
c
 100  CALL error(where, 1)

      RETURN
      END
