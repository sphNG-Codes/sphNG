	PROGRAM TESTMU

c--Converts the binary tables into ascii format. Requires
c  the binary tables to have been created first. The data
c  is put into fort.30, whilst the density axis is given
c  in fort.31, and the temperature axis in fort.32.

c--Written 21/09/2010 by Ben Ayliffe
	
	INTEGER i, j
	
        INTEGER tgmxu,tgmxrh,status,nrho, ntm
	REAL*8 lrho, ltm
	PARAMETER (tgmxu=1799)
	PARAMETER (tgmxrh=4601)
	REAL*8, ALLOCATABLE :: tgtable(:,:)
	
	ALLOCATE( tgtable(tgmxu,tgmxrh) , STAT=status)
	IF (status /= 0) print *, 'Failed-Allocate'
	
	OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
      	DO i=1, tgmxu
	   READ(8) (tgtable(i,j), j=1, tgmxrh)
	   WRITE(30,123) (tgtable(i,j), j=1, tgmxrh)
      	ENDDO

      	CLOSE(8)
	DEALLOCATE( tgtable, STAT=status)
	IF (status /= 0) print *, 'Failed-Deallocate'
	
 123	FORMAT (4601(1PE12.5,2X))

	DO nrho=-20000,3000,5
	   lrho=nrho/1000.0
	   WRITE(31,234) lrho
	ENDDO

	DO ntm=0,9000,5
	   ltm=ntm/1000.0
	   WRITE(32,234) ltm
	ENDDO

 234	FORMAT (1PE12.5)

	END
