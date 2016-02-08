	PROGRAM TESTMU

	INCLUDE '../COMMONS/eostbl'
	INCLUDE '../COMMONS/units'
	REAL K

	REAL*4 rho

	umass = 1.0
	udist = 1.0
	udens = 1.0
	utime = 1.0
	uergg = 1.0
	uergcc = 1.0

	PRINT *,"Rho?"
	READ *,rho

	OPEN(UNIT=10,FILE='molmasstbl_lowdens',FORM='unformatted')

	DO I=1,mumxu
	   READ(10) (mutable(i,j), j=1,mumxrh)
	END DO


	OPEN(UNIT=8,FILE='gasttbl_lowdens',FORM='unformatted')
      	DO i=1, tgmxu
	   READ(8) (tgtable(i,j), j=1, tgmxrh)
           PRINT *,tgtable(i,2000),i
      	END DO
	print *,"Done"
      	CLOSE(8)

	PRINT *,"Read table"
	K=7.725
	DO I=1,mumxu+1000
	   K=K+0.005
	   u=10.0**(REAL(K))
	   y2=getcv(rho,u) 
	   y1=get1overmu(rho,u)
!           print *,u,rho,y1
	   WRITE(57,99001)  u/y2,u,y2,1.0/y1
99001	   FORMAT (4(1PE12.5,1X))
	END DO

	CLOSE(10)

	END
