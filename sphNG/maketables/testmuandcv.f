		PROGRAM TESTMU

		INCLUDE '../COMMONS/tgtbl'
		INCLUDE '../COMMONS/mutbl'
		REAL K

		PRINT *,"Rho?"
		READ *,rho

      OPEN(UNIT=10,FILE='molmasstbl',FORM='unformatted')

		DO I=1,mumxu
          READ(10) (mutable(i,j), j=1,mumxrh)
      END DO


		OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
      	DO i=1, tgmxu
      		READ(8) (tgtable(i,j), j=1, tgmxrh)
c				PRINT *,tgtable(i,2000),i
      	ENDDO
		print *,"Done"
      	CLOSE(8)

	   PRINT *,"Read table"
      K=7.725
      DO I=1,mumxu+1000
         K=K+0.005
         u=10.0**(REAL(K))
		   y2=getcv(rho,u) 
         y1=get1overmu(rho,u)
!			print *,u,rho,y1
         WRITE(57,*)  u/y2,1.0/y1,y2
      ENDDO

	

      CLOSE(10)

 		END
