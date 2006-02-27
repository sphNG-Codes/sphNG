		PROGRAM TESTKAPPA

		INCLUDE '../COMMONS/tgtbl'
		INCLUDE '../COMMONS/optbl'

      OPEN(UNIT=8,FILE='/home/scw/tables/opacitytbl')
      	DO i=1, opmxtg
      		READ(8,*) (optable(i,j), j=1, opmxrh)
      	ENDDO
      CLOSE(8)

      OPEN(UNIT=8,FILE='/home/scw/tables/gasttbl',FORM='unformatted')
      	DO i=1, tgmxu
      		READ(8) (tgtable(i,j), j=1, tgmxrh)
      	ENDDO
      CLOSE(8)
 
		PRINT *,"Enter rho:"
		READ *,rho
		startu = 7.725
		DO K=1,1300
			startu = startu + 0.005
			u = 10.0**DBLE(startu)
			cv = getcv(rho,u)
			rk = getkappa(u,cv,rho)
			WRITE(30,*) u/cv,rk
		ENDDO

		END
