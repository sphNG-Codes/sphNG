	PROGRAM TESTKAPPA

	INCLUDE '../COMMONS/tgtbl'
	INCLUDE '../COMMONS/optbl'
	INCLUDE '../COMMONS/units'

	REAL*4 rho

	umass = 1.0
        udist = 1.0
        udens = 1.0
        utime = 1.0
        uergg = 1.0
        uergcc = 1.0

	OPEN(UNIT=8,FILE='/home/mrbate/tables/opacitytbl')
      	DO i=1, opmxtg
	   READ(8,*) (optable(i,j), j=1, opmxrh)
      	ENDDO
	CLOSE(8)

	OPEN(UNIT=8,FILE='/home/mrbate/tables/gasttbl',
     &   FORM='unformatted')
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
