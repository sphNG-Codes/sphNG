      SUBROUTINE prout(where)
c************************************************************
c                                                           *
c  This routine prints out all interesting quantities at    *
c     the present time                                      *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/angm'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/latti'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/trans'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/out1'
      INCLUDE 'COMMONS/out2'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/task'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/physeos'
      INCLUDE 'COMMONS/outmhd'

      CHARACTER*24 sentenc
      CHARACTER*7 where
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine prout')
c
c--Check for output type
c
      IF (where(1:6).EQ.'inform') THEN
c
c--Write integration output :
c
c--Number of cycle
c
c         IF (nbuild.EQ.1 .OR. MOD(ncount, nstep).EQ.0 .OR. 
         IF (MOD(ncount, nstep).EQ.0 .OR. 
     &                                        iptcreat.EQ.1) THEN
            WRITE (iprint, 99002, ERR=100) irec, MOD(ncount, nstep)+1
         ELSE
            WRITE (iprint, 99202, ERR=100) irec, MOD(ncount, nstep)+1
         ENDIF
99002    FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I2,
     &           '  W R I T T E N   O N   D I S K <------', /)
99202    FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I2,
     &           ' <------', /)
c
c--Write present time
c
         WRITE (iprint, 99003, ERR=100) gt
99003    FORMAT (1X, 'TIME  : ', 1PE16.10)
         tcomp = SQRT((3 * pi) / (32 * rhozero))
         WRITE (iprint, 99203, ERR=100) gt/tcomp
99203    FORMAT (1X, 'Free Fall Time  : ', 1PE16.10)
c
c--Energies + total angular momentum
c
         total = tkin + tgrav
         IF (igrp.NE.0) total = total + tterm + tmag + trad
         WRITE (iprint, 99004, ERR=100) total,tkin,trotz,trotx,tgrav,
     &                               tterm,tmag,trad,angto,totmom
99004    FORMAT (' General properties of system : ', /,
     &           ' total energy                       : ', 1PE14.5, /,
     &           ' kinetic energy                     : ', 1PE14.5, /,
     &           ' rotational energy around z         : ', 1PE14.5, /,
     &           ' rotational energy around x         : ', 1PE14.5, /,
     &           ' potential energy                   : ', 1PE14.5, /,
     &           ' internal energy                    : ', 1PE14.5, /,
     &           ' magnetic energy                    : ', 1PE14.5, /,
     &           ' radiation energy                   : ', 1PE14.5, /,
     &           ' total angular momentum             : ', 1PE14.5, /,
     &           ' total linear momentum              : ', 1PE14.5)
         
         IF (ABS(tgrav).GT.tiny) THEN
            alph = tterm / ABS(tgrav)
            betl = trotz / ABS(tgrav)
            betr = trotx / ABS(tgrav)
            ajeans = 1 / alph
         ELSE
            alph = 0.
            betl = 0.
            betr = 0.
            ajeans = 0.
         ENDIF
         WRITE (iprint, 88001, ERR=100) alph, betl, betr, ajeans
88001    FORMAT (' evolutionary energy parmeters : ', /,
     &           '                        alpha  : ', 1PE14.5, /,
     &           '       (z)       beta parallel : ', 1PE14.5, /,
     &           '       (x)       beta perpend  : ', 1PE14.5, /,
     &           '                 Jeans number  : ', 1PE14.5)
c
c--Escapors
c
         WRITE (iprint, 99005, ERR=100) escap
99005    FORMAT (' escapors (total mass)         : ', 1PE14.5)
c
c--Accretion
c
         WRITE (iprint, 88002, ERR=100) nactive
88002    FORMAT (' number of active particles    : ', I8)
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9 .OR. ibound.EQ.100) THEN
            WRITE (iprint, 88003, ERR=100) nreassign
88003       FORMAT (' number of reassigned part.    : ', I8)
            WRITE (iprint, 88004, ERR=100) naccrete
88004       FORMAT (' number of accreted part.      : ', I8)
            WRITE (iprint, 88005, ERR=100) nkill
88005       FORMAT (' number of killed particles    : ', I8)
         ENDIF
         IF (iaccevol.EQ.'v' .OR. iaccevol.EQ.'s') THEN
            WRITE (iprint, 88006, ERR=100)
88006       FORMAT (' point mass data:')
            DO i = 1, nptmass
               ipt = listpm(i)
               WRITE (iprint, 88007, ERR=100) i, iphase(ipt), 
     &              xyzmh(5,ipt)
88007          FORMAT ('  point mass: ',I2,' type: ',I1,
     &                 ' hacc: ',1PE12.3)
            END DO
            WRITE (iprint,88008) rmax
88008       FORMAT (' new boundary radius: ',1PE12.5)
         ENDIF
c
c--Magnetic field parameters	(computed in mhdparams)
c
         IF (imhd.EQ.idim) THEN    
	    WRITE (iprint, 66001, ERR=100) betamhdmin,betamhdmax,
     &           betamhdav,divBmax,divBav,curlBmax,curlBav,div2curlBmax,
     &           div2curlBav,omegamhdmax,omegamhdav,omegtol,fracdivBok,
     &           fluxtot,crosshel
66001       FORMAT (/,' Magnetic field parameters : ', /,
     &   '          plasma beta  min : ', 1PE14.5,'  max : ', 1PE14.5,
     &   ' mean : ', 1PE14.5, /,     
     &   '                 div B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '                curl B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '      (div B)/(curl B) max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '              divB*h/B max : ', 1PE14.5,
     &   ' mean : ', 1PE14.5,' frac < ',1PE9.2,' : ',0PF6.2,'%',//,
     &   ' total magnetic flux  (int nabla.B dV) : ', 1PE14.5,/,
     &   ' total cross helicity (int v.B dV)     : ', 1PE14.5)
         ENDIF
c
c--Object no 1
c
         WRITE (iprint, 99006, ERR=100) n1, cmx1, cmy1, cmz1, vcmx1,
     &        vcmy1, vcmz1, hmi1, hma1, dmax1,
     &        zmax1, romean1, romax1, rocen1,
     &        valphamin1, valphamax1, tgmean1, tgmax1, tgcen1
99006    FORMAT (/, ' Object number 1 (', I8, ' particles ) : ', /,
     &           ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,
     &           '  z :', 1PE14.5, /, ' velocity cm    vx :', 1PE14.5,
     &           ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,
     &           ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' max. dist. cm   r :', 1PE14.5, '  z :', 1PE14.5, /,
     &           ' density      mean :', 1PE14.5, ' max:', 1PE14.5, 
     &           ' cen:', 1PE14.5, /,
     &           ' visc. switch  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' temperature  mean :', 1PE14.5, ' max:', 1PE14.5,
     &           ' cen:', 1PE14.5, /
     &           )
         IF (imhd.EQ.idim) THEN
            WRITE (iprint, 99206) Bmin, Bmean, Bmax
99206       FORMAT (' mag field    min  :', 1PE14.5, 'mean:', 1PE14.5,
     &              ' max:', 1PE14.5,/)         
         ENDIF
         IF (encal.EQ.'r') THEN
            WRITE (iprint, 99106, ERR=100) trmean1, trmax1, trcen1
99106       FORMAT (' temper. rad. mean :', 1PE14.5, ' max:', 
     &           1PE14.5, ' cen:', 1PE14.5, /)
         ENDIF

c
c--Object no 2
c
         IF (n2.NE.0) THEN
            WRITE (iprint, 99007, ERR=100) n2, cmx2, cmy2, cmz2, vcmx2,
     &             vcmy2, vcmz2, hmi2, hma2, dmax2, zmax2, romean2,
     &             romax2, rocen2, valphamin2, valphamax2
99007       FORMAT (/, ' Object number 2 (', I8, ' particles ) : ', /,
     &              ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,
     &              '  z :', 1PE14.5, /, ' velocity cm    vx :',
     &              1PE14.5, ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,
     &              ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5,
     &              /, ' max. dist. cm   r :', 1PE14.5, '  z :',
     &              1PE14.5, /, ' density      mean :', 1PE14.5,
     &              ' max:', 1PE14.5,' cen:', 1PE14.5, /,
     &           ' visc. switch  min :', 1PE14.5, ' min:', 1PE14.5, /)
         ENDIF
c
c--Write transfer output
c
      ELSEIF (where(1:5).EQ.'trans') THEN
         WRITE (iprint, 99008) ibegin, iend, file1, file2
99008    FORMAT (//, ' dumps no ', I4, ' to ', I4, ' copied from file ',
     &           A7, ' to file ', A7, /)
         IF (ichang.EQ.0) WRITE (iprint, 99009)
99009    FORMAT (' no change done. ')
         IF (ichang.EQ.1) WRITE (iprint, 99010) frac, energc
99010    FORMAT (' the inner ', F5.3,
     &           ' of the total number of particles',
     &           ' have had their internal energy multiplied by: ',
     &           1PE12.5)
         IF (ichang.EQ.2) THEN
            WRITE (iprint, 99010) frac, energc
            WRITE (iprint, 99011) vexpan*udist/utime, rnorm*udist
99011       FORMAT (' a homologous expansion of ', 1PE12.5, ' cm/s',
     &              ' at ', 1PE12.5, ' cm has been added')
         ENDIF
c
c--Write output for initialisation
c
      ELSEIF (where(1:6).EQ.'newrun') THEN
         IF (varsta.EQ.'entropy') sentenc = 'specific entropy'
         IF (varsta.EQ.'intener') sentenc = 'specific internal energy' 

         IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
            WRITE (iprint, 99012)
99012       FORMAT (' New initial conditions defined from scratch')
            WRITE (iprint, 99013) idist
99013       FORMAT (
     &           ' Particles distributed according to distribution no :'
     &           , I2)
            WRITE (iprint, 99014) sentenc, thermal
99014       FORMAT (' Value of ', A24, ' : ', 1PE12.5)

            IF (encal.EQ.'i') WRITE (iprint, 99100)
99100       FORMAT (' Isothermal equation of state')
            IF (encal.EQ.'t') WRITE (iprint, 99777)
99777       FORMAT (' Thermal instability equation of state')
            IF (encal.EQ.'a') WRITE (iprint, 99103)
99103       FORMAT (' Adiabatic equation of state')
            IF (encal.EQ.'p') WRITE (iprint, 99101) gamma
99101       FORMAT (' Polytropic equation of state gamma is ',1PE12.5)
            IF (encal.EQ.'v') WRITE (iprint, 99102)
99102       FORMAT (' Gamma variable equation of state')
            IF (encal.EQ.'x') WRITE (iprint, 99104)
99104       FORMAT (' Physical equation of state')
            IF (encal.EQ.'r') WRITE (iprint, 99105)
99105       FORMAT (' Implicit Radiation/Adiabatic equation of state')

            IF (idist.GE.1) WRITE (iprint, 99015) xlmax, xlmin, ylmax,
     &                               ylmin, zlmax, zlmin
99015       FORMAT (/,' Particles are in volume given by :', /,
     &              ' xmax : ', 1PE12.5, ' xmin : ', 1PE12.5, /,
     &              ' ymax : ', 1PE12.5, ' ymin : ', 1PE12.5, /,
     &              ' zmax : ', 1PE12.5, ' zmin : ', 1PE12.5)
         ENDIF

         IF (what(1:5).EQ.'exist' .OR. what.EQ.'e') THEN
            WRITE (iprint, 99016) file1, file2
99016       FORMAT (//, ' new initial conditions made from file ', A7,
     &              ' and ', A7)
            WRITE (iprint, 99017) n1, xx1, yy1, zz1, vvx1, vvy1, vvz1
99017       FORMAT (/, ' first   object is made of ', I4,
     &              ' particles : ', /, '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99018) n2, xx2, yy2, zz2, vvx2, vvy2, vvz2
99018       FORMAT (/, ' second  object is made of ', I4,
     &              ' particles : ', /, '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99019) angto
99019       FORMAT (/, ' total angular momentum of the system : ',
     &              1PE12.5)
         ENDIF
      ENDIF
      GOTO 200
c
c--Handle errors during writing
c
 100  where = 'prout'
      CALL error(where, 1)

 200  CALL FLUSH (iprint)
      RETURN
      END
