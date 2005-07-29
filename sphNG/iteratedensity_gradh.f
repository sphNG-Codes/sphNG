      SUBROUTINE iterate_density(jneigh,npart,xyzmh,vxyzu,
     &           nlst_in,nlst_end,list,itime)
c************************************************************
c                                                           *
c  Subroutine to calculate the density and smoothing        *
c  length self consistently as described in Price(2004)     *
c                                                           *
c  Uses Newton-Raphson iterations                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'
      INCLUDE 'COMMONS/physcon'

      INTEGER maxits,nneighmax,neides
      PARAMETER(maxitsnr=100,maxits=100)
      PARAMETER(nneighmax=300)
      REAL tol,third,hfact
      PARAMETER(hfact=1.2) !!!,neides=INT(32./3.*pi*hfact**3))
      PARAMETER(tol=1.e-3,third=1./3.)

      DIMENSION xyzmh(5,idim),vxyzu(4,idim),list(idim)
      DIMENSION newlist(idim), iredolist(idim), ifakelist(idim)

      INTEGER jneigh,ncalc
      REAL fsx,fsy,fsz,epot
      REAL pmassi,rhoi,omegai
      DIMENSION hminbii(idim),hmaxbii(idim)
      LOGICAL RedoNeighbours,BiSection
      
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/btree'

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine iterate_density')
c
c--Initialise
c
      ncalc = nlst_end - nlst_in + 1
      jneigh = 0
      its = 0
      RedoNeighbours = .false.
      nlst_beg = nlst_in
      nlst_fin = nlst_end
      irestrict = 0
      BiSection = .false.

      DO i=nlst_beg,nlst_fin
         j = list(i)
         iredolist(i) = j
         gradhs(1,j) = 0.
         gradhs(2,j) = 0.
      ENDDO
c
c--We calculate the density self-consistently with the smoothing length
c  by taking a Newton-Raphson iteration
      
      DO WHILE ((ncalc.GT.0).AND.(its.LE.maxits))
         its = its + 1
c
c--recalculate the neighbour lists if necessary
c  (I have made this a direct call to treef so that we only recalculate
c   neighbours on particles in the redo list)
c
c         WRITE(iprint,*) 'Density iteration ',its,' ncalc = ',ncalc
c         IF (ncalc.LE.5 .AND. ncalc.GT.0) THEN
c            WRITE(iprint,*) 'redolist = ',
c     &           (iredolist(i),' nneigh=',nneigh(i),';',i=1,ncalc)
c         ENDIF

         IF (RedoNeighbours) THEN
            WRITE(iprint,*) 'Recalculating neighbours...'
            DO j=nlst_beg,nlst_fin
               i = iredolist(j)
               CALL treef(i,npart,xyzmh,acc,0,fsx,fsy,fsz,epot)
            ENDDO
         ENDIF
c
c--calculate density using current h value
c
         CALL densityi(jneigh,npart,xyzmh,vxyzu,
     &               nlst_beg,nlst_fin,iredolist,itime)
c
c--calculate gradient terms, take Newton-Raphson iteration
c
         ncalc = 0
         RedoNeighbours = .false.
         DO j=nlst_beg,nlst_fin
            i = iredolist(j)
            pmassi = xyzmh(4,i)
            hi = xyzmh(5,i)
            rhoi = pmassi/(hi/hfact)**3
            dhdrhoi = -hi/(3.*rhoi)
            omegai = 1. - dhdrhoi*gradhs(1,i)
            IF (omegai.GT.tiny) gradhs(1,i) = 1./omegai
            gradhs(2,i) = dhdrhoi*gradhs(2,i)
            
            print*,i,'rho=',rho(i),rhoi,'gradh= ',gradhs(1,i)
c
c--Newton Raphson
c            
            func = rhoi - rho(i)
            IF (.not. BiSection) THEN
               dfdh = omegai/dhdrhoi
!               hnew = hi - func/dfdh
!               IF (hnew.le.0. .or. omegai.LT.tiny) THEN
!                  WRITE(iprint,*) 'using fixed point ',i,hnew,omegai
                  hnew = hfact*(pmassi/rho(i))**third
!               ENDIF
            
            ELSE
               IF (func .lt. 0.) THEN
                  hmaxbii(i) = hi
               ELSE
                  hminbii(i) = hi 
               ENDIF
               hnew = 0.5*(hminbii(i) + hmaxbii(i))
c               WRITE(iprint,*) i,'bisection, h = ',hnew,func
            ENDIF
c
c--work out whether particle is converged, if not add to list to redo
c           
            IF (abs(hnew-hi)/hi.gt.tol 
     &          .OR. omegai .LT.0. .OR. nneigh(i).GT.nneighmax) THEN
c
c--don't let number of neighbours get too big
c
               IF (hnew.GT.hi .OR. nneigh(i).GT.nneighmax) THEN
                  IF (nneigh(i).GT.nneighmax) THEN
c--stop iterations on this particle
                     irestrict = irestrict + 1
c                     WRITE(iprint,*) 'restricting h growth on part ', 
c     &            i,'h=',hi,hnew,'rho=',rho(i),'neigh = ',nneigh(i)
                     hi = 0.95*hi
                     xyzmh(5,i) = hi
                     ifakelist(1) = i
c--calculate density/gradh with this h and then exit
                     CALL densityi(jneigh,npart,xyzmh,vxyzu,
     &                    1,1,ifakelist,itime)
                     rhoi = pmassi/(hi/hfact)**3
                     dhdrhoi = -hi/(3.*rhoi)
                     gradhs(1,i) = 1.
                     gradhs(2,i) = 1.
                  ELSE
                     RedoNeighbours = .true.                  
                     ncalc = ncalc + 1
                     newlist(ncalc) = i
                     xyzmh(5,i) = hnew
                  ENDIF
               ELSE
                  ncalc = ncalc + 1
                  newlist(ncalc) = i
                  xyzmh(5,i) = hnew
               ENDIF

               IF (its.GE.(maxits-1)) THEN
                  WRITE(iprint,*) i,'h(i),hnew = ',hi,hnew,
     &             ' error = ',abs(func)/rho(i),abs(hnew-hi)/hi
               ENDIF
            ENDIF
         ENDDO

         nlst_beg = 1
         nlst_fin = ncalc
         DO i=1,ncalc
            iredolist(i) = newlist(i)
         ENDDO
c
c--copy h, grad h to ghost particles
c
         DO i=npart+1,npart+nghost
            j = ireal(i)
            rho(i) = rho(j)
            xyzmh(5,i) = xyzmh(5,j)
            gradhs(1,i) = gradhs(1,j)
            gradhs(2,i) = gradhs(2,j)
         ENDDO
c
c--switch to Bisection if not converging
c
         IF (its.EQ.maxitsnr) THEN
            WRITE(iprint,*) 'Newton-Raphson not converging:'
            WRITE(iprint,*) 'starting Bisection, ncalc=',ncalc
            DO i=1,ncalc
               j = iredolist(i)
               hminbii(j) = 0.5*xyzmh(5,j)
               !--set hmax so that h(j) is the initial guess
               hmaxbii(j) = 2.*xyzmh(5,j) - hminbii(j)
               WRITE(iprint,*) j,'h=',xyzmh(5,j),' hmin = ',hminbii(j),
     &                           ' hmax = ',hmaxbii(j)
            ENDDO
            BiSection = .true.
         ENDIF

      ENDDO
c
c--if *still* not converged something very wrong
c
      IF (its.GE.maxits) THEN
         WRITE(iprint,*) 'WARNING: DENSITY NOT CONVERGED on ',ncalc, 
     &                   ' PARTICLES'
         IF (ncalc.le.10) THEN
            WRITE(iprint,*) (iredolist(i),gradhs(1,i),
     &                       xyzmh(5,i),rho(i),i=1,ncalc)
         ENDIF
         CALL wdump(idisk1)
         CALL quit
      ENDIF
c
c--divide divv and curlv by gradh now that we know it
c
      DO j=nlst_in,nlst_end
         i = list(j)
         divv(i) = divv(i)*gradhs(1,i)
         curlv(i) = curlv(i)*gradhs(1,i)
      ENDDO
      
      WRITE(iprint,*) 'Finished density, its = ',its,
     &     ' on ',nlst_end-nlst_in+1,' particles'
      IF (irestrict.GT.0) WRITE(iprint,*) 
     & 'restricted h growth on ',irestrict,' particles'

      IF (itrace.EQ.'all') WRITE (iprint,600)
  600 FORMAT ('exit subroutine iterate_density')

      RETURN
      END
