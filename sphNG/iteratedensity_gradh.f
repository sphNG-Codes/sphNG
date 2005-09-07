      SUBROUTINE iterate_density(npart,xyzmh,vxyzu,dha,
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
      PARAMETER(maxits=100)
      PARAMETER(nneighmax=300)
      REAL tol,third,hfact
      PARAMETER(hfact=1.2) !!!,neides=INT(32./3.*pi*hfact**3))
      PARAMETER(tol=1.e-3,third=1./3.)

      REAL*4 dha(2,idim)
      DIMENSION xyzmh(5,idim),vxyzu(4,idim),list(idim)
      DIMENSION newlist(idim), iredolist(idim), ifakelist(idim)

      INTEGER ncalc,nredone
      REAL fsx,fsy,fsz,epot
      REAL pmassi,rhoi,omegai
      DIMENSION h_old(idim)
      LOGICAL NeighboursChanged
      INTEGER iRedoNeighbours(idim)
      
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/useles'

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine iterate_density')
c
c--Initialise
c
      ncalc = nlst_end - nlst_in + 1
      ncalctot = 0
      its = 0
      NeighboursChanged = .false.
      nlst_beg = nlst_in
      nlst_fin = nlst_end
      irestrict = 0
      nredone = 0

      DO i=nlst_beg,nlst_fin
         j = list(i)
         iredolist(i) = j
         gradhs(1,j) = 0.
         gradhs(2,j) = 0.
         h_old(j) = xyzmh(5,j)
         iRedoNeighbours(j) = 0
      ENDDO
c
c--We calculate the density self-consistently with the smoothing length
c  by taking a Newton-Raphson iteration
      
      DO WHILE ((ncalc.GT.0).AND.(its.LE.maxits))
         its = its + 1
         ncalctot = ncalctot + ncalc
c
c--recalculate the neighbour lists if necessary
c  (I have made this a direct call to treef so that we only recalculate
c   neighbours on particles in the redo list)
c
c         WRITE(iprint,*) 'Density iteration ',its,' ncalc = ',ncalc

         DO j=nlst_beg,nlst_fin
            i = iredolist(j)
            IF (iRedoNeighbours(i).EQ.1) THEN
               CALL treef(i,npart,xyzmh,acc,0,fsx,fsy,fsz,epot)
               NeighboursChanged = .true.
               nredone = nredone + 1
            ENDIF
         ENDDO
c
c--calculate density using current h value
c
         CALL densityi(npart,xyzmh,vxyzu,
     &               nlst_beg,nlst_fin,iredolist,itime)
c
c--calculate gradient terms, take Newton-Raphson iteration
c
         ncalc = 0
         DO j=nlst_beg,nlst_fin
            i = iredolist(j)
            iRedoNeighbours(i) = 0
            pmassi = xyzmh(4,i)
            hi = xyzmh(5,i)
            rhoi = pmassi/(hi/hfact)**3
            dhdrhoi = -hi/(3.*rhoi)
            omegai = 1. - dhdrhoi*gradhs(1,i)
            gradhs(1,i) = 1./omegai
            gradhs(2,i) = dhdrhoi*gradhs(2,i)
            
            !!print*,i,'rho=',rho(i),rhoi,'gradh= ',gradhs(1,i)
c
c--Newton Raphson
c            
            func = rhoi - rho(i)
            dfdh = omegai/dhdrhoi
            hnew = hi - func/dfdh
            IF (hnew.le.0. .or. omegai.LT.tiny) THEN
               WRITE(iprint,*) 'using fixed point ',i,hnew,omegai
               hnew = hfact*(pmassi/rho(i))**third
            ENDIF
c
c--work out whether particle is converged, if not add to list to redo
c           
            IF (abs(hnew-hi)/h_old(i).gt.tol 
     &          .OR. omegai .LT.0. .OR. nneigh(i).GT.nneighmax) THEN
c
c--don't let number of neighbours get too big
c
               IF (hnew.GT.h_old(i) .OR. nneigh(i).GT.nneighmax) THEN
                  IF (nneigh(i).GT.nneighmax) THEN
c--stop iterations on this particle
                     irestrict = irestrict + 1
c                     WRITE(iprint,*) 'restricting h growth on part ', 
c     &            i,'h=',hi,hnew,'rho=',rho(i),'neigh = ',nneigh(i)
                     hi = 0.95*hi
                     xyzmh(5,i) = hi
                     ifakelist(1) = i
c--calculate density/gradh with this h and then exit
                     CALL densityi(npart,xyzmh,vxyzu,
     &                    1,1,ifakelist,itime)
                     rhoi = pmassi/(hi/hfact)**3
                     dhdrhoi = -hi/(3.*rhoi)
                     gradhs(1,i) = 1.
                     gradhs(2,i) = 0.
                     dha(1,i) = 0. !!dhdrhoi*(-divv(i))
                  ELSE
                     iRedoNeighbours(i) = 1                  
                     ncalc = ncalc + 1
                     newlist(ncalc) = i
                     xyzmh(5,i) = hnew
                  ENDIF
               ELSE
                  ncalc = ncalc + 1
                  newlist(ncalc) = i
                  xyzmh(5,i) = hnew
                  dha(1,i) = dhdrhoi*(-divv(i))
c
c--neighbour statistics for output
c
                  IF (icall.EQ.3) THEN
                     numneigh = nneigh(i)
                     inmin = MIN(inmin,numneigh)
                     inmax = MAX(inmax,numneigh)
                     inminsy = MIN(inminsy,numneigh)
                     inmaxsy = MAX(inmaxsy,numneigh)
                     IF (hnew.LT.hmin .AND. numneigh.GT.neimin)
     &                  ioutmin = ioutmin + 1
                     IF (numneigh.GT.neimax) ioutsup = ioutsup + 1
                     IF (numneigh.LT.neimin) ioutinf = ioutinf + 1
                  ENDIF
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

      ENDDO
c
c--if not converged something very wrong
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
      
      IF (NeighboursChanged) THEN
         WRITE(iprint,*) 'Recalculating all neighbours...',nredone
         DO i=1,npart
            CALL treef(i,npart,xyzmh,acc,0,fsx,fsy,fsz,epot)
         ENDDO
      ENDIF
      WRITE(iprint,*) 'density its ',its,' total= ',
     &     ncalctot,' on ',nlst_end-nlst_in+1,' particles'
      IF (irestrict.GT.0) WRITE(iprint,*) 
     & 'restricted h growth on ',irestrict,' particles'

      IF (itrace.EQ.'all') WRITE (iprint,600)
  600 FORMAT ('exit subroutine iterate_density')

      RETURN
      END
