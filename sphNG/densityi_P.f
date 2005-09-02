      SUBROUTINE densityi (npart,xyzmh,vxyzu,
     &            nlst_in,nlst_end,list,itime)
c************************************************************
c                                                           *
c  Subroutine to compute the density and the velocity       *
c     divergence at ipart particle location. This           *
c     subroutine uses the binary tree algorithm to locate   *
c     neighbors.                                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim), list(idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/sort'

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityi')
c
c--Initialise
c
      rhonext = 0.
      icreate = 0
      radcrit2 = radcrit*radcrit

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,list,divv,curlv)
C$OMP& shared(hmax,xyzmh,vxyzu,pr,vsound,rho)
C$OMP& shared(nneigh,neighb,neighover,selfnormkernel)
C$OMP& shared(cnormk,radkernel,dvtable,wij,grwij)
C$OMP& shared(listpm,iphase)
C$OMP& shared(iprint,nptmass,iptmass,radcrit2,iorig)
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi)
C$OMP& private(pmassj,hmean,hmean21,hmean31,hmean41)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur)
C$OMP& reduction(MAX:rhonext)
      DO 50 n = nlst_in, nlst_end
         ipart = list(n)

         IF (iphase(ipart).GE.1) GOTO 50

         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         pmassi = xyzmh(4,ipart)
         hi = xyzmh(5,ipart)
         vxi = vxyzu(1,ipart)
         vyi = vxyzu(2,ipart)
         vzi = vxyzu(3,ipart)

         rhoi = 0.
         divvi = 0.
         curlvxi = 0.
         curlvyi = 0.
         curlvzi = 0.

         DO 40 k = 1, nneigh(ipart)
            IF (k.GE.nlmax) THEN
               IF (-neighb(nlmax,ipart).GT.noverflow .OR. 
     &              -neighb(nlmax,ipart).LE.0) THEN
                  WRITE (*,*) 'ERROR - neighb(nlmax,ipart) ',
     &                 -neighb(nlmax,ipart)
                  CALL quit
               ENDIF
               IF (k-nlmax+1.GT.nlovermax) THEN
                  WRITE (*,*) 'ERROR - k-nlmax+1 ',k-nlmax+1
                  CALL quit
               ENDIF

               j = neighover(k-nlmax+1,ABS(neighb(nlmax,ipart)))

               IF (j.LE.0 .OR. j.GT.idim) THEN
                  WRITE (*,*) 'ERROR - j densityi_P ',j
                  CALL quit
               ENDIF
            ELSE
               j = neighb(k,ipart)
            ENDIF

            IF (iphase(ipart).EQ.-1 .OR. iphase(j).EQ.-1) GOTO 40
            IF (iphase(j).GE.1) THEN
               WRITE(iprint,*)'ERROR - ptmass in neighbour list ', j,
     &              iorig(j)
               CALL quit
            ENDIF
c
c--Use mean h
c
            pmassj = xyzmh(4,j)
            hj = xyzmh(5,j)
            hmean = 0.5*(hi + hj)
            hmean21 = 1./(hmean*hmean)
            hmean31 = hmean21/hmean
            hmean41 = hmean21*hmean21

            dx = xi - xyzmh(1,j)
            dy = yi - xyzmh(2,j)
            dz = zi - xyzmh(3,j)

            dvx = vxi - vxyzu(1,j)
            dvy = vyi - vxyzu(2,j)
            dvz = vzi - vxyzu(3,j)

            rij2 = dx*dx + dy*dy + dz*dz + tiny
            v2 = rij2*hmean21

            IF (v2.LT.radkernel**2) THEN
               rij1 = SQRT(rij2)
c
c--Get kernel quantities from interpolation in table
c
               index = v2/dvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dwdx = (wij(index1) - wij(index))/dvtable
               wtij = (wij(index) + dwdx*dxx)*hmean31
               dgrwdx = (grwij(index1) - grwij(index))/dvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hmean41/rij1
c
c--Compute density
c
               rhoi = rhoi + pmassj*wtij
c
c--Velocity divergence times density
c
               projv = grwtij*(dvx*dx + dvy*dy + dvz*dz)
               divvi = divvi - pmassj*projv
c
c--Velocity curl in 3D times density
c
               procurlvz = grwtij*(dvy*dx - dvx*dy)
               procurlvy = grwtij*(dvx*dz - dvz*dx)
               procurlvx = grwtij*(dvz*dy - dvy*dz)

               curlvxi = curlvxi - pmassj*procurlvx
               curlvyi = curlvyi - pmassj*procurlvy
               curlvzi = curlvzi - pmassj*procurlvz
            ENDIF
 40      CONTINUE

         rho(ipart) = cnormk*(rhoi + selfnormkernel*pmassi/(hi*hi*hi))
         divv(ipart) = cnormk*divvi
         curlv(ipart) = cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)
c
c--Pressure and sound velocity from ideal gas law...
c
         CALL eospg(ipart, vxyzu, rho, pr, vsound)
c
c--Find particle with highest density outside radcrit of point mass
c
         IF (iptmass.NE.0) THEN
            DO l = 1, nptmass
               iptcur = listpm(l)
               IF ( (xi - xyzmh(1,iptcur))**2 + 
     &              (yi - xyzmh(2,iptcur))**2 + 
     &              (zi - xyzmh(3,iptcur))**2 .LT.radcrit2) GOTO 50
            END DO
            rhonext = MAX(rhonext, rho(ipart))
         ENDIF

 50   CONTINUE
C$OMP END PARALLEL DO
c
c--Possible to create a point mass
c
      IF (rhonext.GT.rhocrea .AND. nptmass.LT.iptdim .AND.
     &                     iptmass.NE.0 .AND. icall.EQ.3) THEN
c
c--Find particle with highest density outside radcrit of point mass
c
         DO n = nlst_in, nlst_end
            ipart = list(n)
            IF (rho(ipart).EQ.rhonext) THEN
               irhonex = ipart
            ENDIF
         END DO
c
c--Make sure that all neighbours of point mass candidate are being
c     done on this time step. Otherwise, not possible to accrete
c     them to form a point mass and it may create a point mass without
c     accreting many particles!
c
         IF ((2.0*xyzmh(5,irhonex)).LT.hacc) THEN
            WRITE(iprint,*)'Ptmass creation passed h ',xyzmh(5,irhonex)

            CALL getneigh(irhonex,npart,xyzmh(5,irhonex),xyzmh,nlist,
     &           iptneigh,nearl)

            iokay = 1
            DO n = 1, nlist
               j = nearl(n)
               IF (it0(j).NE.itime) iokay = 0
            END DO
c
c--Set creation flag to true. Other tests done in accrete.f
c
            IF (iokay.EQ.1) THEN
               icreate = 1 
               WRITE(iprint,*) ' and all particles on step'
            ELSE
               WRITE(iprint,*) ' but not all particles on step'
            ENDIF
         ELSE
            WRITE(iprint,*) 'Ptmass creation failed on h ',
     &           xyzmh(5,irhonex)
         ENDIF
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityi')

      RETURN
      END
