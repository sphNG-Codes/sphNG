      SUBROUTINE divBiterate(dtmax,nlst_in,nlst_end,npart,list,
     &     xyzmh,rho,Bxyz)
c          ,moresweep,nit,error)
  
      INCLUDE 'idim'

      DIMENSION list(idim)
      DIMENSION Bxyz(3,imhd)
      DIMENSION xyzmh(5,idim)
      REAL*4 rho(idim)

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/divcurlB'

      PARAMETER (nswmax = 200000)
      PARAMETER (eta = - 1.0E-2)
      PARAMETER (weight = 1.0/1.2**3)

      REAL Bxyznew(3,imhd)
      REAL dBxyz(3,imhd)

      INTEGER nosweep
      REAL dW,dr
      LOGICAL*1 moresweep
      REAL dx,dy,dz

      PARAMETER (icompactmax=100*idim)
      DIMENSION vari(7,idim),varij(4,icompactmax)
      DIMENSION ivar(2,idim),ijvar(icompactmax)

      PARAMETER (ntests=10)
      REAL xmaxerrold(ntests)

      REAL maxerrE2
c
c--Set up constants in Code units
c
      tolerance = 1.0e-3
      dtimax = dtmax/imaxstep
c
c--Set errors to zero for iteration start
c
      numoscillations = 0
      numequal = 0
      numcomp = 0
      ipos = 1
      DO i = 1, ntests
         xmaxerrold(i) = 0.0
      END DO
c
c--Make compact list of neighbours
c
      ihasghostcount = 0
      ihasghost = 0
      icompact = 0
      Bxyzmax = 0.0
      DO n = nlst_in, nlst_end
        i = list(n)
        rneigh = radkernel*xyzmh(5,i)

        CALL getneighi(i,xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),rneigh,
     &       numneighi,neighlist,xyzmh)

        ivar(1,n) = numneighi
        ivar(2,n) = icompact
c
c--NOTE: n loop cannot be parallised because of this loop:
c
        DO k = 1, numneighi
           icompact = icompact + 1
           IF (icompact.GT.icompactmax) THEN
              WRITE (*,*) 'ERROR - compact not large enough'
              STOP
           ENDIF
           ijvar(icompact) = neighlist(k)
           IF(i.EQ.ijvar(icompact)) THEN
              PRINT *,"TRAPIMPL: Particle interacting with itself"
              PRINT *,"TRAPIMPL: Error. Bye."
              STOP
           END IF
        END DO
      END DO

      totalmagenergy = 0.

C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,list)
C$OMP& shared(icall,dtimax,dtmax,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij)
C$OMP& shared(xyzmh,dvtable,grwij,cnormk)
C$OMP& shared(nlst0,ihasghost)
C$OMP& shared(Bxyznew,Bxyz,dBxyz)
C$OMP& private(n,i,j,k,icompact)
C$OMP& private(rxyi,ryzi,rxzi,rxxi,ryyi,rzzi)
C$OMP& private(dti,dx,dy,dz,dtn)
C$OMP& private(rij2,rij,rij1,dr,pmj,rhoj,hi,hi21,hi41)
C$OMP& private(v2,v,index,dxx,index1,dgrwdx,grwtij,dW)
C$OMP& private(pmjdWrij1rhoj,runix,runiy,runiz,denom1)
C$OMP& private(dBxi,dByi,dBzi,projB)
C$OMP& reduction(+:ihasghostcount)
C$OMP& reduction(+:totalmagenergy)
C$OMP& reduction(MAX:Bxyzmax)
c
c--Copy arrays for all particles
c
C$OMP DO SCHEDULE(static)
      DO n = nlst_in, nlst_end
      	 i = list(n)

         Bxyznew(1,i) = Bxyz(1,i)
         Bxyznew(2,i) = Bxyz(2,i)
         Bxyznew(3,i) = Bxyz(3,i)

         Bxyzmax = MAX(Bxyzmax,Bxyz(1,i))
         Bxyzmax = MAX(Bxyzmax,Bxyz(2,i))
         Bxyzmax = MAX(Bxyzmax,Bxyz(3,i))

         totalmagenergy = totalmagenergy + xyzmh(4,i)*(Bxyz(1,i)**2 + 
     &        Bxyz(2,i)**2 + Bxyz(3,i)**2)/rho(i)
      END DO
C$OMP END DO
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(static)
      DO n = nlst_in, nlst_end
         i = list(n)

         IF(hasghost(i)) THEN
            ihasghostcount = ihasghostcount + 1
         ENDIF

         rxyi = 0.
         rxzi = 0.
         ryzi = 0.
         rxxi = 0.
         ryyi = 0.
         rzzi = 0.

         dBxi = 0.
         dByi = 0.
         dBzi = 0.

         hi = xyzmh(5,i)
         hi21 = 1./(hi*hi)
         hi41 = hi21*hi21

         IF (icall.EQ.1) THEN
            IF (isteps(i).EQ.0) THEN
               dti = dtmax*1.0d-12
            ELSE
               dti = dtimax*isteps(i)*1.0d-12
            ENDIF
         ELSEIF (n.LE.nlst0) THEN
            dti = dtimax*isteps(i)
         ELSE
            dti = dtimax/2.0*isteps(i)
         ENDIF
         dtn = dti*eta
         neireal = 0

         DO k = 1, ivar(1,n)
            icompact = ivar(2,n) + k
            j = ijvar(icompact)

            dx = xyzmh(1,i) - xyzmh(1,j)
            dy = xyzmh(2,i) - xyzmh(2,j)
            dz = xyzmh(3,i) - xyzmh(3,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
            dr = rij

            pmj = xyzmh(4,j)
            rhoj = rho(j)

            v2 = rij2*hi21
            v = rij/hi
            index = v2/dvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable
            dgrwdx = (grwij(index1) - grwij(index))/dvtable
            grwtij = (grwij(index) + dgrwdx*dxx)*hi41
            dW = grwtij * cnormk

            pmjdWrij1rhoj = pmj*dW*rij1/rhoj
            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1

            rxyi = rxyi + pmjdWrij1rhoj*runix*runiy
            ryzi = ryzi + pmjdWrij1rhoj*runiy*runiz
            rxzi = rxzi + pmjdWrij1rhoj*runix*runiz
            rxxi = rxxi + pmjdWrij1rhoj*(5.0*runix*runix - 1.0)
            ryyi = ryyi + pmjdWrij1rhoj*(5.0*runiy*runiy - 1.0)
            rzzi = rzzi + pmjdWrij1rhoj*(5.0*runiz*runiz - 1.0)

            varij(1,icompact) = pmjdWrij1rhoj
            varij(2,icompact) = runix
            varij(3,icompact) = runiy
            varij(4,icompact) = runiz

            projB = (Bxyz(1,i)-Bxyz(1,j))*runix + 
     &      (Bxyz(2,i)-Bxyz(2,j))*runiy + (Bxyz(3,i)-Bxyz(3,j))*runiz

            dBxi = dBxi + pmjdWrij1rhoj*(5.0*projB*runix - 
     &           (Bxyz(1,i) - Bxyz(1,j)))
            dByi = dByi + pmjdWrij1rhoj*(5.0*projB*runiy - 
     &           (Bxyz(2,i) - Bxyz(2,j)))
cc            dByi = dByi - pmj*dW/rhoj*projB
cc            dByi = dByi - weight*projB*dW/hi41/hi
            dBzi = dBzi + pmjdWrij1rhoj*(5.0*projB*runiz - 
     &           (Bxyz(3,i) - Bxyz(3,j)))

         END DO

c         dBxyz(1,i) = dtn*dBxi
c         dBxyz(2,i) = dtn*dByi
c         dBxyz(3,i) = dtn*dBzi
         dBxyz(1,i) = dBxi
         dBxyz(2,i) = dByi
         dBxyz(3,i) = dBzi

         rxyi = 5.0*dtn*rxyi
         ryzi = 5.0*dtn*ryzi
         rxzi = 5.0*dtn*rxzi
         rxxi = 1.0 - dtn*rxxi
         ryyi = 1.0 - dtn*ryyi
         rzzi = 1.0 - dtn*rzzi

         denom1 = - 1.0/(2.0*rxyi*ryzi*rxzi - rxxi*ryyi*rzzi + 
     &        rxxi*ryzi**2 + ryyi*rxzi**2 + rzzi*rxyi**2)

         vari(1,n) = (ryyi*rzzi - ryzi**2)*denom1
         vari(2,n) = (rxyi*rzzi + rxzi*ryzi)*denom1
         vari(3,n) = (rxzi*ryyi + rxyi*ryzi)*denom1
         vari(4,n) = (rxxi*rzzi - rxzi**2)*denom1
         vari(5,n) = (ryzi*rxxi + rxyi*rxzi)*denom1
         vari(6,n) = (ryyi*rxxi - rxyi**2)*denom1
         vari(7,n) = dtn

c         IF (n.LT.100) THEN
c            print *,vari(1,n),vari(2,n),vari(3,n),vari(4,n),dtn
c         ENDIF
      END DO
C$OMP END DO
C$OMP SINGLE
c      print *, 'entry single',ekcle(2,22),npart,nghost,ihasghostcount

      IF (ihasghostcount.GE.1) ihasghost = 1
C$OMP END SINGLE
C$OMP DO SCHEDULE(static)
      DO i = npart + 1, npart + nghost*ihasghost
         j = ireal(i)
         DO k = 1, 7
c            dvdx(k,i) = dvdx(k,j)
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
c
c--Begin iterating
c
c      print *, 'begin iterations',ekcle(2,22)
      DO nosweep = 1, nswmax
c      print *, 'it ',nosweep
c
c--Set error to zero for this iteration    
c
         maxerrE2 = 0.0
c
c--Calculate fluxlimiter values without using separate subroutine
c
C$OMP PARALLEL default(none)
C$OMP& shared(list,vari,ivar,varij2,ijvar,varij)
C$OMP& shared(nosweep,Bxyz,Bxyznew,Bxyzmax)
C$OMP& shared(npart,nghost,ihasghost,ireal,nlst_in,nlst_end)
C$OMP& shared(origEU,moresweep,xyzmh,gamma,alpha,beta,Rg,gmw)
C$OMP& shared(udens,radconst,nlstall,iflag,icall,ifsvi)
C$OMP& private(i,j,k,n,projB,projBx,projBy,projBz,dtn)
C$OMP& private(icompact,pmjdWrij1rhoj,runix,runiy,runiz)
C$OMP& private(moresweep2,Bxnew,Bynew,Bznew)
C$OMP& reduction(MAX:maxerrE2)
c
c--Particle I loop
c
C$OMP DO SCHEDULE(static)
         DO n = nlst_in, nlst_end
            i = list(n)

c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) i,n

            projBx = 0.
            projBy = 0.
            projBz = 0.

            dtn = vari(7,n)
c
c--All the neighbours loop
c
            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)

               pmjdWrij1rhoj = varij(1,icompact)
               runix = varij(2,icompact)
               runiy = varij(3,icompact)
               runiz = varij(4,icompact)

               projB = Bxyznew(1,j)*runix + Bxyznew(2,j)*runiy + 
     &              Bxyznew(3,j)*runiz

               projBx = pmjdWrij1rhoj*Bxyznew(1,j) - 5.0*projB*runix
               projBy = pmjdWrij1rhoj*Bxyznew(2,j) - 5.0*projB*runiy
               projBz = pmjdWrij1rhoj*Bxyznew(3,j) - 5.0*projB*runiz

            END DO              !J-loop
c
c--Calculate new B field
c
            Bxnew = (Bxyz(1,i) + dtn*projBx)*vari(1,n) + 
     &           (Bxyz(2,i) + dtn*projBy)*vari(2,n) +
     &           (Bxyz(3,i) + dtn*projBz)*vari(3,n)
            Bynew = (Bxyz(2,i) + dtn*projBy)*vari(4,n) + 
     &           (Bxyz(1,i) + dtn*projBx)*vari(2,n) +
     &           (Bxyz(3,i) + dtn*projBz)*vari(5,n)
            Bznew = (Bxyz(3,i) + dtn*projBz)*vari(6,n) + 
     &           (Bxyz(1,i) + dtn*projBx)*vari(3,n) +
     &           (Bxyz(2,i) + dtn*projBy)*vari(5,n)

c            print *,n,Bxyz(1,i),Bxyz(2,i),Bxyz(3,i),Bxnew,Bynew,Bznew
c
c--And the error is...
c
            IF (.TRUE.) THEN
               maxerrE2 = MAX(maxerrE2,
     &              ABS((Bxyznew(1,i) - Bxnew) /Bxyzmax))
               maxerrE2 = MAX(maxerrE2,
     &              ABS((Bxyznew(2,i) - Bynew) /Bxyzmax))
               maxerrE2 = MAX(maxerrE2,
     &              ABS((Bxyznew(3,i) - Bznew) /Bxyzmax))
            ELSE
c               maxerrE2 = MAX(maxerrE2,maxerrE2)
            ENDIF
c 
c--Copy values
c
            Bxyznew(1,i) = Bxnew
            Bxyznew(2,i) = Bynew
            Bxyznew(3,i) = Bznew

         END DO ! I-loop
C$OMP END DO
C$OMP DO SCHEDULE(static)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
         END DO
C$OMP END DO
C$OMP END PARALLEL
         IF (MOD(nosweep,1).EQ.0) THEN
            PRINT *,"DIVB: Finished iteration ",
     $        nosweep," of ",nswmax," (",maxerrE2,")"
         ENDIF
c
c--The actual test
c
         IF(maxerrE2.LE.tolerance) THEN
c            PRINT *,"Complete with ",nosweep," iterations"
            GOTO 150
         ENDIF
c
c--Test for oscillations
c
c         xmaxerrtot = MAX(maxerrE2,maxerrU2)
c         IF (xmaxerrtot.GT.xmaxerrold(1) .AND. nosweep.GT.10) THEN
c            numoscillations = numoscillations + 1
c            IF (numoscillations.GT.100) THEN
c               PRINT *,"GSIMPL: Oscillating ",xmaxerrtot,
c     &              (xmaxerrold(ii),ii=1,ntests)
c               moresweep = .TRUE.
c               RETURN
c            ENDIF
c         ENDIF
c
c--Test for convergence to non-zero value (incorrect minimum)
c     Must have equal value at least twice to stop detecting up and down
c     as non-convergence
c
c         GOTO 333

c         DO itest = ntests,1,-1
c
cc            GOTO 332
c
c            IF (nosweep.GT.10) THEN
c               IF (xmaxerrtot.GT.0.99999*xmaxerrold(itest).AND.
c    &              xmaxerrtot.LT.1.00001*xmaxerrold(itest)) THEN
c                  DO iii = 1, numcomp
c                     IF (xmaxerrtot.GT.0.99999*xmaxerrcomp(iii).AND.
c     &                    xmaxerrtot.LT.1.00001*xmaxerrcomp(iii)) THEN
c                        PRINT *,"GSIMPL: Non-convergence ",xmaxerrtot,
c     &                    xmaxerrcomp(iii),(xmaxerrold(ii),ii=1,ntests)
c                        moresweep = .TRUE.
c                        RETURN
c                     ENDIF
c                  END DO
c                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
c                  xmaxerrcomp(numcomp) = xmaxerrtot
c                  PRINT *,"GSIMPL: Almost Non-convergence ",
c     &                 xmaxerrtot,(xmaxerrold(ii),ii=1,ntests)
c               ENDIF
c            ENDIF
c 332        CONTINUE
c            IF (itest.NE.1) THEN
c               xmaxerrold(itest) = xmaxerrold(itest-1)
c            ELSE
c               xmaxerrold(itest) = xmaxerrtot
c            ENDIF
c         END DO
c 333     CONTINUE
      END DO ! Iterations loop
c
c--Maximum number of iterations reached
c
      PRINT *,"divBiterate: Warning. Maximum iterations reached"
      moresweep = .TRUE.
      RETURN
      STOP
c
c--Output success
c
 150  PRINT *,"Succeeded with ",nosweep," iterations ",maxerrE2
      nit = nosweep
      error = maxerrE2
c
c--And that done, return everything to ASS
c
      inum = 6660
      print *,'Initial field ',Bxyz(1,inum),Bxyz(2,inum),Bxyz(3,inum)
      print *,'Euler field ',Bxyz(1,inum)+dBxyz(1,inum),
     &     Bxyz(2,inum)+dBxyz(2,inum),
     &     Bxyz(3,inum)+dBxyz(3,inum)
      print *,'Delta field ',dBxyz(1,inum),
     &     dBxyz(2,inum),
     &     dBxyz(3,inum)
      print *,'Impl Delta ',Bxyznew(1,inum)-Bxyz(1,inum),
     &     Bxyznew(2,inum)-Bxyz(2,inum),Bxyznew(3,inum)-Bxyz(3,inum)

      print *,'Initial magnetic energy = ',totalmagenergy
      totalmagenergy = 0.

cC$OMP PARALLEL DO SCHEDULE(static) default(none)
cC$OMP& shared(nlst_in,nlst_end,list)
cC$OMP& shared(Bxyznew,Bxyz,xyzmh,rho,dBxyz)
cC$OMP& private(n,i,j,k,icompact)
cC$OMP& reduction(+:totalmagenergy)
c
c--Copy arrays for all particles
c
      DO n = nlst_in, nlst_end
         i = list(n)
         
         if (abs(xyzmh(2,i)).lt.0.1 .AND. abs(xyzmh(3,i)).lt.0.1) then
            print*,i,' x = ',xyzmh(1,i),
     &               ' divB = ',divcurlB(1,i),dBxyz(2,i)
            write (33,*) xyzmh(1,i),Bxyz(1,i),divcurlB(1,i),
     &           dBxyz(1,i),dBxyz(2,i),dBxyz(3,i)
         endif

c         Bxyz(1,i) = Bxyznew(1,i)
c         Bxyz(2,i) = Bxyznew(2,i)
c         Bxyz(3,i) = Bxyznew(3,i)

         Bxyz(1,i) = Bxyz(1,i) + dBxyz(1,i)
         Bxyz(2,i) = Bxyz(2,i) + dBxyz(2,i)
         Bxyz(3,i) = Bxyz(3,i) + dBxyz(3,i)

         totalmagenergy = totalmagenergy + xyzmh(4,i)*(Bxyz(1,i)**2 +
     &        Bxyz(2,i)**2 + Bxyz(3,i)**2)/rho(i)
      END DO
cC$OMP END PARALLEL DO

      print *,'Final magnetic energy = ',totalmagenergy
      print *,'New field ',Bxyz(1,inum),Bxyz(2,inum),Bxyz(3,inum)

      stop

      RETURN

      END
