      SUBROUTINE indexx3merge(nactive,listp,xyzmh,next1p,next2p,next3p)
c************************************************************
c                                                           *
c  Subroutine by M. Bate (15/09/06).                        *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION listp(nactive), next1p(nactive), next2p(nactive), 
     &     next3p(nactive)
      DIMENSION xyzmh(5,mmax2)

      PARAMETER (nmaxthreads = 64)
      DIMENSION indexs(nmaxthreads),indexe(nmaxthreads),
     &     ioffset(nmaxthreads)
      DIMENSION ipoint(nmaxthreads)

      INTEGER inumthreads
      INTEGER OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/treecom_P'

      nhalf = nactive/2
      nrest = nactive - nhalf

c      CALL getused(time1)

c
c--Note: Needs to copy from one list to another when merging together sub-lists
c     To save memory, it sorts originally into
c     key for next1
c     ihash for next2
c     nay for next3
c
C$OMP PARALLEL default(none)
C$OMP& shared(nactive,nhalf,nrest,imax)
C$OMP& shared(listp,xyzmh,next1p,next2p,next3p,indexs,indexe,ioffset)
C$OMP& shared(key,ihash,nay)
C$OMP& private(inumthreads,i,istart,inumber,xnext,jposold,ilistnumber)
C$OMP& private(j,k,indexj,inext,xnow,imove,jpos,iswap,ipoint)
C$OMP& private(listpnow,listpnext)
c
c--Sort the list in NUM_THREADS sections
c
      inumthreads = OMP_GET_NUM_THREADS()
c      print *,'num threads', inumthreads

c      inumthreads = 4
      IF(inumthreads.GT.nmaxthreads) THEN
         WRITE (*,*) 'ERROR - inumthreads.GT.nmaxthreads'
         CALL quit
      ENDIF
c
c--Produce an x-sorted list
c
cC$OMP SINGLE
c      print *,'Threads ',inumthreads
cC$OMP END SINGLE      

C$OMP DO SCHEDULE(static)
      DO i = 1, inumthreads
         istart = 1 + (i-1)*(nactive/inumthreads)
         inumber = nactive/inumthreads
         indexs(i) = istart
         ioffset(i) = istart - 1
         IF (i.EQ.inumthreads) inumber = nactive - istart + 1
         indexe(i) = istart + inumber - 1
         CALL indexx3unique(inumber, listp(istart), xyzmh,1,key(istart))
      END DO
C$OMP END DO
c
c--Produce a y-sorted list
c
cC$OMP SINGLE
c      print *,'Sorted 1 '
c      do i = 1, nactive
c         write (51,*) listp(key(i)),xyzmh(1,listp(key(i)))
c      end do
c      do i = 1, inumthreads
c         do j = indexs(i),indexe(i)
c            write (51,*) listp(key(j)+ioffset(i)),
c     &           xyzmh(1,listp(key(j)+ioffset(i)))
c         end do
c      end do
c      do i = 1, 10
c         print *,listp(key(i)),listp(key(i+(nactive/inumthreads)))
c      end do
cC$OMP END SINGLE      

C$OMP DO SCHEDULE(static)
      DO i = 1, inumthreads
         istart = 1 + (i-1)*(nactive/inumthreads)
         inumber = nactive/inumthreads
         IF (i.EQ.inumthreads) inumber = nactive - istart + 1
         CALL indexx3unique(inumber,listp(istart),xyzmh,2,ihash(istart))
      END DO
C$OMP END DO
c
c--Produce a z-sorted list
c
cC$OMP SINGLE
c      print *,'Sorted 2 '
cC$OMP END SINGLE      

C$OMP DO SCHEDULE(static)
      DO i = 1, inumthreads
         istart = 1 + (i-1)*(nactive/inumthreads)
         inumber = nactive/inumthreads
         IF (i.EQ.inumthreads) inumber = nactive - istart + 1
         CALL indexx3unique(inumber, listp(istart), xyzmh,3,nay(istart))
      END DO
C$OMP END DO
c
c--Merge sorted sections into completely sorted lists
c     This is parallel over 6 processors
c
cC$OMP SINGLE
c      print *,'Sorted 3 '
cC$OMP END SINGLE      

      IF (inumthreads.GT.1) THEN
C$OMP DO SCHEDULE(static)
      DO i = 0, 5
c         print *,'Doing ',i
         jposold = 0
         ilistnumber = MOD(i,3)
c         print *,'Doing list ',ilistnumber
         IF (ilistnumber.EQ.0) THEN
            IF (MOD(i,2).EQ.0) THEN
c
c--Sort from start of x-list
c
               DO j = 1, inumthreads
                  ipoint(j) = indexs(j)
               END DO
               DO k = 1, nhalf
c                  print *,'k is ',k,nhalf
                  xnext = +1.0E+30
                  listpnext = imax
                  imove = imax
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = key(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(1,listpnow)
c                        print *,'inext, xnow ',inext,xnow,xnext,imove
                        IF (xnow.LT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.LT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
c                           print *,'update ',xnext,imove,jpos
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
c                  if (k.LT.100) 
c     &                 print *,'swapping ',imove,k,xnext,jpos
                  IF (imove.EQ.imax) THEN
                     print *,'imove.EQ.imax'
                     CALL quit
                  ENDIF
                  next1p(k) = imove
                  ipoint(jpos) = ipoint(jpos) + 1
                  IF (ipoint(jpos).GT.indexe(jpos)) ipoint(jpos) = 0
               END DO
c
c--Sort from end of x-list
c
            ELSE
               DO j = 1, inumthreads
                  ipoint(j) = indexe(j)
               END DO
               DO k = nactive, nactive-nrest+1, -1
                  xnext = -1.0E+30
                  listpnext = 0
                  imove = 0
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = key(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(1,listpnow)
                        IF (xnow.GT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.GT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
                  IF (imove.EQ.0) THEN
                     print *,'imove.EQ.0'
                     CALL quit
                  ENDIF
                  next1p(k) = imove
                  ipoint(jpos) = ipoint(jpos) - 1
                  IF (ipoint(jpos).LT.indexs(jpos)) ipoint(jpos) = 0
               END DO
            ENDIF
         ELSEIF(ilistnumber.EQ.1) THEN
            IF (MOD(i,2).EQ.0) THEN
c
c--Sort from start of y-list
c
               DO j = 1, inumthreads
                  ipoint(j) = indexs(j)
               END DO
               DO k = 1, nhalf
                  xnext = 1.0E+30
                  listpnext = imax
                  imove = imax
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = ihash(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(2,listpnow)
                        IF (xnow.LT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.LT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
                  IF (imove.EQ.imax) THEN
                     print *,'imove.EQ.imax'
                     CALL quit
                  ENDIF
                  next2p(k) = imove
                  ipoint(jpos) = ipoint(jpos) + 1
                  IF (ipoint(jpos).GT.indexe(jpos)) ipoint(jpos) = 0
               END DO
c
c--Sort from end of y-list
c
            ELSE
               DO j = 1, inumthreads
                  ipoint(j) = indexe(j)
               END DO
               DO k = nactive, nactive-nrest+1, -1
                  xnext = -1.0E+30
                  listpnext = 0
                  imove = 0
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = ihash(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(2,listpnow)
                        IF (xnow.GT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.GT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
                  IF (imove.EQ.0) THEN
                     print *,'imove.EQ.0'
                     CALL quit
                  ENDIF
                  next2p(k) = imove
                  ipoint(jpos) = ipoint(jpos) - 1
                  IF (ipoint(jpos).LT.indexs(jpos)) ipoint(jpos) = 0
               END DO
            ENDIF
         ELSE
            IF (MOD(i,2).EQ.0) THEN
c
c--Sort from start of z-list
c
               DO j = 1, inumthreads
                  ipoint(j) = indexs(j)
               END DO
               DO k = 1, nhalf
                  xnext = 1.0E+30
                  listpnext = imax
                  imove = imax
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = nay(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(3,listpnow)
                        IF (xnow.LT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.LT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
                  IF (imove.EQ.imax) THEN
                     print *,'imove.EQ.imax'
                     CALL quit
                  ENDIF
                  next3p(k) = imove
                  ipoint(jpos) = ipoint(jpos) + 1
                  IF (ipoint(jpos).GT.indexe(jpos)) ipoint(jpos) = 0
               END DO
c
c--Sort from end of z-list
c
            ELSE
               DO j = 1, inumthreads
                  ipoint(j) = indexe(j)
               END DO
               DO k = nactive, nactive-nrest+1, -1
                  xnext = -1.0E+30
                  listpnext = 0
                  imove = 0
                  DO j = 1, inumthreads
c
c--Find particle that should be next in the list
c
                     indexj = ipoint(j)
                     IF (indexj.NE.0) THEN
                        inext = nay(indexj) + ioffset(j)
                        listpnow = listp(inext)
                        xnow = xyzmh(3,listpnow)
                        IF (xnow.GT.xnext .OR.
     &                   xnow.EQ.xnext .AND. listpnow.GT.listpnext) THEN
                           listpnext = listpnow
                           xnext = xnow
                           imove = inext
                           jpos = j
                        ENDIF
                     ENDIF
                  END DO
c
c--Put this particle in the correct place (and exchange out the one that is
c     already there)
c               
                  IF (imove.EQ.0) THEN
                     print *,'imove.EQ.0'
                     CALL quit
                  ENDIF
                  next3p(k) = imove
                  ipoint(jpos) = ipoint(jpos) - 1
                  IF (ipoint(jpos).LT.indexs(jpos)) ipoint(jpos) = 0
               END DO
            ENDIF
         ENDIF
      END DO
C$OMP END DO
      ENDIF
C$OMP END PARALLEL

c      CALL getused(time2)

c      print *,'Time ',time2-time1

      RETURN
      END
