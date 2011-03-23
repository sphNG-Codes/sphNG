      SUBROUTINE ASS(nlst_in,nlst_end,list,dt,itime,
     $     npart,ntot,xyzmh,vxyzu,ekcle,trho,dedxyz,alphaMM)

      INCLUDE 'idim'
      
      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      REAL*4 trho(idim),alphaMM(idim)
      DIMENSION list(idim)
      DIMENSION dedxyz(3,iradtrans)

      RETURN

      END
