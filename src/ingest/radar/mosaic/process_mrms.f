
      subroutine process_mrms(ni,nj,nk)

      real ref_3d(ni,nj,nk)
      real ref3d_column(nk+2,ni*nj)

!     ref3d_column=-999.0
!     numref=0
!     do j=2,nlat-1
!     do i=2,nlon-1
!     numlvl=0
!     do k=1,maxlvl
!       if(abs(ref0(i,j,k)) < 888.0 ) numlvl=numlvl+1
!     enddo ! k
!     if(numlvl > 0 ) then
!       numref=numref+1
!       ref3d_column(1,numref)=float(i)
!       ref3d_column(2,numref)=float(j)
!       do k=1,maxlvl
!          ref3d_column(2+k,numref)=ref0(i,j,k)
!       enddo
!     endif
!     enddo ! i
!     enddo ! j

!     sample file location
!     /scratch2/portfolios/bmc/rtrr/rapdev3/cycle/2015032715/obsprd/refingsi.dat

      maxlvl = nk
      open(10,file='./'//'refingsi.dat',form='unformatted')
        read(10) maxlvl,nj,ni,numref,i1,i2
        write(*,*) 'dump out results', numref, 'out of', nj,ni
        read(10) ((ref3d_column(k,i),k=1,maxlvl+2),i=1,numref)
      close(10)

      return
      end
