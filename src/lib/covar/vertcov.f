c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine vert_cov(i4time,v_inno
     +                          ,num_of_ens
     +                          ,ix,iy,iz,isolevel
     +                          ,covert)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c c history:
c date         name          action
c --------     -----------   -------------------------------------------
c 03/27/2007   ok-yeon kim   created. 
c ----------------------------------------------------------------------

      implicit none

      integer       ix,iy,iz
      integer       num_of_ens,iens
      integer       i,j,k,l,m,n,ii,jj
      integer       ilength,jlength
      integer     i4time
      integer     isolevel(iz)
      real        v_inno(ix,iy,iz,20)
      real        amean(ix,iy,iz)
      real        dev(ix,iy,iz,20)
      real        var(ix,iy,iz)
      real        covert(iz,iz,ix,iy)

c--------------- end of diclaration -----------------------------------


c average of innovation at each point  !!
      amean=0.
      do n=1,num_of_ens
        do k=1,iz
          do j=1,iy
            do i=1,ix
              amean(i,j,k)=amean(i,j,k)+v_inno(i,j,k,n)
            enddo 
          enddo
        enddo
      enddo

      do k=1,iz
        do j=1,iy
          do i=1,ix
            amean(i,j,k)=amean(i,j,k)/float(num_of_ens)
          enddo
        enddo
      enddo


c deviation, variance at each point  !!
      var=0.
      do k=1,iz
        do j=1,iy  
          do i=1,ix
            do n=1,num_of_ens
               dev(i,j,k,n)=v_inno(i,j,k,n)-amean(i,j,k)  ! deviation
               var(i,j,k)=var(i,j,k)+(dev(i,j,k,n)**2.)   ! variance
            enddo
          enddo
        enddo
      enddo


c covariance in the vertical direction  !!

      iens=num_of_ens-1            
      do k=1,iz
        do l=1,iz
          do j=1,iy
            do i=1,ix

            do n=1,num_of_ens
              covert(l,k,i,j)=covert(l,k,i,j)
     +                         +(dev(i,j,k,n)*dev(i,j,l,n))
            enddo     ! n
              covert(l,k,i,j)=covert(l,k,i,j)/float(iens)

            enddo   ! i
          enddo     ! j
        enddo       ! l
      enddo         ! k


      return
      end
