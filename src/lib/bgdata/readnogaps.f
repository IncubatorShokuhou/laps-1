      subroutine read_nogaps(lun,nx,ny,nz
     .               ,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     .               ,ht,tp,sh,uw,vw
     .               ,ht_sfc,pr_sfc,sh_sfc,tp_sfc
     .               ,uw_sfc,vw_sfc,mslp
     .               ,istatus)
c
      implicit none
c
c     integer nx,ny,nz,i,j,k,l,it
      integer nx,ny,nz,i,j,k,l
      integer lun,istatus
      integer nvarsmax,nvars
     .       ,nshl
     .       ,nlevs(nvarsmax)
     .       ,ivarcoord(nvarsmax)
     .       ,ivarid(nvarsmax)
c
      real   ht(nx,ny,nz),     !nogaps height (m)
     .       tp(nx,ny,nz),     !nogaps temperature (k)
     .       sh(nx,ny,nz),     !nogaps specific humidity (kg/kg) 
     .       uw(nx,ny,nz),     !nogaps u-wind (m/s)
     .       vw(nx,ny,nz),     !nogaps v-wind (m/s)
     .       dummy(nx,ny,nz)

      real   ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,sh_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

c
c     character*9   fname
c     character*4   af
c
c     real   xe,esat
      real   esat
      common /estab/esat(15000:45000)
c_______________________________________________________________________________

      istatus = 1
c
c *** open nogaps file.
c nogaps file already open, only need to read it (j. smart 9-4-98)
      print*
      print*,'read nogaps 3-d variables'
c nvar =1
      do k=1,nz
         read(lun,err=50) ((tp(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'read u'
c nvar=2
      do k=1,nz
         read(lun,err=50) ((uw(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'read v'
c nvar=3
      do k=1,nz
         read(lun,err=50) ((vw(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'read td'
c nvar=4
      nshl=nlevs(4)
      do k=1,nshl
         read(lun,err=50) ((sh(i,j,k),i=1,nx),j=1,ny) !read in as dew point depression.
      enddo
      do k=nshl+1,nz
      do j=1,ny
      do i=1,nx
         sh(i,j,k)=-99999.0
      enddo
      enddo
      enddo

c     print*,'read ht'
      do k=1,nz
         read(lun,err=50) ((ht(i,j,k),i=1,nx),j=1,ny)
      enddo
c
c read nogaps sfc variables
c
      print*,'read sfc variables'
      read(lun,err=50) ((tp_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((uw_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((vw_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((sh_sfc(i,j),i=1,nx),j=1,ny) !dew point depression 
      read(lun,err=50) ((mslp(i,j),i=1,nx),j=1,ny)
c
c nvar = 12
c currently nogaps does not provide sfc pressure. also, nogaps only
c has 10 variables in the file presently (9-4-98).
c
      if(.false.)then
      do l=12,nvars
        if(ivarid(l).eq.1.and.ivarcoord(l).eq.1)then
           read(lun,err=50) ((pr_sfc(i,j),i=1,nx),j=1,ny)
           goto  188
        else
           do k=1,nlevs(l)
              read(lun,err=50)((dummy(i,j,k),i=1,nx),j=1,ny)
           enddo
        endif
      enddo
      print*,'did not find sfc p data!'

188   continue
      endif
c
      close(lun)
c
c *** convert dew point to specific humidity.
c this code is now in readdgprep.f (j. smart 9-2-98)
c     print*,'convert dew point depression to td'
c
c return dewpoint temperature in "sh" arrays
      do k=1,nshl
      do j=1,ny
      do i=1,nx
         sh(i,j,k)=tp(i,j,k)-sh(i,j,k)
      enddo
      enddo
      enddo

      do j=1,ny
      do i=1,nx
         sh_sfc(i,j)=tp_sfc(i,j)-sh_sfc(i,j)
      enddo
      enddo
c
      istatus=0
      print*

      return
c
50    print*,'error reading nogaps files'
      return

990   continue
      print *,'error finding nogaps file.'
c
      end
