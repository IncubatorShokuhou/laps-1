module input_bg_obs

   use prmtrs_stmas
   use generaltools
   use readobserves, only: rdradrobs, rdbufrobs, rdbufrobs_xie, allocatob, &
                           dealoctob, nobsmax, obp, obs, obe, nst, oba, rdobstest, rdlapsrdr
   use read_backgrd, only: rdlapsbkg, rdbckgrnd, alloctbkg, dealctbkg, bk0, c00, d00, x00, y00, p00, &
                           dx0, dy0, z00, dz0, dt0, heightu, heightl, rdbkgtest

   implicit none
   integer, allocatable :: gridmask(:, :, :, :, :)
   real, allocatable :: obradius(:, :)

   public bkgrndobs
   private readnmlst, bkgmemalc, obsmemalc, memoryalc, memoryrls, addbkgrnd, getbkgrnd_new

!**************************************************
!comment:
!   this module is used by main.f90 to get informations of background and observations.
!   subroutines:
!      bkgrndobs: main routine of this module, get in informations of background and observations.
!      readnmlst: read name list file.
!      bkgmemalc: allocate memory for background arrays.
!      obsmemalc: allocate observation memory and save the abservation informations.
!      memoryalc: get informations from the model or some initial files, allocate memory for the relevant arrays.
!      memoryrls: memory release.
!      getbkgrnd: interpolate the background field to analysis grid points.
!      addbkgrnd: set background as observations at the grid points where observations can not affect.
!      getbkgrnd_new: modified from getbkgrnd.
!
!   arrays:
!      gridmask: mask to define whethere the grid point is affected by observations.
!      obradius: this array defines the distance in which the field can be affect by the observations.
!**************************************************

contains

!doc==================================================================
!
!>
!! this is a function of input_bg_obs module and included in
!! input_bg_obs.f90 for reading in stmas 3d analysis namelist.
!! all of the namelist variables are defined in prmtrs_stmas.f90.
!!
!! \author yuanfu xie
!!
!! \b history \n
!! created: dec. 2013
!!
!endoc================================================================

!doc==================================================================
!
!>
!! this is a function of input_bg_obs module and included in
!! input_bg_obs.f90 for reading in stmas 3d analysis namelist.
!! all of the namelist variables are defined in prmtrs_stmas.f90.
!!
!! \author yuanfu xie
!!
!! \b history \n
!! created: dec. 2013
!!
!endoc================================================================

   subroutine read_namelist

      implicit none

      namelist /number_states/ numstat

      namelist /stmas3d/ ifbkgnd, ifbound, ifpcdnt, &
         penal0x, penal0y, penal0z, penal0t, &
         pnlt0pu, pnlt0pv, fnstgrd, &
         numdims, numgrid, maxgrid, &
         u_cmpnnt, v_cmpnnt, w_cmpnnt, pressure, temprtur, humidity, &
         raincont, snowcont, grapcont, cloudice, cloudwat, &
         cosstep, midgrid, finstep, pnlt0hy, taul_hy, &
         endhylv, endgslv

      character(len=256) ::filename
      integer :: n, nm, ns, ierr, limgrid(2)
      real :: ratio

      ! get namelist file from laps static:
      call get_directory('static', filename, n)
      filename(n:n + 11) = '/stmas3d.nl'

      ! open file for read:
      open (11, file=filename(1:n + 11))
      read (11, nml=number_states, iostat=ierr)

      ! allocate memory:
      allocate (sl0(numstat), penal0x(numstat), penal0y(numstat), penal0z(numstat), &
                penal0t(numstat), penal_x(numstat), penal_y(numstat), penal_z(numstat), &
                penal_t(numstat), obradius(maxdims, numstat), stat=ierr)

      read (11, nml=stmas3d, iostat=ierr)
      close (11)

      ! default scaling:
      sl0 = 1.0

      ! radar data influence radius: option later for readin from namelist
      do ns = 1, numstat
         obradius(1:2, ns) = 200000.0
         obradius(3, ns) = 50000.0
         obradius(4, ns) = 0.0
      end do

      ! initial grid positions:
      oripstn = 0

      ! coordinate indices:
      xsl = 1
      ysl = 2
      psl = 3
      csl = 4
      dsl = 5

      ! unit conversion to meters:
      xytrans = 1.0
      z_trans = 1.0

      ! for testing: if_test = 1; otherwise,
      if_test = 0

      ! multigrid v cycle option: 1 full cycle (coarse to fine repeated);
      !                           0 half cycle (one coarse to fine)
      ifrepet = 0
      itrepet = 0

      ! threshold values qc observation in the vertical and time:
      limit_3 = 10
      limit_4 = 1

      ! initializing stmas grids:
      ngptobs = 2**numdims        ! number of observation grid indices
      nallobs = 0                ! total number of all observations

      ! multigrid setup:
      if (maxgrid(1) .eq. 0) then

         ! using a default multigrid setup based on the fcstgrd:
         fnstgrd = 4       ! default levels of multigrid

         ! get laps fcstgrd:
         call get_grid_dim_xy(fcstgrd(1), fcstgrd(2), ierr)
         call get_laps_dimensions(fcstgrd(3), ierr)
         fcstgrd(4) = 3    ! hardcode for now

         ! for x and y directions:
         ! limit of 401 maxgrid in x and y directions:
         ratio = float(fcstgrd(1) - 1)/float(fcstgrd(2) - 1)
         limgrid(1) = min(fcstgrd(1), 401)
         limgrid(2) = min(fcstgrd(2), 401)
         if (ratio .ge. 1.0) then
            limgrid(2) = int((limgrid(1) - 1)/ratio) + 1
         else
            limgrid(1) = int((limgrid(2) - 1)*ratio) + 1
         end if
         do n = 1, 2
            numgrid(n) = int((limgrid(n) - 1)/2**(fnstgrd - 1))
            maxgrid(n) = 2**(fnstgrd - 1)*numgrid(n) + 1
            numgrid(n) = numgrid(n) + 1
         end do

         ! for z direction:
         nm = 0
         numgrid(3) = fcstgrd(3) - 1
         do n = 1, fnstgrd
            if (mod(numgrid(3), 2) .eq. 0) then
               numgrid(3) = numgrid(3)/2
               nm = nm + 1
            elseif (nm .eq. 0) then
               print *, 'currently, the number of analysis vertical levels must be odd!'
               stop
            else
               exit      ! use current numgrid(3) to start multigrid
            end if
         end do
         maxgrid(3) = numgrid(3)*2**(nm) + 1
         numgrid(3) = numgrid(3) + 1

         ! for t direction:
         numgrid(4) = 2
         maxgrid(4) = 3

      else

         ! using maxgrid to setup multigrid:
         do n = 1, numdims
            if (maxgrid(n) .gt. 1 .and. numgrid(n) .gt. 1) then
               if (mod(maxgrid(n) - 1, numgrid(n) - 1) .eq. 0) then
                  nm = (maxgrid(n) - 1)/(numgrid(n) - 1)
                  ns = 1
                  do while (nm .ge. 2)
                     if (mod(nm, 2) .eq. 0) then
                        nm = nm/2
                        ns = ns + 1
                     else
                        print *, 'maxgrid should be (numgrid-1)*2**n+1'
                        stop
                     end if
                  end do
               else
                  print *, 'maxgrid should be (numgrid-1)*2**n+1'
                  stop
               end if
               if (ns .gt. fnstgrd) then
                  maxgrid(n) = (numgrid(n) - 1)*2**(fnstgrd - 1) + 1
               end if
            end if
         end do

      end if ! end multigrid setup

      ! maxgrid in time is the same as final analysis:
      fcstgrd(4) = maxgrid(4)

      inigrid = numgrid          ! save initial start multigrid numbers

      ! initial vertical temporarl gridspacing
      grdspac(3:4) = 0.0
      if (maxgrid(3) .gt. 1) grdspac(3) = (maxgrid(3) - 1)/float(numgrid(3) - 1)

      print *, ''
      print *, 'stmas namelist has been read with'
      write (*, 1) numgrid, maxgrid
1     format(' numgrid: ', 4i4, ' maxgrid: ', 4i4)
      print *, ''

   end subroutine read_namelist

! include 'laps_configs.f90'
   subroutine laps_config

!********************************************************************
!  read in laps configuration parameters.
!
!  history: jan. 2009 by yuanfu xie.
!
!           modified dec. 2013 by yuanfu xie for p_sfc_f used for tpw
!           calculation
!********************************************************************

      use prmtrs_stmas

      implicit none

      integer :: n, istate

      ! system time:
      call get_systime(lapsi4t, lapsast, istate)

      ! laps cycle time:
      call get_laps_cycle_time(icycle, istate)

      ! i4time of jan 1, 1970
      call cv_asc_i4time('700010000', i4t_gps)

      ! missing value in real:
      call get_r_missing_data(rmissing, istate)
      ! bad surface data flag:
      call get_sfc_badflag(badsfcdt, istate)

      ! spatial x y dimensions:
      call get_grid_dim_xy(fcstgrd(1), fcstgrd(2), istate)
      call get_laps_dimensions(fcstgrd(3), istate)
      print *, (fcstgrd(n), n=1, numdims)

      ! laps array:
      grdspac(1) = ((fcstgrd(1) - 1)*1.0)/((numgrid(1) - 1)*1.0)
      grdspac(2) = ((fcstgrd(2) - 1)*1.0)/((numgrid(2) - 1)*1.0)
      if (numgrid(4) .ge. 2) then
         grdspac(4) = (fcstgrd(4) - 1.0)/(numgrid(4) - 1.0)
      else
         grdspac(4) = 0.0
      end if
      allocate (gridmask(fcstgrd(1), fcstgrd(2), fcstgrd(3), fcstgrd(4), numstat), stat=istate)
      if (istate .ne. 0) stop 'gridmask allocate wrong'

      allocate (z_fcstgd(fcstgrd(3)), stat=istate)
      if (istate .ne. 0) stop 'z_fcstgd allocate wrong'
      allocate (z_maxgid(maxgrid(3)), stat=istate)
      if (istate .ne. 0) stop 'z_maxgid allocate wrong'

      call alloctbkg
      call allocatob

      ! laps lat/lon/topography:
      call read_static_grid(fcstgrd(1), fcstgrd(2), 'lat', latitude, istate)
      if (istate .ne. 1) then
         write (6, *) 'rdlapsrdr: error get laps lat'
         stop
      end if
      call read_static_grid(fcstgrd(1), fcstgrd(2), 'lon', longitud, istate)
      if (istate .ne. 1) then
         write (6, *) 'rdlapsrdr: error get laps lon'
         stop
      end if
      call read_static_grid(fcstgrd(1), fcstgrd(2), 'avg', topogrph, istate)
      if (istate .ne. 1) then
         write (6, *) 'rdlapsrdr: error get laps avg'
         stop
      end if

      ! allocate memory for surface pressure:
      allocate (p_sfc_f(fcstgrd(1), fcstgrd(2)), stat=istate)

   end subroutine laps_config

   subroutine bkgrndobs
!*************************************************
! main routine of this preparation code
! history: september 2007, coded by wei li.
!*************************************************
      implicit none

      integer      :: istate

      print *, 'readnmlst'
      ! call readnmlst
      call read_namelist

      print *, 'bkgmemalc'
      call bkgmemalc
      print *, 'rdbckgrnd'
      if (if_test .ne. 1) then
         ! laps configuration:
         call laps_config
         !  call rdbckgrnd
         call rdlapsbkg
      else
         open (11, file='fort.11', status='old', action='read')
         read (11, *) fcstgrd(1:numdims)
         close (11)
         print *, 'memoryalc'
         call memoryalc
         call rdbkgtest
      end if
      print *, 'getbkgrnd_new'
      call getbkgrnd_new
      print *, 'rdbufrobs'
      if (if_test .ne. 1) then
         ! call open_lapsprd_file(tmgobs_channel,lapsi4t,'tmg',istate)
         ! call open_lapsprd_file(pigobs_channel,lapsi4t,'pig',istate)
         call rdbufrobs_xie
         ! close(tmgobs_channel)
         ! close(pigobs_channel)
         !  call rdradrobs
         ! call rdlapsrdr
         call read_laps_radar  ! switch to a new routine using less memory by yuanfu
         ! call gpswdelay
      else
         call rdobstest
      end if
      if (ifbkgnd .eq. 1) then
         print *, 'addbkgrnd'
         call addbkgrnd
      end if
      print *, 'obsmemalc'
      call obsmemalc
      print *, 'memoryrls'
      call memoryrls
      return
   end subroutine bkgrndobs

   subroutine readnmlst
!*************************************************
! read name list
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: s, n, er, nu, ln, nm, ns
      character(len=200) :: dr
      character(len=256) :: st

      ! get namelist from static:
      call get_directory('static', st, n)
! replacing the file name from stmas3d.nl to stmas3d.txt, and also increase the number by 1. n+10 to n+11. hj 6/28/2011
!  dr = st(1:n)//'stmas3d.nl'
      dr = st(1:n)//'stmas3d.txt'
! --------------------
      nu = 2
      !open(nu,file='namelist.txt',action='read',status='old')
! replacing the file name from stmas3d.nl to stmas3d.txt, and also increase the number by 1. n+10 to n+11. hj 6/28/2011
!  open(nu,file=dr(1:n+10),action='read',status='old')
      open (nu, file=dr(1:n + 11), action='read', status='old')

      read (nu, *) ifbkgnd
      read (nu, *) ifbound
      read (nu, *) ifpcdnt

      if (ifpcdnt .eq. 1) then
         print *, 'ifpcdnt=', ifpcdnt, 'for presure coordinate'
      elseif (ifpcdnt .eq. 2) then
         print *, 'ifpcdnt=', ifpcdnt, 'for height coordinate'
      elseif (ifpcdnt .eq. 0) then
         print *, 'ifpcdnt=', ifpcdnt, 'for sigma coordinate'
      end if

      read (nu, *) numstat
      allocate (sl0(numstat), stat=er)
      if (er .ne. 0) stop 'sl0 allocate wrong'
      do s = 1, numstat
         read (nu, *) sl0(s)
      end do
      allocate (penal0x(numstat), stat=er)
      if (er .ne. 0) stop 'penal0x allocate wrong'
      allocate (penal0y(numstat), stat=er)
      if (er .ne. 0) stop 'penal0y allocate wrong'
      allocate (penal0z(numstat), stat=er)
      if (er .ne. 0) stop 'penal0z allocate wrong'
      allocate (penal0t(numstat), stat=er)
      if (er .ne. 0) stop 'penal0t allocate wrong'
      allocate (penal_x(numstat), stat=er)
      if (er .ne. 0) stop 'penal_x allocate wrong'
      allocate (penal_y(numstat), stat=er)
      if (er .ne. 0) stop 'penal_y allocate wrong'
      allocate (penal_z(numstat), stat=er)
      if (er .ne. 0) stop 'penal_z allocate wrong'
      allocate (penal_t(numstat), stat=er)
      if (er .ne. 0) stop 'penal_t allocate wrong'
      do s = 1, numstat
         read (nu, *) penal0x(s)
         read (nu, *) penal0y(s)
         read (nu, *) penal0z(s)
         read (nu, *) penal0t(s)
         print *, 'smoothing var: ', s, ' with: ', penal0x(s), penal0y(s), penal0z(s), penal0t(s)
      end do
      read (nu, *) pnlt0pu
      read (nu, *) pnlt0pv
      read (nu, *) numdims
      read (nu, *) numgrid(1)
      read (nu, *) numgrid(2)
      read (nu, *) numgrid(3)
      read (nu, *) numgrid(4)
      read (nu, *) fcstgrd(4)                ! for laps ingest: read into laps time frame for temporal analysis yuanfu
      read (nu, *) grdspac(4)
      read (nu, *) oripstn(1)
      read (nu, *) oripstn(2)
      read (nu, *) oripstn(3)
      read (nu, *) oripstn(4)
      read (nu, *) maxgrid(1)
      read (nu, *) maxgrid(2)
      read (nu, *) maxgrid(3)
      read (nu, *) maxgrid(4)
      read (nu, *) fnstgrd
      allocate (obradius(maxdims, numstat), stat=er)
      if (er .ne. 0) stop 'obradius allocate wrong'
      do s = 1, numstat
         do n = 1, maxdims
            read (nu, *) obradius(n, s)
         end do
      end do
      read (nu, *) u_cmpnnt
      read (nu, *) v_cmpnnt
      read (nu, *) w_cmpnnt
      read (nu, *) pressure         ! 'pressure' is for z coordinate, for presure coordinate it is height
      read (nu, *) temprtur
      read (nu, *) humidity
      if (numstat .gt. 5) &
         read (nu, *) rour_cmpnnt  ! added by shuyuan20100721 for density(rou)*rain water mixing ratio
      if (numstat .gt. 6) &
         read (nu, *) rous_cmpnnt  ! added by shuyuan20100721 for density(rou)*snow water mixing ratio
      read (nu, *) xsl
      read (nu, *) ysl
      read (nu, *) psl
      read (nu, *) csl
      read (nu, *) dsl
      read (nu, *) cosstep
      read (nu, *) midgrid
      read (nu, *) finstep

      read (nu, *) xytrans             ! coefficient used to translate the x and y coordinate to meters
      read (nu, *) z_trans             ! coefficient used to translate the z coordinate to meters
      read (nu, *) if_test             ! whether run for test case, 1 is for the test case
      read (nu, *) ifrepet             ! whethere run the mutigrid frame in monotonous or repeatedly, 0 for monotonous, 1 for repeatedly
      read (nu, *) itrepet             ! the times to repeat
      if (ifrepet .eq. 0) itrepet = 0

      read (nu, *) pnlt0hy             ! initial hydrostatic condition penalty coefficent
      read (nu, *) taul_hy             ! reduction coefficent of hydrostatic condition penalty term
      read (nu, *) endhylv             ! the grid level after which the hydrostatic condition penalty term is ommitted
      read (nu, *) endgslv             ! the grid level after which the geostrophic balance penalty term is ommitted

      print *, 'hhh: ', humidity, midgrid, endgslv, endhylv
      read (nu, *) limit_3             ! the limitation of height or pressure to decide in which range the observation is aviable.
      read (nu, *) limit_4             ! the limitation of time to decide in which range the observation is aviable.

      close (nu)
      ngptobs = 2**numdims
      nallobs = 0                     ! initializing the number of observation, addied by zhongjie he
      grdspac(3) = 0.
      if (maxgrid(3) .ne. 1) grdspac(3) = (maxgrid(3) - 1.0)/(numgrid(3) - 1.0)      ! by zhongjie he

      do n = 1, numdims
         if (maxgrid(n) .gt. 1 .and. numgrid(n) .gt. 1) then
            if (mod(maxgrid(n) - 1, numgrid(n) - 1) .eq. 0) then
               nm = (maxgrid(n) - 1)/(numgrid(n) - 1)
               ns = 1
               do while (nm .ge. 2)
                  if (mod(nm, 2) .eq. 0) then
                     nm = nm/2
                     ns = ns + 1
                  else
                     print *, 'maxgrid should be (numgrid-1)*2**n+1'
                     stop
                  end if
               end do
            else
               print *, 'maxgrid should be (numgrid-1)*2**n+1'
               stop
            end if
            if (ns .gt. fnstgrd) then
               maxgrid(n) = (numgrid(n) - 1)*2**(fnstgrd - 1) + 1
            end if
         end if
      end do

      do n = 1, numdims
         inigrid(n) = numgrid(n)
      end do

      return
   end subroutine readnmlst

   subroutine bkgmemalc
!*************************************************
! allocate memory for background array
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: er
! --------------------
      ! two additional spaces of grdbkgd0 for saving lower and upper bounds of sh:
      allocate (grdbkgd0(maxgrid(1), maxgrid(2), maxgrid(3), maxgrid(4), numstat + 2), stat=er)
      if (er .ne. 0) stop 'grdbkgd0 allocate wrong'
      allocate (xx0(maxgrid(1), maxgrid(2)), stat=er)
      if (er .ne. 0) stop 'xx0 allocate wrong'
      allocate (yy0(maxgrid(1), maxgrid(2)), stat=er)
      if (er .ne. 0) stop 'yy0 allocate wrong'
      allocate (cr0(maxgrid(1), maxgrid(2)), stat=er)
      if (er .ne. 0) stop 'cr0 allocate wrong'
      allocate (dg0(maxgrid(1), maxgrid(2)), stat=er)
      if (er .ne. 0) stop 'dg0 allocate wrong'
      allocate (dn0(maxgrid(1), maxgrid(2), maxgrid(3), maxgrid(4)), stat=er)
      if (er .ne. 0) stop 'dn0 allocate wrong'
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then
         allocate (zz0(maxgrid(1), maxgrid(2), maxgrid(3), maxgrid(4)), stat=er)
         if (er .ne. 0) stop 'zz0 allocate wrong'
         allocate (zzb(maxgrid(3)), stat=er)
         if (er .ne. 0) stop 'zzb allocate wrong'
      elseif (ifpcdnt .eq. 1) then
         allocate (pp0(maxgrid(3)), ppm(maxgrid(3)), stat=er)
         if (er .ne. 0) stop 'pp0 and ppm allocate wrong'
      end if
      return
   end subroutine bkgmemalc

   subroutine obsmemalc
!*************************************************
! allocate observation memory and read in data and scale
! history: august 2007, coded by wei li.
!        : march 2008, modified by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: n, o, er, ln, s, no
      character(len=200) :: dr
! --------------------
      if (nallobs .eq. 0) return
      allocate (obspostn(numdims, nallobs), stat=er)
      if (er .ne. 0) stop 'obspostn allocate wrong'
      allocate (obscoeff(ngptobs, nallobs), stat=er)
      if (er .ne. 0) stop 'obscoeff allocate wrong'
      allocate (obsidxpc(maxdims, nallobs), stat=er)
      if (er .ne. 0) stop 'obsidxpc allocate wrong'
      allocate (obsvalue(nallobs), stat=er)
      if (er .ne. 0) stop 'obsvalue allocate wrong'
      allocate (obserror(nallobs), stat=er)
      if (er .ne. 0) stop 'obserror allocate wrong'
!jhui
      allocate (nobstat(numstat + 3), stat=er)
      if (er .ne. 0) stop 'nobstat allocate wrong'
      allocate (obseinf1(nallobs), stat=er)
      if (er .ne. 0) stop 'obseinf1 allocate wrong'
      allocate (obseinf2(nallobs), stat=er)
      if (er .ne. 0) stop 'obseinf2 allocate wrong'
!  allocate(obseinf3(nallobs),stat=er)
!  if(er.ne.0)stop 'obseinf3 allocate wrong'

      o = 0
      do s = 1, numstat + 3
         do no = 1, nst(s)
            o = o + 1
            do n = 1, numdims
               obspostn(n, o) = obp(n, no, s)
            end do
            obsvalue(o) = obs(no, s)
            obserror(o) = obe(no, s)
         end do
         nobstat(s) = nst(s)
      end do
      s = numstat + 1
      do no = 1, nst(s)
         obseinf1(no) = oba(no, 1)
         obseinf2(no) = oba(no, 2)
      end do
!jhui
!  s=numstat+3
!  do no=1,nst(s)
!   write(306,*) obs(no,s)
!  enddo

      return
   end subroutine obsmemalc

   subroutine memoryalc
!*************************************************
! memory allocate
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: n, er
      integer :: st                ! status
! --------------------

      grdspac(1) = ((fcstgrd(1) - 1)*1.0)/((numgrid(1) - 1)*1.0)
      grdspac(2) = ((fcstgrd(2) - 1)*1.0)/((numgrid(2) - 1)*1.0)
      if (numgrid(4) .ge. 2) then
         grdspac(4) = (fcstgrd(4) - 1.0)/(numgrid(4) - 1.0)
      else
         grdspac(4) = 0.0
      end if
      allocate (gridmask(fcstgrd(1), fcstgrd(2), maxgrid(3), fcstgrd(4), numstat), stat=er) !!!!!attention
      if (er .ne. 0) stop 'gridmask allocate wrong'

      allocate (z_fcstgd(fcstgrd(3)), stat=er)
      if (er .ne. 0) stop 'z_fcstgd allocate wrong'
      allocate (z_maxgid(maxgrid(3)), stat=er)
      if (er .ne. 0) stop 'z_maxgid allocate wrong'

      call alloctbkg
      call allocatob
      return
   end subroutine memoryalc

   subroutine memoryrls
!*************************************************
! memory release
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
      deallocate (gridmask)
      deallocate (obradius)
      call dealctbkg
      call dealoctob
      return
   end subroutine memoryrls

   subroutine getbkgrnd_new
!*************************************************
! get data of background for analysis, adapt to pressure height and sigma coordinate.
! history: september 2007, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, ln, nm
      real     :: xb(maxgrid(1)), yb(maxgrid(2)), zb(maxgrid(3)), tb(maxgrid(4))
      real     :: xf(fcstgrd(1)), yf(fcstgrd(2)), zf(fcstgrd(3)), tf(fcstgrd(4))
      integer  :: fg(maxdims), mg(maxdims)
      real     :: z1(1), t1(1), z2(1), t2(1), dx, dy, dt
      character(len=200) :: dr

      real     :: sh          ! coefficient to translate the zz0 to length with unit of meters.
      logical  :: lexist      ! used to decide whether the level data file is exist
      integer :: clock_count0, clock_rate, clock_max, clock_count1

      integer :: init_timer, ishow_timer
! --------------------

      sh = z_trans

      do i = 1, fcstgrd(1)
         xf(i) = (i - 1)*1.0d0
      end do
      do j = 1, fcstgrd(2)
         yf(j) = (j - 1)*1.0d0
      end do
      dx = ((fcstgrd(1) - 1)*1.0d0)/((maxgrid(1) - 1)*1.0d0)
      do i = 1, maxgrid(1)
         xb(i) = xf(1) + (i - 1)*dx
      end do
      dy = ((fcstgrd(2) - 1)*1.0d0)/((maxgrid(2) - 1)*1.0d0)
      do j = 1, maxgrid(2)
         yb(j) = yf(1) + (j - 1)*dy
      end do

      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
         do k = 1, fcstgrd(3)
            zf(k) = z00(k)
         end do
      else                                         ! for pressure coordinate
         do k = 1, fcstgrd(3)
            zf(k) = p00(k)
         end do
      end if

!  open(2,file='p_level.dat',status='old',action='read')
!    do k=1,maxgrid(3)
!      read(2,*)zb(k)
!    enddo
!  close(2)
      if (maxgrid(3) .lt. fcstgrd(3)) then
         print *, 'read the vertical levels from data file: level.dat'
         inquire (file='level.dat', exist=lexist)
         if (lexist) then
            open (2, file='level.dat', status='old', action='read')
            do k = 1, maxgrid(3)
               read (2, *) zb(k)
            end do
            close (2)
         else
            print *, 'the file level.dat does not exist'
            stop
         end if
      elseif (maxgrid(3) .gt. 2*fcstgrd(3) - 1) then
         print *, 'error! maxgrid(3) should be smaller than 2*fcstgrd(3)-1!'
         stop
      elseif (maxgrid(3) .eq. fcstgrd(3)) then
         do k = 1, maxgrid(3)
            zb(k) = zf(k)
         end do
      else
         nm = maxgrid(3) - fcstgrd(3)
         do k = 1, nm
            zb(2*k - 1) = zf(k)
            zb(2*k) = 0.5*(zf(k) + zf(k + 1))
         end do
         do k = nm + 1, fcstgrd(3)
            zb(k + nm) = zf(k)
         end do
!======
!    do k=1,fcstgrd(3)-1
!      zb(k)=zf(k)
!    enddo
!    do k=fcstgrd(3),maxgrid(3)
!      zb(k)=zf(fcstgrd(3)-1)+(zf(fcstgrd(3))-zf(fcstgrd(3)-1))/float(maxgrid(3)-fcstgrd(3)+1)*(k-fcstgrd(3)+1)
!    enddo
!======
      end if

!  do k=1,maxgrid(3)
!    zb(k)=zf(1)+(zf(fcstgrd(3))-zf(1))/float(maxgrid(3)-1)*(k-1)
!  enddo

!===============
      do k = 1, maxgrid(3)
         z_maxgid(k) = zb(k)
      end do
!===============

      do k = 1, maxgrid(3)
         if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
            zzb(k) = zb(k)
         else                                          ! for pressure coordinate
            pp0(k) = zb(k)
         end if
      end do
      ! multigrid pressure coordinate:
      i = (maxgrid(3) - 1)/(numgrid(3) - 1)
      do k = 1, numgrid(3)
         ppm(k) = pp0((k - 1)*i + 1)
      end do

      if (ifpcdnt .eq. 0) then                ! for sigma coordinate
         do t = 1, maxgrid(4)
         do j = 1, maxgrid(2)
         do i = 1, maxgrid(1)
         do k = 1, maxgrid(3)
            zz0(i, j, k, t) = (zb(k)*(heightu(i, j) - heightl(i, j)) + heightl(i, j))*sh   ! translate the unit to meters
         end do
         end do
         end do
         end do
      elseif (ifpcdnt .eq. 2) then            ! for height coordinate
         do t = 1, maxgrid(4)
         do j = 1, maxgrid(2)
         do i = 1, maxgrid(1)
         do k = 1, maxgrid(3)
            zz0(i, j, k, t) = zb(k)
         end do
         end do
         end do
         end do
      end if

!=======================================

      do t = 1, fcstgrd(4)
         tf(t) = (t - 1)*1.0
      end do
      if (maxgrid(4) .ge. 2) dt = ((fcstgrd(4) - 1)*1.0)/((maxgrid(4) - 1)*1.0)
      do t = 1, maxgrid(4)
         tb(t) = tf(1) + (t - 1)*dt
      end do

!  call system_clock(clock_count0,clock_rate,clock_max)
!  call fcst2bkgd(numdims,ngptobs,numstat,fcstgrd, &
!                 xf,yf,zf,tf,maxgrid,xb,yb,zb,tb,bk0,grdbkgd0)

      print *, 'starting interpolating background onto maxgrid: ', init_timer()
      if (uniform .eq. 0) then
         call bkgtofine(numstat, fcstgrd, xf, yf, zf, tf, maxgrid, xb, yb, zb, tb, bk0, grdbkgd0)
      else ! use yuanfu's uniform interpolation routine:
         do s = 1, numstat
         do t = 1, fcstgrd(4)
            call uniform_interpolation3(fcstgrd, maxgrid, bk0(:, :, :, t, s), grdbkgd0(:, :, :, t, s))
         end do
         end do
      end if
      print *, 'ending : interpolating background onto maxgrid: ', ishow_timer()

!  call system_clock(clock_count1,clock_rate,clock_max)
!  write ( *, '(a)' ) '  system_clock with integer arguments reports: in subroutine of <fcst2bkgd> '
!  write ( *, '(a,i12)' ) '    the current clock count is    ', clock_count1-clock_count0
!  write ( *, '(a,i12)' ) '    the clock count per second is ', clock_rate
!  write ( *, '(a,i12)' ) '    the maximum clock count is    ', clock_max

      ! yuanfu remove s00 array as it was set to 1 and it was passed to dn0
      ! call bkgtofine(1,fcstgrd,xf,yf,zf,tf,maxgrid,xb,yb,zb,tb,s00,dn0)
      dn0 = 1.0 ! density on the multigrid. it is supposed to use real density but constant for now.

      fg(1) = fcstgrd(1)
      fg(2) = fcstgrd(2)
      fg(3) = 1
      fg(4) = 1
      z1(1) = 0.0
      t1(1) = 0.0
      mg(1) = maxgrid(1)
      mg(2) = maxgrid(2)
      mg(3) = 1
      mg(4) = 1
      z2(1) = 0.0
      t2(1) = 0.0

!  call fcst2bkgd(2,4,1,fg,xf,yf,z1,t1,mg,xb,yb,z2,t2,c00,cr0)
      if (uniform .eq. 0) then
         call bkgtofine(1, fg, xf, yf, z1, t1, mg, xb, yb, z2, t2, c00, cr0)
      else
         call uniform_interpolation2(fg, mg, c00, cr0)
      end if

!  call fcst2bkgd(2,4,1,fg,xf,yf,z1,t1,mg,xb,yb,z2,t2,d00,dg0)
      ! yuanfu skips this call for dg0 since d00 is set to zero, a constant in read_backgrd:
      ! call bkgtofine(1,fg,xf,yf,z1,t1,mg,xb,yb,z2,t2,d00,dg0)
      dg0 = 0.0

      dx = ((fcstgrd(1) - 1)*dx0)/((maxgrid(1) - 1)*1.0d0)
      dy = ((fcstgrd(2) - 1)*dy0)/((maxgrid(2) - 1)*1.0d0)
      if (maxgrid(4) .ge. 2) dt = ((fcstgrd(4) - 1)*dt0)/((maxgrid(4) - 1)*1.0)
      do i = 1, maxgrid(1)
      do j = 1, maxgrid(2)
         xx0(i, j) = 0.0d0 + (i - 1)*dx
         yy0(i, j) = 0.0d0 + (j - 1)*dy
      end do
      end do

      return
   end subroutine getbkgrnd_new

   subroutine addbkgrnd
!*************************************************
! set background as observations at the grid points where observations can not affect
! history: october 2007, coded by wei li.
!          march 2008, modified by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, n, rx, ry, rt, np(maxdims)
      real     :: rh, rz, oe(numstat), oc(numdims), p, pp(maxgrid(3))
      integer  :: ss, os                                             ! added by zhongjie he
! --------------------
      if (nallobs .eq. 0) return
      do s = 1, numstat
         do t = 1, fcstgrd(4)
         do k = 1, maxgrid(3)      !!!!!!!attention, due to special character of vertical coordinate
         do j = 1, fcstgrd(2)
         do i = 1, fcstgrd(1)
            gridmask(i, j, k, t, s) = 1
         end do
         end do
         end do
         end do
      end do
      do s = 1, numstat
         oe(s) = 0.0d0
      end do
      do s = 1, numstat
         do o = 1, nst(s)
            do n = 1, numdims
               oc(n) = obp(n, o, s) + 1.0
            end do
            if (ifpcdnt .eq. 1) then       ! for pressure coordinate
               do k = 1, maxgrid(3)
                  pp(k) = pp0(k)
               end do
               if (maxgrid(3) .ge. 2) then
                  k = int(oc(3))
                  p = (oc(3) - k)*(pp(k + 1) - pp(k)) + pp(k)
               end if
            else                      ! for sigma and height coordinate
               do k = 1, maxgrid(3)
                  pp(k) = z00(k)
               end do
               if (maxgrid(3) .ge. 2) then
                  k = int(oc(3))
                  p = (oc(3) - k)*(pp(k + 1) - pp(k)) + pp(k)
               end if
            end if
            rx = 0
            ry = 0
            rt = 0
            if (fcstgrd(1) .ge. 2) rx = int(obradius(1, s)/dx0) + 1
            if (fcstgrd(2) .ge. 2) ry = int(obradius(2, s)/dy0) + 1
            if (fcstgrd(4) .ge. 2) rt = int(obradius(4, s)/dt0) + 1
            do t = max0(int(oc(4)) - rt, 1), min0(int(oc(4)) + rt + 1, fcstgrd(4))
            do k = 1, maxgrid(3)      !!!!!!!attention, due to special character of vertical coordinate
            do j = max0(int(oc(2)) - ry, 1), min0(int(oc(2)) + ry + 1, fcstgrd(2))
            do i = max0(int(oc(1)) - rx, 1), min0(int(oc(1)) + rx + 1, fcstgrd(1))
               rh = sqrt(((oc(1) - i)*dx0)**2 + ((oc(2) - j)*dy0)**2)
               if (ifpcdnt .eq. 1 .or. ifpcdnt .eq. 2) then       ! for pressure and height coordinate
                  rz = abs((pp(k) - p))
               elseif (ifpcdnt .eq. 0) then                     ! for sigma coordinate
                  rz = abs((pp(k) - p))*(heightu(i, j) - heightl(i, j))
               end if
               if (rh .le. obradius(1, s) .and. rz .le. obradius(3, s)) gridmask(i, j, k, t, s) = 0
            end do
            end do
            end do
            end do
            oe(s) = max(oe(s), obe(o, s)*1.0)
         end do
      end do

      o = nallobs
      do s = 1, numstat
         do t = 1, fcstgrd(4)
         do k = 1, maxgrid(3)      !!!!!!!attention, due to special character of vertical coordinate
         do j = 1, fcstgrd(2), 3
         do i = 1, fcstgrd(1), 3
            if (gridmask(i, j, k, t, s) .eq. 1 .and. nst(s) .ge. 1) then
               np(1) = i - 1
               np(2) = j - 1
               np(3) = k - 1
               np(4) = t - 1
               o = o + 1
               nst(s) = nst(s) + 1
               if (o .gt. nobsmax) stop 'number of observations exceeded'
               do n = 1, numdims
                  obp(n, nst(s), s) = np(n)
               end do
               obs(nst(s), s) = 0.0d0
               obe(nst(s), s) = oe(s)
            end if
         end do
         end do
         end do
         end do
      end do
      nallobs = o

      if (.true.) then                      ! the following is added by zhongjie he
         do os = numstat + 1, numstat + 2
         do o = 1, nst(os)
            do n = 1, numdims
               oc(n) = obp(n, o, os) + 1.0
            end do
            if (ifpcdnt .eq. 1) then        ! for pressure coordinate
               do k = 1, maxgrid(3)
                  pp(k) = pp0(k)
               end do
               if (maxgrid(3) .ge. 2) then
                  k = int(oc(3))
                  p = (oc(3) - k)*(pp(k + 1) - pp(k)) + pp(k)
               end if
            else                       ! for sigma and height coordinate
               do k = 1, maxgrid(3)
                  pp(k) = z00(k)
               end do
               if (maxgrid(3) .ge. 2) then
                  k = int(oc(3))
                  p = (oc(3) - k)*(pp(k + 1) - pp(k)) + pp(k)
               end if
            end if
            rx = 0
            ry = 0
            rt = 0
            if (fcstgrd(1) .ge. 2) rx = int(max(obradius(1, u_cmpnnt), obradius(1, v_cmpnnt))/dx0) + 1
            if (fcstgrd(2) .ge. 2) ry = int(max(obradius(2, u_cmpnnt), obradius(2, v_cmpnnt))/dy0) + 1
            if (fcstgrd(4) .ge. 2) rt = int(max(obradius(4, u_cmpnnt), obradius(4, v_cmpnnt))/dt0) + 1
            do t = max0(int(oc(4)) - rt, 1), min0(int(oc(4)) + rt + 1, fcstgrd(4))
            do k = 1, maxgrid(3)      !!!!!!!attention, due to special character of vertical coordinate
            do j = max0(int(oc(2)) - ry, 1), min0(int(oc(2)) + ry + 1, fcstgrd(2))
            do i = max0(int(oc(1)) - rx, 1), min0(int(oc(1)) + rx + 1, fcstgrd(1))
               rh = sqrt(((oc(1) - i)*dx0)**2 + ((oc(2) - j)*dy0)**2)
               if (ifpcdnt .eq. 1 .or. ifpcdnt .eq. 2) then       ! for pressure and height coordinate
                  rz = abs((pp(k) - p))
               elseif (ifpcdnt .eq. 0) then                     ! for sigma coordinate
                  rz = abs((pp(k) - p))*(heightu(i, j) - heightl(i, j))
               end if
               do ss = 1, 2
                  if (ss .eq. 1) s = u_cmpnnt
                  if (ss .eq. 2) s = v_cmpnnt
                  if (rh .le. obradius(1, s) .and. rz .le. obradius(3, s)) gridmask(i, j, k, t, s) = 0
               end do
            end do
            end do
            end do
            end do
            oe(u_cmpnnt) = max(oe(u_cmpnnt), obe(o, os)*1.0)
            oe(v_cmpnnt) = max(oe(v_cmpnnt), obe(o, os)*1.0)
         end do
         end do

         o = nallobs
         if (nst(numstat + 1) + nst(numstat + 2) .ge. 1) then
            do t = 1, fcstgrd(4)
            do k = 1, maxgrid(3)      !!!!!!!attention, due to special character of vertical coordinate
            do j = 1, fcstgrd(2), 3
            do i = 1, fcstgrd(1), 3
               do ss = 1, 2
                  if (ss .eq. 1) s = u_cmpnnt
                  if (ss .eq. 2) s = v_cmpnnt
                  if (gridmask(i, j, k, t, s) .eq. 1) then
                     np(1) = i - 1
                     np(2) = j - 1
                     np(3) = k - 1
                     np(4) = t - 1
                     o = o + 1
                     nst(s) = nst(s) + 1
                     if (o .gt. nobsmax) stop 'number of observations exceeded'
                     do n = 1, numdims
                        obp(n, nst(s), ss) = np(n)
                     end do
                     obs(nst(s), s) = 0.0d0
                     obe(nst(s), s) = oe(s)
                  end if
               end do
            end do
            end do
            end do
            end do
         end if
         nallobs = o
      end if
      return

      print *, 'after adding background, nallobs=', nallobs

   end subroutine addbkgrnd

end module input_bg_obs
