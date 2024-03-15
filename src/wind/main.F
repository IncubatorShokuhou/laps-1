cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis

      program laps_wind_driver

      use mem_namelist, only: read_namelist_laps
      use mem_namelist, only: nx_l,ny_l,nk_laps,r_missing_data
     1                         ,iwrite_output

      use mem_grid,only:lat,lon,topo
#ifdef usempi
     1                 ,npes,rank,istartiend,recvcounts,displs
#endif

!      use mem_wind, wind%uanl->uanl,wind%vanl->vanl,wind%wanl->wanl
!     +          , wind%uanl_sfcitrp->uanl_sfcitrp
!     +          , wind%vanl_sfcitrp->vanl_sfcitrp
      use mem_wind

      implicit none

!     general declarations
      integer n_3d_fields,len_dir,ntmin,ntmax,n_meso,n_sao,n_pirep
      integer istatus,istatus_alloc,istat_alloc,istat_ht,istat_wind
      integer istat_lw3
      parameter (n_3d_fields = 3)

      character*3 var_3d(n_3d_fields) 
      character*9 a9_time
      character*3 var_2d
      character*31 ext

!     laps grid dimensions

      real, allocatable, dimension(:,:,:) :: heights_3d
      !real, allocatable, dimension(:,:,:) :: uanl
      !real, allocatable, dimension(:,:,:) :: vanl
      !real, allocatable, dimension(:,:) :: uanl_sfcitrp
      !real, allocatable, dimension(:,:) :: vanl_sfcitrp

!     wind parameters

      character*150 static_dir,filename
!     logical l_use_raob, l_use_cdw, l_use_radial_vel

      integer i4time_sys

#ifdef usempi
!     mpi stuff
      include 'mpif.h'
      integer :: ierr,n,flag,count
      integer :: status(mpi_status_size)
      integer :: pointsperpe,extra,adjust

      integer :: max_obs,imax,jmax,kmax,ncnt_total
      integer :: idelt,jdelt,kdelt
      real    :: r_missing_data_local,r0_vert_grid

      integer,allocatable :: ibuf(:)
      real   ,allocatable :: rbuf(:)
      real   ,allocatable,dimension(:,:,:,:) :: wt_b_local,local_var
      real   ,allocatable :: fnorm_lut(:,:)
      real :: buffer !argument to mpi_gatherv ignored on non-root

      include 'barnesob.inc'
      include 'windparms.inc'

      type (barnesob),allocatable :: obs_barnes(:)

!     set up the mpi tasks and spit process 0 from the rest
      call mpi_init( ierr )
      if( ierr /= 0) then
        print*,'mpi_init error',ierr
        stop
      endif
      call mpi_comm_rank( mpi_comm_world, rank, ierr )
      if( ierr /= 0) then
        print*,'mpi_comm_rank error',ierr
        stop
      endif
      call mpi_comm_size( mpi_comm_world, npes, ierr )
      if( ierr /= 0) then
        print*,'mpi_comm_size error',ierr,rank
        stop
      endif
      print *, 'process ', rank, ' of ', npes, ' is alive'
      allocate(istartiend(2,0:npes-1), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' error: could not allocate istartiend',istat_alloc,rank
        stop
      endif
      allocate(recvcounts(  0:npes-1), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' error: could not allocate recvcounts',istat_alloc,rank
        stop
      endif
      allocate(displs    (  0:npes-1), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
        print*,' error: could not allocate displs',istat_alloc,rank
        stop
      endif
      if(rank > 0) then
       call mpi_bcast(istartiend,2*npes,mpi_integer,0,mpi_comm_world,
     1                                                          ierr)
       if( ierr /= 0) then
         print*,'mpi_bcast istartiend error',ierr,rank
         stop
       endif
       print*,'rank,istart,iend',rank,istartiend(1,rank)
     1                               ,istartiend(2,rank)
       call mpi_bcast(displs,npes,mpi_integer,0,mpi_comm_world,ierr)
       if( ierr /= 0) then
         print*,'mpi_bcast displs error',ierr,rank
         stop
       endif
       call mpi_bcast(recvcounts,npes,mpi_integer,0,mpi_comm_world,ierr)
       if( ierr /= 0) then
         print*,'mpi_bcast recvcounts error',ierr,rank
         stop
       endif
       print*,'process',rank,' is waiting'
       flag = 1
       do while(flag==1) 
         call mpi_bcast(flag,1,mpi_integer,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_bcast flag error',ierr,rank
           stop
         endif
         if(flag==-1) exit
         print*,'process',rank,' received go-ahead',flag,ierr
!        receive and unpack integer buffers
         allocate(ibuf(10), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' error: could not allocate ibuf(10)',istat_alloc,rank
           stop
         endif
         call mpi_bcast(ibuf,10,mpi_integer,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_bcast ibuf error',ierr,rank
           stop
         endif
         if(n_var /= ibuf(1)) then
           print*,'error in main.f: n_var /= ibuf(1)',n_var,ibuf(1)
           stop
         endif
         max_obs    = ibuf(2)
         imax       = ibuf(3)
         jmax       = ibuf(4)
         kmax       = ibuf(5)
         i4time_sys = ibuf(6)
         ncnt_total = ibuf(7)
         idelt      = ibuf(8)
         jdelt      = ibuf(9)
         kdelt      = ibuf(10)
         deallocate  (ibuf)
         count      = 4*ncnt_total
         allocate(ibuf(count), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'error could not allocate ibuf',istat_alloc,rank,count
           stop
         endif
         call mpi_bcast(ibuf,count,mpi_integer,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_bcast ibuf 1 error',ierr,count,rank
           stop
         endif
         allocate(obs_barnes(max_obs), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' error: could not allocate abs_barnes',istat_alloc,
     1                         rank,max_obs
           stop
         endif
         do n=1,ncnt_total
           obs_barnes(n)%i      = ibuf(4*(n-1)+1)
           obs_barnes(n)%j      = ibuf(4*(n-1)+2)
           obs_barnes(n)%k      = ibuf(4*(n-1)+3)
           obs_barnes(n)%i4time = ibuf(4*(n-1)+4)
         enddo
         deallocate(ibuf)
!        receive an unpack real buffer
         count = (2+max_var)*ncnt_total+2
         allocate(rbuf(count), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'error could not allocate rbuf',istat_alloc,rank,count
           stop
         endif
         call mpi_bcast(rbuf,count,mpi_real,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_bcast rbuf error',ierr,count,rank
           stop
         endif
         r_missing_data_local = rbuf(1)
         r0_vert_grid         = rbuf(2)
         count                = 2
         do n=1,ncnt_total
           obs_barnes(n)%value(1:max_var)=rbuf(count+1:count+max_var)
           obs_barnes(n)%vert_rad_rat    =rbuf(count+1      +max_var)
           obs_barnes(n)%weight          =rbuf(count+2      +max_var)
           count                         =count+2+max_var
         enddo
         deallocate(rbuf)

!        receive arrays
         allocate(wt_b_local(imax,jmax,kmax,n_var), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' error: could not allocate wt_b_local in main',
     1              istat_alloc
           print*,rank,imax,jmax,kmax,n_var
           stop
         endif
         allocate(fnorm_lut(-(imax-1):(imax-1),-(jmax-1):(jmax-1)), 
     1                                           stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,' error: could not allocate fnorm_lut in main',
     1              istat_alloc,rank,imax,jmax
           stop
         endif
         call mpi_scatterv(buffer,recvcounts,displs,mpi_real,wt_b_local,
     1                  recvcounts(rank),mpi_real,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
          print*,'mpi_scatterv error in main',ierr,rank,recvcounts(rank)
          stop
         endif
         count = (2*(imax-1)+1)*(2*(jmax-1)+1)
         call mpi_bcast(fnorm_lut,count,mpi_real,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_bcast fnorm_lut error in main',ierr,count,rank
           stop
         endif
         allocate(local_var(istartiend(1,rank):istartiend(2,rank),jmax,
     1                                  kmax,n_var), stat=istat_alloc )
         if(istat_alloc .ne. 0)then
           print*,'error: could not allocate local_var in main',
     1            istat_alloc,rank
           print*,istartiend(1,rank),istartiend(2,rank),jmax,kmax,n_var
           stop
         endif
         call analyze_weights_mpi(
     1        n_var,max_obs,obs_barnes                                    ! i
     1       ,imax,jmax,kmax                                              ! i
     1       ,wt_b_local                                                  ! i
     1       ,i4time_sys                                                  ! i
     1       ,r_missing_data_local,r0_vert_grid                           ! i
     1       ,ncnt_total,idelt,jdelt,kdelt                                ! i
     1       ,fnorm_lut                                                   ! i
     1       ,istartiend(1,rank),istartiend(2,rank)                       ! i
     1       ,local_var)                                                  ! o
         call mpi_gatherv(local_var,recvcounts(rank),mpi_real,buffer,
     1               recvcounts,displs,mpi_real,0,mpi_comm_world,ierr)
         if( ierr /= 0) then
           print*,'mpi_gatherv error',ierr,rank,recvcounts(rank)
           stop
         endif

         deallocate(wt_b_local,fnorm_lut,obs_barnes,local_var)
       enddo !while
       call mpi_finalize(ierr)
       stop
      endif !rank > 0
      print*,'process',rank,' is proceeding'
#endif

!     obtain global parameters

!     read global parameters into module memory structure
      call get_directory('static',static_dir,len_dir)
      filename = static_dir(1:len_dir)//'/nest7grid.parms'
      call read_namelist_laps('lapsparms',filename)

!     read wind parameters into module memory structure
      filename = static_dir(1:len_dir)//'/wind.nl'
      call read_namelist_laps('wind',filename)

csms$create_decomp(grid_dh, <nx_l>, <0>)
csms$serial(<r_missing_data, out>:default=ignore)  begin              
#ifdef usempi
!     calculate the distribution of points in the i direction on each processor
      pointsperpe = nx_l/npes
      if(pointsperpe < 1) then
        print*,'error in main.f'
        print*,'rank,pointsperpe,nx_l,npes',rank,pointsperpe,nx_l,npes
        call flush(6)
        stop
      endif
      extra = nx_l - pointsperpe*npes
      if(extra > 0) then
        adjust=1
        extra = extra-1
      else
        adjust = 0
      endif
      istartiend(1,0) = 1
      istartiend(2,0) = pointsperpe+adjust
      do n=1,npes-1
        if(extra > 0) then
          adjust=1
          extra = extra-1
        else
          adjust = 0
        endif
        istartiend(1,n) = istartiend(2,n-1)+1
        istartiend(2,n) = istartiend(1,n)+pointsperpe-1+adjust
      enddo
      do n=0,npes-1
        print*,'distribution:',n,istartiend(1,n),istartiend(2,n),nx_l
      enddo
      call mpi_bcast(istartiend,2*npes,mpi_integer,0,mpi_comm_world,
     1                                                         ierr)
      if(ierr/=0) then
        print*,'bcast istartiend error',rank,2*npes,ierr
        stop
      endif
      do n=0,npes-1
        recvcounts(n)=(istartiend(2,n)-istartiend(1,n)+1)*ny_l*nk_laps
     1                                                           *n_var
      enddo
      displs(0) = 0
      do n=1,npes-1
        displs(n) = displs(n-1) + recvcounts(n-1)
      enddo
      call mpi_bcast(displs,npes,mpi_integer,0,mpi_comm_world,ierr)
      if( ierr /= 0) then
        print*,'mpi_bcast displs',ierr,npes,rank
        stop
      endif
      call mpi_bcast(recvcounts,npes,mpi_integer,0,mpi_comm_world,ierr)
      if( ierr /= 0) then
        print*,'mpi_bcast recvcounts',ierr,npes,rank
        stop
      endif
#endif

      ntmin = -1
      ntmax = +1

	call get_meso_sao_pirep(n_meso,n_sao,n_pirep,istatus)
	if (istatus .ne. 1) then
	   write (6,*) 'error obtaining n_meso, n_sao, and n_pirep'
	   stop
	endif

!     allocate static arrays (lat, lon, topo)
      allocate( lat(nx_l,ny_l), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate lat'
          stop
      endif

      allocate( lon(nx_l,ny_l), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate lon'
          stop
      endif

      allocate( topo(nx_l,ny_l), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate topo'
          stop
      endif

!     read static arrays (lat, lon, topo)
      call read_static_grid(nx_l,ny_l,'lat',lat,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error getting laps lat'
          stop
      endif

      call read_static_grid(nx_l,ny_l,'lon',lon,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error getting laps lon'
          stop
      endif

      call read_static_grid(nx_l,ny_l,'avg',topo,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error getting laps topo'
          stop
      endif

!     get system time
      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999
      write(6,*)' systime = ',a9_time

!     allocate global arrays
      allocate( heights_3d(nx_l,ny_l,nk_laps), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate heights_3d'
          stop
      endif

      allocate( wind%uanl(nx_l,ny_l,nk_laps), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate uanl'
          stop
      endif

      allocate( wind%vanl(nx_l,ny_l,nk_laps), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate vanl'
          stop
      endif

      allocate( wind%uanl_sfcitrp(nx_l,ny_l), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate uanl_sfcitrp'
          stop
      endif

      allocate( wind%vanl_sfcitrp(nx_l,ny_l), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate vanl_sfcitrp'
          stop
      endif

csms$serial end

      write(6,*)
      write(6,*)' getting 3_d height analysis from model:'

      istat_ht = 0
      var_2d = 'ht'
      call get_modelfg_3d(i4time_sys,var_2d
     1                             ,nx_l,ny_l,nk_laps
     1                             ,heights_3d,istat_ht)

      if(istat_ht .ne. 1)then
          write(6,*)
     1    ' aborting from laps wind anal - error in getting mdl heights'
          goto 999
      endif

      call lapswind_anal  (i4time_sys,
     1            ntmin,ntmax,                                     ! i 
     1            n_meso,n_sao,                                    ! i
     1            heights_3d,                                      ! i
     1            istat_wind)                                      ! o

      if(istat_wind .eq. 1 .and. iwrite_output .ge. 0)then
          ext = 'lw3'
          call write_wind_output(i4time_sys,ext,var_3d                
     1                            ,wind%uanl,wind%vanl                   ! i
     1                            ,wind%wanl                             ! i
     1                            ,wind%uanl_sfcitrp,wind%vanl_sfcitrp   ! i
     1                            ,nx_l,ny_l,nk_laps,n_3d_fields         ! i
     1                            ,r_missing_data                        ! i
     1                            ,istat_lw3)
      endif

!     if (allocated(wind%wanl) ) deallocate(wind%wanl)

999   continue
#ifdef usempi
      flag=-1
      call mpi_bcast(flag,1,mpi_integer,0,mpi_comm_world,ierr)
      if( ierr /= 0) then
        print*,'mpi_bcast flag final error',ierr,count,rank
        stop
      endif
      call mpi_finalize(ierr)
      stop
#endif

csms$exit
      end

