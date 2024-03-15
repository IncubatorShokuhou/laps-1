cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis

        subroutine get_mother_domain(ni,nj,lat,lon,istatus)
c
cdoc    reads in lat and lon fields for nested domains. only wrfsi
c
c       2003 john smart

        integer ni,nj                 ! input

        real   lat(ni,nj)             ! output
        real   lon(ni,nj)             ! output
        real   topo(ni,nj)            ! output
        real   rlaps_land_frac(ni,nj) ! output
        real   grid_spacing_m

        character*2   cpid
        character*3   var
        character*150 directory
        character*150 cstatic_name
        character*31  ext
        character*10  units
        character*125 comment

        include 'grid_fname.cmn'

        write(6,*)' reading in lat/lon for nest = ',nest

        ext = grid_fnam_common

!       get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        call s_len(grid_fnam_common,leng)
        call get_parent_id(iparent_id,istatus)
        write(cpid,'(i2.2)')iparent_id
        cstatic_name=grid_fnam_common(1:leng)//'.d'//cpid
        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        if(grid_fnam_common(1:leng).ne.'wrfsi')then
           print*,'error: get_domain_nest only works for wrfsi'
           stop
        endif
        var = 'lat'
        call rd_laps_static(directory,cstatic_name,ni,nj,1,var
     1,units,comment,lat,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading lac field'
            return
        endif

        var = 'lon'
        call rd_laps_static(directory,cstatic_name,ni,nj,1,var
     1,units,comment,lon,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading loc field'
            return
        endif

        return
        end

c------------------------------------------------------------------------
        subroutine get_laps_domain(ni,nj,grid_fnam,lat,lon,topo,istatus)
c
cdoc    reads in lat/lon/topo/landfrac fields
c
c       1994 steve albers

        integer ni,nj               ! input

        real lat(ni,nj)             ! output
        real lon(ni,nj)             ! output
        real topo(ni,nj)            ! output
        real rlaps_land_frac(ni,nj) ! local

        character*(*) grid_fnam       ! dummy

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)

        return
        end


        subroutine get_domain_laps(ni,nj,grid_fnam,lat,lon,topo
     1                             ,grid_spacing_cen_m,istatus)
c
cdoc    reads in lat/lon/topo fields with grid spacing
c
c       1994 steve albers
c
        integer ni,nj               ! input

        real lat(ni,nj)             ! output
        real lon(ni,nj)             ! output
        real topo(ni,nj)            ! output
        real rlaps_land_frac(ni,nj) ! local

        character*(*) grid_fnam       ! dummy

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)

        return
        end



        subroutine get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)
c
cdoc    reads in lat/lon/topo/landfrac fields with grid spacing
c
c       1994 steve albers

        integer ni,nj               ! input

        real lat(ni,nj)             ! output
        real lon(ni,nj)             ! output
        real topo(ni,nj)            ! output
        real rlaps_land_frac(ni,nj) ! output

        character*2   cnest
        character*3   var
        character*150 directory
        character*31  ext
        character*10  units
        character*125 comment

        character*6 c6_maproj

        include 'grid_fname.cmn'

        write(6,*)'    reading in lat/lon/topo/land frac '

        ext = grid_fnam_common

!       get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        call s_len(grid_fnam_common,leng)
        if(grid_fnam_common(1:leng).eq.'wrfsi')then
           write(cnest,'(i2.2)')nest
           ext=grid_fnam_common(1:leng)//'.d'//cnest
        else
           ext=grid_fnam_common
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        var = 'lat'
        if(grid_fnam_common(1:leng).eq.'wrfsi')then
           var = 'lac'
        endif
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                 ,lat,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading lat field'
            return
        endif

        var = 'lon'
        if(grid_fnam_common(1:leng).eq.'wrfsi')then
	   var = 'loc'
        endif
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,lon,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading lon field'
            return
        endif

        var = 'avg'
        if(grid_fnam_common(1:leng).eq.'wrfsi')then
	   var = 'avc'
        endif

        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,topo,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading avg (topo) field'
            return
        endif

c the staggered grid may have r_missing_data top and right row/column.
        call array_minmax(topo,ni,nj,rmin,rmax,r_missing_data)
        if(rmin .lt. -1000. .or. rmax .gt. 9000. .and.
     1     rmax.ne.r_missing_data)then
           write(6,*)' error, topo range out of bounds',rmin,rmax
           istatus = 0
           return
        endif

        var = 'ldf'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                 ,rlaps_land_frac,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading ldf (land fraction) field'
            return
        endif

c       write(6,*)' lat/lon corner > ',lat(   1,   1),lon(   1,   1)
c       write(6,*)' lat/lon corner > ',lat(   1,nj),lon(   1,nj)
c       write(6,*)' lat/lon corner > ',lat(ni,   1),lon(ni,   1)
c       write(6,*)' lat/lon corner > ',lat(ni,nj),lon(ni,nj)


        do i = 1,ni
        do j = 1,nj
            if(rlaps_land_frac(i,j) .le. 0.5)then
                rlaps_land_frac(i,j) = 0.           ! water
            else
                rlaps_land_frac(i,j) = 1.           ! land
            endif
        enddo
        enddo

        call check_domain(lat,lon,ni,nj,grid_spacing_m,5,istat_chk)
        if(istat_chk .ne. 1)then
            write(6,*)' warning or error in check_domain'
        endif

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading map projection'
            return
        endif

!       test for a conformal map projection
        if(c6_maproj .ne. 'latlon' .and. c6_maproj .ne. 'icshdr')then 
            icen = ni/2 + 1
            jcen = nj/2 + 1
            call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)       
     1                                  ,grid_spacing_cen_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error return from get_grid_spacing_actual'       
                return
            endif

        else
            write(6,*)' non-conformal map projection: ',c6_maproj
            write(6,*)' set grid_spacing_cen_m to parameter value'
            grid_spacing_cen_m = grid_spacing_m

        endif

        return
        end

        subroutine read_static_grid(ni,nj,var,static_grid,istatus)

cdoc    reads an arbitrary static field given an input 'var' string

c       2000    steve albers

        integer ni,nj                       ! input
        character*(*) var                     ! input
        real static_grid(ni,nj)             ! output

        character*150 directory
        character*31  ext
        character*10  units
        character*125 comment

        include 'grid_fname.cmn'

        write(6,*)'    reading in static grid: ',var

        ext = grid_fnam_common

!       get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                     ,static_grid,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading field ',var
            return
        endif

        return
        end

      subroutine force_get_laps_config(new_dataroot,istatus)

!     update the internal dataroot to supersede what was in $laps_data_root

      include 'lapsparms.cmn'
      include 'grid_fname.cmn'
      character*(*) new_dataroot
      integer istatus
      iflag_lapsparms=0
      generic_data_root = new_dataroot
      call get_laps_config_sub('nest7grid',istatus)
      return
      end

      subroutine get_laps_config(grid_fnam,istatus)
c
c higher level wrapper to allow get_laps_config to use new wrf si
c namelist with hardwire grid_fnam (as is the case in most of the repository).
c
      implicit none

      integer istatus
      character*(*) grid_fnam

      include 'grid_fname.cmn'

      grid_fnam_common = grid_fnam

      call get_config(istatus)
      if(istatus .ne. 1)then
         print*,'error returned from get_config'
         return
      endif
      
      return
      end
c --------------------------------------------------------------------

      subroutine get_config(istatus)

        implicit none

        integer max_files
        parameter (max_files = 500)

        character*80  domain_name
        integer       len_grid_fnam

        include   'grid_fname.cmn'

        character c_filenames(max_files)*255
        character cfilespec*255
        character directory*200
        character generic_data_root_laps*200
        character generic_data_root_wrfsi*200
        character exe*200
        integer   i,j,jj
        integer   lend,lens,len_root
        integer   len_root_laps
        integer   len_root_wrfsi
        integer   numoffiles
        integer   narg,iargc
        integer   istatus

        integer   idir_len
        data      idir_len/0/
        save      idir_len
c
c laps_data_root is either environment variable or command line input.
c command line overrides env variable if it exists.
c
c lets be sure only one "data_root" is on. if both, then we should
c terminate since this is ambiguous.
c
        istatus = 0

        call getenv('laps_data_root',generic_data_root_laps)
        call s_len(generic_data_root_laps,len_root_laps)
        call getenv('moad_dataroot',generic_data_root_wrfsi)
        call s_len(generic_data_root_wrfsi,len_root_wrfsi)

        if(len_root_laps.gt.0 .and. len_root_wrfsi.gt.0)then
           print*
           print*,'---------------------------------------------------'
           print*,'*** error: two dataroots set => ambiguous info '
           print*,'*** laps_data_root and moad_dataroot are both set.'
           print*,'*** unable to proceed. please fix your enviroment'
           print*,'---------------------------------------------------'
           print*

           print*,'len_root_laps = ',len_root_laps
           if(len_root_laps.gt.0)print*,'laps generic dataroot',
     +generic_data_root_laps(1:len_root_laps)
           print*,'len_root_wrfsi = ',len_root_wrfsi
           if(len_root_wrfsi.gt.0)print*,'moad generic dataroot',
     +generic_data_root_wrfsi(1:len_root_wrfsi)

           return
        endif

        if(idir_len.eq.0) then
           narg = iargc()
           if(narg.gt.0)then
              call getarg(0,exe)
              call getarg(narg,generic_data_root)
              call s_len(generic_data_root,len_root)
              if(generic_data_root(len_root-1:len_root).eq.'.x'.or.
     1           generic_data_root(len_root-3:len_root).eq.'.exe' .or.
     1           index(exe,'laps2grib.exe') .ge. 1)
     1        then 
                 print*,'not a typical dataroot on command line'
                 print*,'trying laps_data_root environment variable'
                 call getenv('laps_data_root',generic_data_root)
                 call get_directory_length(generic_data_root
     1                 ,len_root)
                 if(len_root.eq.0)then
                    call getenv('moad_dataroot',generic_data_root)
                    call get_directory_length(generic_data_root
     1,len_root)
                 endif
              endif
           else
              call getenv('laps_data_root',generic_data_root)
              call get_directory_length(generic_data_root,len_root)
              if(len_root.eq.0)then
                 call getenv('moad_dataroot',generic_data_root)
                 call get_directory_length(generic_data_root,len_root)
              endif
           endif
           if(len_root.eq.0)then
              call s_len(generic_data_root,len_root)
              if(len_root.eq.0)then
                 print*,'use either command line or env variable for ',
     +'system dataroot'
                 print*,'stop in get_config'
                 stop
              endif
           endif

           call find_domain_name(generic_data_root,domain_name,istatus)
           if(istatus.ne.1)then
              print*,'error returned from find_domain_name'
              return
           endif
           call s_len(grid_fnam_common,len_grid_fnam)
           idir_len = len_grid_fnam

        endif


c use the "grid namelist" to load lapsparms.cmn with appropriate values.
        if(grid_fnam_common(1:9).eq.'nest7grid')then

c****  because laps does not "nest", we hardwire it as 1 for now *****
c****  --------------------------------------------------------- *****
           nest=1
           call get_laps_config_sub(grid_fnam_common,istatus)
           if(istatus.ne.1)then
              print*,'error returned from get_laps_config_sub'
              return
           endif
        elseif(grid_fnam_common(1:5).eq.'wrfsi')then
!          commented out by sa as a quick fix with respect to parameter name
!          changes.
!          call get_wrfsi_config(nest,istatus)
           if(istatus.ne.1)then
              print*,'error returned from get_wrfsi_config'
              return
           endif
        else
           print*,'unknown grid filename spec'
           return
        endif

        istatus = 1

        return
        end

        subroutine get_laps_config_sub(grid_fnam,istatus)

!       1992 steve albers

!       read in parameters from parameter file. the first time the routine
!       is called, parameters are read in from the namelist. subsequent
!       calls will return these parameters (now stored in common), without
!       rereading the namelist. this is more efficient if the routine is called
!       many times.

        include 'lapsparms.cmn'

        character*(*) grid_fnam   ! input (warning: trailing blanks won't work)
        character*150  directory
        character*150 static_dir,filename
        character*31  ext
        character*8 a8
        character*200 tempchar

        integer init
        data init/0/
        save init

        include 'grid_fname.cmn'

        if(init.eq.1 .and. iflag_lapsparms .eq. 1) then ! data already read in
!          print *, 'it works'
!          goto999                                           ! normal return
!          init = 1
           istatus = 1
           return
        else
           iflag_lapsparms = 0 ! initializes variable
        endif

!       while we are here, let's put the grid name into the common area
        grid_fnam_common = grid_fnam  ! used in get_directory to modify
                                      ! extension based on the grid domain

!       get the location of the parameter directory
        call s_len(grid_fnam,len_grid_fnam)
c        len_dir = index(grid_fnam_common,'/',.true.)
        call get_directory_length(grid_fnam_common,len_dir)

        if(len_dir.gt.0) then
           ext='nest7grid'
!          directory = generic_data_root(1:len_root)//'static/'
!          call s_len(directory,len_dir)
        else
           ext=grid_fnam
!          directory = generic_data_root(1:len_root)//'static/'
!          call s_len(directory,len_dir)
        endif

        fdda_model_source = '         '

!       this recursive call may not be needed if we assign the directory
!       just above.
        call get_directory(ext,directory,len_dir)

! this is laps specific for nest7grid.parms ... the laps namelist file.
        call s_len(ext,len_ext)
        if(directory(len_dir:len_dir).ne.'/') then
          tempchar = directory(1:len_dir)//'/'//ext(1:len_ext)//'.parms'
        else
          tempchar = directory(1:len_dir)//ext(1:len_ext)//'.parms'
        endif
        call s_len(tempchar,len_dir)
 

        min_to_wait_for_metars=10

!       read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

        if(iflag_lapsparms .ne. 1)then                      ! error return
            goto910                                    

        else                                                    ! normal return
            write(6,*)tempchar(1:len_dir),
     +            ' get_laps_config - parameters read in ok'
            goto999                                     

        endif

910     write(6,*)' read error in get_laps_config'              ! error return
        write(6,*)' check runtime parameter file ',tempchar
        iflag_lapsparms = 0
        istatus = 0
        return

920     write(6,*)' read error in get_laps_config'              ! error return
        write(6,*)' truncated runtime parameter file ',tempchar
        iflag_lapsparms = 0
        istatus = 0
        return

999     continue                                                ! normal return

        init = 1
        istatus = 1
        return

        end
c ----------------------------------------------------------------------------
      subroutine find_domain_name(data_root, domain_name,
     +                            istatus)

        implicit none

        integer max_files
        parameter (max_files = 500)

        character*(*)  domain_name
        character*(*)  data_root
        integer   istatus

        character grid_fnam*80
        integer   len_grid_fnam, len_gfc
        integer   lend,lens,len_root
        character c_filenames(max_files)*255
        character cfilespec*255
        integer   numoffiles, idir_len, lgfc
        integer   i, j, jj
        logical   lfound_name

        include 'grid_fname.cmn'

        istatus = 0

        call s_len(grid_fnam_common,lgfc)
        if(lgfc.gt.0)then
           domain_name = grid_fnam_common
           istatus = 1
           return
        endif

c  determine the grid_fnam by searching the directory for
c  known namelists files used by this system
c (currently nest7grid.parms and wrfsi.nl)
c  assume the namelist file is always in data_root/static.

        lfound_name = .false.
        i=1
        call s_len(data_root,len_root)
        cfilespec=data_root(1:len_root)//'/static/*'
        call get_file_names(cfilespec,numoffiles,c_filenames,
     +max_files,istatus)
        do while (i.le.numoffiles.and.(.not.lfound_name))
          call get_directory_length(c_filenames(i),lend)
          call s_len(c_filenames(i),lens)
          if(c_filenames(i)(lend+1:lens).eq.'nest7grid.parms'
     +  .or. c_filenames(i)(lend+1:lens).eq.'wrfsi.nl')then
              jj = 0
              do j=lend,lens
                 jj=jj+1
                 if(c_filenames(i)(j:j).eq.'.')then
                    grid_fnam = c_filenames(i)(lend+1:lend+jj-2)
                 endif
              enddo
              lfound_name = .true.
           endif
           i=i+1
        enddo

        if(.not. lfound_name .or. numoffiles .eq. 0)then
          print*, 'error in getting domain name from ',
     +           data_root(1:len_root),'/static/'
          print*, 'lfound_name = ',lfound_name
          print*, 'numoffiles = ',numoffiles
          istatus = 0
          return
        endif

        call s_len(grid_fnam,len_grid_fnam)
        if (len_grid_fnam .gt. 0) then

          if ((grid_fnam(1:len_grid_fnam) .eq. 'nest7grid') .or.
     +        (grid_fnam(1:len_grid_fnam) .eq. 'wrfsi')) then
            istatus = 1
c     set the grid fname common block value
            len_gfc = len(grid_fnam_common)
            if (len_grid_fnam .le. len_gfc) then
              grid_fnam_common = grid_fnam
              domain_name = grid_fnam
            else
              print*, 'error length of domain name exceeds size of'
     +               ,' grid_fnam_common: '
              print*, 'name= ', grid_fnam
              istatus = 0
            endif
          else
            print*, 'error in getting domain name from ',
     +             data_root(1:len_root),'/static/'
            istatus = 0
          endif
        else
          print*, 'error in getting domain name from ',
     +           data_root(1:len_root),'/static/'
          istatus = 0
        endif

        return
        end
c ---------------------------------------------------------------------
      subroutine get_mother_dims(nx_mother,ny_mother,istatus)
      include 'lapsparms.cmn' ! standard_longitude
      include 'grid_fname.cmn'! grid_fnam_common
      include 'wrf_horzgrid.cmn' ! entire hgridspec namelist section

!     this routine accesses the standard_longitude variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then 
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      if(parent_id(nest).ne.1)then
         nx_mother=(domain_origin_uri(parent_id(nest))-
     +           domain_origin_lli(parent_id(nest)))*
     +           ratio_to_parent(parent_id(nest))+1
         ny_mother=(domain_origin_urj(parent_id(nest))-
     +           domain_origin_llj(parent_id(nest)))*
     +           ratio_to_parent(parent_id(nest))+1
      else
         nx_mother=xdim
         ny_mother=ydim
      endif 

      return
      end

c ---------------------------------------------------------------------
      subroutine get_standard_longitude(std_lon_ret,istatus)

      include 'lapsparms.cmn' ! standard_longitude
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the standard_longitude variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      std_lon_ret = standard_longitude

      istatus = 1
      return
      end


      subroutine get_grid_spacing(grid_spacing_m_ret,istatus)

      include 'lapsparms.cmn' ! grid_spacing_m
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      grid_spacing_m_ret = grid_spacing_m

      istatus = 1
      return
      end


      subroutine get_grid_center(grid_cen_lat_ret,grid_cen_lon_ret
     1                          ,istatus)

      include 'lapsparms.cmn' ! grid_cen_lat,grid_cen_lon
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      grid_cen_lat_ret = grid_cen_lat
      grid_cen_lon_ret = grid_cen_lon

      istatus = 1
      return
      end

      subroutine get_standard_latitude(std_lat,istatus)

      include 'lapsparms.cmn' ! standard_latitude
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the standard_latitude variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      std_lat = standard_latitude

      istatus = 1
      return
      end


      subroutine get_standard_latitudes(std_lat1,std_lat2,istatus)

      include 'lapsparms.cmn' ! standard_latitude, standard_latitude2
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the standard_latitude variables from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      std_lat1 = standard_latitude
      std_lat2 = standard_latitude2

      istatus = 1
      return
      end

      subroutine get_vertical_grid(vert_grid,istatus)

      include 'lapsparms.cmn' ! standard_latitude, standard_latitude2
      include 'grid_fname.cmn'! grid_fnam_common

      character*40 vert_grid

!     this routine accesses the vertical_grid variables from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      vert_grid = vertical_grid

      istatus = 1
      return
      end


      subroutine get_maxstns(maxstns_ret,istatus)

      include 'lapsparms.cmn' ! maxstns
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the maxstns variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      maxstns_ret = maxstns

      istatus = 1
      return
      end


      subroutine get_c8_project(c8_project_ret,istatus)

      include 'lapsparms.cmn' ! c8_project
      include 'grid_fname.cmn'! grid_fnam_common

      character*8 c8_project_ret

!     this routine accesses the c8_project_common variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' warning, get_laps_config not successfully called'       
          return
      endif

      c8_project_ret = c8_project

      istatus = 1
      return
      end

      subroutine get_c8_blpfmt(c8_blpfmt_ret,istatus)

      include 'lapsparms.cmn' ! c8_blpfmt
      include 'grid_fname.cmn'! grid_fnam_common

      character*8 c8_blpfmt_ret

!     this routine accesses the c8_blpfmt_common variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' warning, get_laps_config not successfully called'       
          return
      endif

      c8_blpfmt_ret = c8_blpfmt

      istatus = 1
      return
      end

c------------------------------------------------------------
      subroutine get_raddat_type(c_raddat_type_ret,istatus)
      include 'lapsparms.cmn' ! c_raddat_type
      include 'grid_fname.cmn'! grid_fnam_common

      character*3 c_raddat_type_ret

!     this routine accesses the c8_project_common variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' warning, get_laps_config not successfully called'
          return
      endif

      c_raddat_type_ret = c_raddat_type

      istatus = 1
      return
      end
c----------------------------------------------------------

      subroutine get_c6_maproj(c6_maproj_ret,istatus)

      include 'lapsparms.cmn' ! c6_maproj
      include 'grid_fname.cmn'! grid_fnam_common

      character*6 c6_maproj_ret

!     this routine accesses the c6_maproj variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      c6_maproj_ret = c6_maproj

      istatus = 1
      return
      end


      subroutine get_c80_description(c80_description_ret,istatus)

      include 'lapsparms.cmn' ! c80_description
      include 'grid_fname.cmn'! grid_fnam_common

      character*80 c80_description_ret
!     this routine accesses the c80_description variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      c80_description_ret = c80_description

      istatus = 1
      return
      end


      subroutine get_laps_dimensions(nk,istatus)

      include 'lapsparms.cmn'              ! nk_laps
      include 'grid_fname.cmn'             ! grid_fnam_common

      integer init
      data init/0/
      save init

!     this routine accesses the nk variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms .ne. 1 .or. init .eq. 0)then 
          call get_laps_config(grid_fnam_common,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1            ' error, get_laps_config not successfully called'      
              return
          endif
          init = 1
      endif

      if(nk_laps .gt. max_lvls)then
          write(6,*)
     1  ' error: number of laps levels exceeds coded limit of ',max_lvls
          write(6,*)' to increase the limit edit lapsparms.for'
          istatus = 0
          return
      endif

      nk = nk_laps

      istatus = 1
      return
      end


      subroutine get_aircraft_time_window(aircraft_time_window_ret
     1                                   ,istatus)       

      include 'lapsparms.cmn'              ! laps_cycle_time
      include 'grid_fname.cmn'             ! grid_fnam_common

!     this routine accesses the aircraft_time_window variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      aircraft_time_window_ret = aircraft_time_window

      istatus = 1
      return
      end

      subroutine get_laps_cycle_time(laps_cycle_time_ret,istatus)

      include 'lapsparms.cmn'              ! laps_cycle_time
      include 'grid_fname.cmn'             ! grid_fnam_common

!     this routine accesses the laps_cycle_time variable from the
!     .parms file via the common block. note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      laps_cycle_time_ret = laps_cycle_time

      istatus = 1
      return
      end

      subroutine get_num_domains(num_domains,istatus)

      include 'lapsparms.cmn'              ! num_domains_cmn
      include 'grid_fname.cmn'             ! grid_fnam_common

!     character*80 grid_fnam

!     this routine accesses the num_domains variable from the
!     namelist file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      num_domains = num_domains_cmn

      istatus = 1
      return
      end
c---------------------------------------------------------------
      subroutine get_domain_origin(lli_orig,llj_orig
     +,uri_orig,urj_orig,istatus)

      include 'lapsparms.cmn'              ! lli_orig_cmn,llj_orig_cmn,uri_orig_cmn,urj_orig_cmn
      include 'grid_fname.cmn'             ! grid_fnam_common

      integer  lli_orig,llj_orig
      integer  uri_orig,urj_orig

!     character*80 grid_fnam

!     this routine accesses the ll/ur i&j_domain_origin variables from the
!     wrfsi.nl namelist file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      if(lli_orig_cmn.eq.0)lli_orig_cmn=1
      if(llj_orig_cmn.eq.0)llj_orig_cmn=1
      if(uri_orig_cmn.eq.0 .or. urj_orig_cmn.eq.0)then
         call get_grid_dim_xy(nx,ny,istatus)
         uri_orig_cmn=nx
         urj_orig_cmn=ny
      endif

      lli_orig=lli_orig_cmn
      llj_orig=llj_orig_cmn
      uri_orig=uri_orig_cmn
      urj_orig=urj_orig_cmn

      istatus = 1
      return
      end
c -------------------------------------------------------------
c
      subroutine get_ratio_to_parent(iratio,istatus)

      include 'lapsparms.cmn'              ! ratio_to_parent_cmn
      include 'grid_fname.cmn'             ! grid_fnam_common

!     character*80 grid_fnam

!     this routine accesses the ll/ur i&j_domain_origin_parent variables from the
!     wrfsi.nl namelist file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      if(ratio_to_parent_cmn.eq.0)ratio_to_parent_cmn=1
      iratio=ratio_to_parent_cmn

      istatus = 1
      return
      end
c---------------------------------------------------------------
c
      subroutine get_parent_id(parent_id,istatus)

      include 'lapsparms.cmn'              ! parent_id_cmn
      include 'grid_fname.cmn'             ! grid_fnam_common

      integer  parent_id

!     character*80 grid_fnam

!     this routine accesses the parent_id variable from the
!     wrfsi.nl namelist file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      if(parent_id_cmn.eq.0)parent_id_cmn=1
      parent_id=parent_id_cmn

      istatus = 1
      return
      end

      subroutine get_grid_dim_xy(nx_l_ret,ny_l_ret,istatus)

      include 'lapsparms.cmn'              ! nx_l, ny_l
      include 'grid_fname.cmn'             ! grid_fnam_common

!     character*80 grid_fnam

!     this routine accesses the nx_l and ny_l variables from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

!     grid_fnam = grid_fnam_common
      print*,'here: get_grid_dim_xy'
      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      nx_l_ret = nx_l
      ny_l_ret = ny_l

      istatus = 1
      return
      end

      subroutine get_topo_parms(silavwt_parm_ret,toptwvl_parm_ret
     1                         ,istatus)

      include 'lapsparms.cmn' ! silavwt, toptwvl
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the silavwt_parm and toptwvl_parm
!     variables from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      silavwt_parm_ret = silavwt_parm
      toptwvl_parm_ret = toptwvl_parm

      istatus = 1
      return
      end

      subroutine get_meso_sao_pirep(n_meso,n_sao,n_pirep_ret,istatus)

      include 'lapsparms.cmn' ! maxstns, n_pirep
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the maxstns and n_pirep
!     variables from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      n_meso  = maxstns
      n_sao   = maxstns
      n_pirep_ret = n_pirep

      istatus = 1
      return
      end

      subroutine get_max_radars (max_radars_ret, istatus)

      include 'lapsparms.cmn' ! max_radars
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the max_radars variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      max_radars_ret = max_radars

      istatus = 1
      return
      end

      subroutine get_max_stations (maxstns_ret, istatus)

      include 'lapsparms.cmn' ! maxstns
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the maxstns variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      maxstns_ret = maxstns

      istatus = 1
      return
      end

      subroutine get_vert_rads (vert_rad_pirep_ret,
     1                          vert_rad_sao_ret,
     1                          vert_rad_meso_ret,
     1                          vert_rad_prof_ret,
     1                          istatus)

      include 'lapsparms.cmn' ! vert_rad_pirep, etc.
      include 'grid_fname.cmn'! grid_fnam_common

      integer vert_rad_pirep_ret
      integer vert_rad_sao_ret
      integer vert_rad_meso_ret
      integer vert_rad_prof_ret

!     this routine accesses the vert_rad_pirep, etc., variables from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      vert_rad_pirep_ret = vert_rad_pirep
      vert_rad_sao_ret = vert_rad_sao
      vert_rad_meso_ret = vert_rad_meso
      vert_rad_prof_ret = vert_rad_prof

      istatus = 1
      return
      end

      subroutine get_r_missing_data(r_missing_data_ret, istatus)

      include 'lapsparms.cmn' ! r_missing_data
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the r_missing_data variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      r_missing_data_ret = r_missing_data

      istatus = 1
      return
      end

      subroutine get_i2_missing_data(i2_missing_data_ret, istatus)

      include 'lapsparms.cmn' ! i2_missing_data
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the i2_missing_data variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      i2_missing_data_ret = i2_missing_data

      istatus = 1
      return
      end


      subroutine get_i_perimeter(i_perimeter_ret, istatus)

      include 'lapsparms.cmn' ! i_perimeter
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the 'i_perimeter' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      i_perimeter_ret = i_perimeter

      istatus = 1
      return
      end


      subroutine get_ref_base(ref_base_ret, istatus)

      include 'lapsparms.cmn' ! ref_base
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the 'ref_base' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      ref_base_ret = ref_base

      istatus = 1
      return
      end

      subroutine get_ref_base_useable(ref_base_useable_ret, istatus)

      include 'lapsparms.cmn' ! ref_base_useable
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the 'ref_base_useable' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      ref_base_useable_ret = ref_base_useable

      istatus = 1
      return
      end

      subroutine get_r_hybrid_first_gate(r_hybrid_first_gate_ret
     1                                 , istatus)       

      include 'lapsparms.cmn' ! r_hybrid_first_gate
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the 'r_hybrid_first_gate' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      r_hybrid_first_gate_ret = r_hybrid_first_gate

      istatus = 1
      return
      end

      subroutine get_max_radar_files(max_radar_files_ret, istatus)       

      include 'lapsparms.cmn' ! max_radar_files
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the 'max_radar_files' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'       
          return
      endif

      max_radar_files_ret = max_radar_files_nl

      istatus = 1
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_l_compress_radar(l_compress_radar_ret, istatus)

      include 'lapsparms.cmn' ! l_compress_radar
      include 'grid_fname.cmn'! grid_fnam_common

      logical l_compress_radar_ret

!     this routine accesses the 'l_compress_radar' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error, get_laps_config not successfully called'       
         return
      endif

      l_compress_radar_ret = l_compress_radar

      istatus = 1
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_l_use_tamdar(l_use_tamdar_ret, istatus)

      include 'lapsparms.cmn' ! l_use_tamdar
      include 'grid_fname.cmn'! grid_fnam_common

      logical l_use_tamdar_ret

!     this routine accesses the 'l_use_tamdar' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error, get_laps_config not successfully called'       
         return
      endif

      l_use_tamdar_ret = l_use_tamdar

      istatus = 1
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_l_3dvar(l_3dvar_ret, istatus)

      include 'lapsparms.cmn' ! l_3dvar
      include 'grid_fname.cmn'! grid_fnam_common

      logical l_3dvar_ret

!     this routine accesses the 'l_3dvar' variable from the
!     .parms file via the common block. note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error, get_laps_config not successfully called'       
         return
      endif

      l_3dvar_ret = l_3dvar

      istatus = 1
      return
      end
c
c-----------------------------------------------------------------------
c
      function ltest_vertical_grid(c_vertical_grid)

!     the input c_vertical_grid is a string you are testing
!     the output (.true. or .false.) states whether the input vertical_grid 
!     name equals the vertical_grid in the namelist

      include 'lapsparms.cmn' ! vertical_grid
      include 'grid_fname.cmn'! grid_fnam_common

      logical ltest_vertical_grid
      character*(*) c_vertical_grid
      character*40  cvgrid_test
      character*40  cvgrid_laps

      integer init
      data init/0/
      save init

      save cvgrid_laps,len_cmn

      if(iflag_lapsparms .ne. 1 .or. init .eq. 0)then 
          call get_laps_config(grid_fnam_common,istatus)

          if(istatus .ne. 1)then
              write(6,*)' ltest_vertical_grid: error detected in '
     1                 ,'calling get_laps_config'
              ltest_vertical_grid = .false.
              return
          endif

          call downcase(vertical_grid,cvgrid_laps)
          call s_len(cvgrid_laps,len_cmn)
          init = 1
      endif

      call downcase(c_vertical_grid,cvgrid_test)
      call s_len(cvgrid_test,len_in)

      if(cvgrid_laps(1:len_cmn).eq.cvgrid_test(1:len_in))then
          ltest_vertical_grid = .true.
      else
          ltest_vertical_grid = .false.
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      function ltest_vertical_grid_lc(c_vertical_grid)

!     the input c_vertical_grid is a string you are testing
!     the output (.true. or .false.) states whether the input vertical_grid 
!     name equals the vertical_grid in the namelist

!     the input vertical grid must be in lower case

      include 'lapsparms.cmn' ! vertical_grid
      include 'grid_fname.cmn'! grid_fnam_common

      logical ltest_vertical_grid_lc
      character*(*) c_vertical_grid
      character*40  cvgrid_test
      character*40  cvgrid_laps

      integer init
      data init/0/
      save init

      save cvgrid_laps,len_cmn

      if(iflag_lapsparms .ne. 1 .or. init .eq. 0)then 
          call get_laps_config(grid_fnam_common,istatus)

          if(istatus .ne. 1)then
              write(6,*)' ltest_vertical_grid_lc: error detected in '
     1                 ,'calling get_laps_config'
              ltest_vertical_grid_lc = .false.
              return
          endif

          call downcase(vertical_grid,cvgrid_laps)
          call s_len(cvgrid_laps,len_cmn)
          init = 1
      endif

      call s_len(c_vertical_grid,len_in)

      if(cvgrid_laps(1:len_cmn).eq.c_vertical_grid(1:len_in))then
          ltest_vertical_grid_lc = .true.
      else
          ltest_vertical_grid_lc = .false.
      endif

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_earth_radius(erad,istatus)

      use mem_namelist, only : earth_radius

      include 'grid_fname.cmn'             ! grid_fnam_common

      integer init
      data init/0/
      save init

      if(iflag_lapsparms .ne. 1 .or. init .eq. 0)then 
          call get_laps_config(grid_fnam_common,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1            ' error, get_laps_config not successfully called'      
              return
          endif
          init = 1
      endif

      erad = earth_radius

      istatus = 1

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_fdda_model_source(fdda_model_source_ret
     1,n_fdda_models,istatus)

      include 'lapsparms.cmn' ! fdda_model_source
      include 'grid_fname.cmn'! grid_fnam_common

      character*(*) fdda_model_source_ret(maxbgmodels)
      integer n_fdda_models,istatus

!     this routine accesses the fdda_model_source variable from the
!     .parms file via the common block. note the variable name in the
!     argument list is different in the calling routine

c     do i=1,maxbgmodels
c        fdda_model_source(i) = ' '
c     enddo
 
      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      n_fdda_models = 0
      do i=1,maxbgmodels
         if(fdda_model_source(i).ne. ' ')then
            n_fdda_models=n_fdda_models+1
            fdda_model_source_ret(n_fdda_models)=fdda_model_source(i)     
         endif
      enddo

      return
      end

      subroutine get_path_to_topo_10m(path_to_topo_10m_ret,istatus)

      include 'lapsparms.cmn' ! path_to_topt_10m
      include 'grid_fname.cmn'

!     this routine accesses the path_to_topt_10m  variable from the
!     .parms file via the common block.

      character*200 path_to_topo_10m_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_topo_10m_ret =  path_to_topt10m

      return
      end

      subroutine get_path_to_topo_30s(path_to_topo_30s_ret,istatus)

      include 'lapsparms.cmn' ! path_to_topt_30s
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_topt_10m  variable from the
!     .parms file via the common block. 

      character*200 path_to_topo_30s_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_topo_30s_ret =  path_to_topt30s

      return
      end

      subroutine get_path_to_soiltype_top(path_to_soiltype_top_30s_ret
     +,istatus)

      include 'lapsparms.cmn' ! path_to_soiltype_top30s
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_soiltype_top30s  variable from the
!     .parms file via the common block.

      character*200 path_to_soiltype_top_30s_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_soiltype_top_30s_ret =  path_to_soiltype_top30s

      return
      end

      subroutine get_path_to_soiltype_bot(path_to_soiltype_bot_30s_ret
     +,istatus)

      include 'lapsparms.cmn' ! path_to_soiltype_bot30s
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_soiltype_bot30s  variable from the
!     .parms file via the common block.

      character*200 path_to_soiltype_bot_30s_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_soiltype_bot_30s_ret =  path_to_soiltype_bot30s

      return
      end

      subroutine get_path_to_landuse_30s(path_to_landuse_30s_ret
     1                                  ,istatus)

      include 'lapsparms.cmn' ! path_to_landuse30s
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_landuse30s  variable from the
!     .parms file via the common block.

      character*200 path_to_landuse_30s_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_landuse_30s_ret =  path_to_landuse30s

      return
      end
c
c-------------------------------------------------------------
      subroutine get_path_to_green_frac(path_to_green_frac_ret
     &,istatus)

      include 'lapsparms.cmn' ! path_to_green_frac
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_green_frac variable from the
!     .parms file via the common block.

      character*200 path_to_green_frac_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_green_frac_ret =  path_to_greenfrac

      return
      end
c
c---------------------------------------------------------------
      subroutine get_path_to_soiltemp_1deg(path_to_soiltemp_1deg_ret
     &,istatus)

      include 'lapsparms.cmn' ! path_to_soiltemp1deg
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_landuse30s  variable from the
!     .parms file via the common block.

      character*200 path_to_soiltemp_1deg_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_soiltemp_1deg_ret =  path_to_soiltemp1deg

      return
      end
c
c---------------------------------------------------------------
      subroutine get_path_to_albedo(path_to_albedo_ret,istatus)

      include 'lapsparms.cmn' ! path_to_albedo
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_albedo  variable from the
!     .parms file via the common block.

      character*200 path_to_albedo_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_albedo_ret =  path_to_albedo

      return
      end
c
c---------------------------------------------------------------
      subroutine get_path_to_maxsnoalb(pathtomaxsnoalb_ret,istatus)

      include 'lapsparms.cmn' ! path_to_maxsnoalb
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_maxsnoalb variable from the
!     .parms file via the common block.

      character*200 pathtomaxsnoalb_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      pathtomaxsnoalb_ret = path_to_maxsnoalb

      return
      end
c
c---------------------------------------------------------------
	subroutine get_path_to_islope(pathtoislope_ret,istatus)

	include 'lapsparms.cmn' ! path_to_islope
	include 'grid_fname.cmn'! grid_fnam_common

!	this routine accesses the path_to_islope variable from the
!	.parms file via the common block.

	character*200 pathtoislope_ret

	call get_laps_config(grid_fnam_common,istatus)

	if(istatus .ne. 1)then
	  write(6,*)' error, get_laps_config not successfully called'
	  return
	endif

        pathtoislope_ret = path_to_islope

	return
	end
c
c---------------------------------------------------------------
      subroutine get_path_to_sst(path_to_sst_ret,istatus)

      include 'lapsparms.cmn' ! path_to_sst
      include 'grid_fname.cmn'! grid_fnam_common

!     this routine accesses the path_to_sst  variable from the
!     .parms file via the common block.

      character*200 path_to_sst_ret

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      path_to_sst_ret =  path_to_sst

      return
      end
c--------------------------------------------------------------
      subroutine get_n_staggers(n_staggers,istatus)
!
! currently we do not carry an "n_staggers" variable in
! nest7grid.parms
!
      include 'lapsparms.cmn' ! path_to_sst
      include 'grid_fname.cmn'! grid_fnam_common

      integer istatus
      integer n_staggers
      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      call s_len(grid_fnam_common,nf)
      if(grid_fnam_common(1:nf).eq.'nest7grid')then
         n_staggers=1
      else
         n_staggers=4
      endif
      return
      end
c--------------------------------------------------------------
      subroutine get_stagger_index(istag,istatus)
!
! currently we do not carry an "n_staggers" variable in
! nest7grid.parms
!
      include 'lapsparms.cmn'           ! for stagger_type
      include 'grid_fname.cmn'          ! grid_fnam_common

      integer istatus
      integer istag
      character*4 stagger_type          !variable soon to be in .cmn

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' error, get_laps_config not successfully called'
          return
      endif

      call s_len(grid_fnam_common,nf)
      if(grid_fnam_common(1:nf).eq.'nest7grid')then
         istag=1
      else
         stagger_type='a-c'             !hardwire testing for now
!        stagger_type='a'
         call s_len(stagger_type,ls)
         if(ls.eq.0 .or. ls.gt.4)then
            print*,'error detected in stagger_type variable'
            if(ls.ne.0)print*,'stagger type = ',stagger_type(1:ls)
            if(ls.eq.0)print*,'stagger type variable has 0 length'
            istatus = 0
            return
         else
            call downcase(stagger_type,stagger_type)
            if(stagger_type(1:ls).eq.'a-c')then
               istag=4
            elseif(stagger_type(1:ls).eq.'a')then
               istag=1
            else
               print*,' cannot set istag in get_stagger_index'
               print*,' stagger_type = ',stagger_type
               istatus = 0
            endif
         endif
      endif

      return
      end
c -----------------------------------------------------------------
      subroutine array_minmax(a,ni,nj,rmin,rmax,r_missing_data)

      real a(ni,nj)

      rmin =  abs(r_missing_data)
      rmax = -abs(r_missing_data)

      do i = 1,ni
      do j = 1,nj
          rmin = min(rmin,a(i,j))
          rmax = max(rmax,a(i,j))
      enddo ! j
      enddo ! i

      return
      end
