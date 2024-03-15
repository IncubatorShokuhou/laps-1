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
c
c

        subroutine get_sat_data(i4time,
     1  i4_sat_window,i4_sat_window_offset,                              ! i
     1  imax,jmax,r_missing_data,                                        ! i
     1  l_use_39,l_use_co2,latency_co2,                                  ! i
     1  lat,lon,                                                         ! i
     1  subpoint_lat_clo_s8a,subpoint_lon_clo_s8a,                       ! o
     1  s8a_k,istat_s8a,comment_s8a,                                     ! o
     1  s3a_k,istat_s3a,comment_s3a,                                     ! o
     1  sst_k,istat_sst,comment_sst,                                     ! o
     1  cldtop_co2_pa_a,cloud_frac_co2_a,istat_co2,lstat_co2_a)          ! o

        use mem_namelist, only: l_use_s8a

!       input
        real lat(imax,jmax)
        real lon(imax,jmax)

!       output
        real s8a_k(imax,jmax)
        real s3a_k(imax,jmax)
        real sst_k(imax,jmax)
        real cldtop_co2_pa_a(imax,jmax)
        real cloud_frac_co2_a(imax,jmax)

        logical lstat_co2_a(imax,jmax)
        logical l_use_39, l_use_co2

!       local
        real pct_pa(imax,jmax)
        real lca(imax,jmax)
        real    subpoint_lat_clo_s8a(imax,jmax)
        real    subpoint_lon_clo_s8a(imax,jmax)
        real    subpoint_lat_clo_s3a(imax,jmax)
        real    subpoint_lon_clo_s3a(imax,jmax)

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character*31 ext
        character var*3,units*10
        character*125 comment_s8a,comment_s3a,comment_sst,comment

        write(6,*)' subroutine get_sat_data...'

        write(6,*)' getting ir satellite data from lvd file'
        ext = lvd_ext
        var = 's8a'
        ilevel = 0
        if(l_use_s8a .eqv. .true.)then
          call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window    
     1                     ,i4time_s8a,lat,lon
     1                     ,subpoint_lat_clo_s8a,subpoint_lon_clo_s8a     ! o 
     1                     ,ext,var,units
     1                     ,comment_s8a,imax,jmax,s8a_k,ilevel
     1                     ,istat_s8a)      
          if(istat_s8a .ne. 1)then
            write(6,*)' no s8a data available, i4_sat_window is: '
     1               ,i4_sat_window
          endif
        else
          write(6,*)' skipping read of s8a satellite data'
          istat_s8a = 0
        endif

!       final qc check on band 8 (11.2 mm) brightness temps
!       hopefully, bad/missing  values were filtered out in the creation of lvd
!       any remaining bad/missing pixels will be evaluated in this qc check
        icount = 0
        do j = 1,jmax
        do i = 1,imax
            if(s8a_k(i,j) .lt. 173. .or. s8a_k(i,j) .gt. 350.
     1                              .or. istat_s8a  .ne. 1      )then
                if(icount .le. 100 .and. istat_s8a .eq. 1)then
                    write(6,*)' bad lvd/s8a satellite brightness '       
     1                       ,'temperature of'
     1                       ,s8a_k(i,j),' at',i,j
                endif
                icount = icount + 1
                s8a_k(i,j) = r_missing_data
            endif
        enddo
        enddo


!       obtain sea surface temperature
        write(6,*)' getting sea sfc temps data from sst file'
        ext = 'sst'
        var = 'sst'
        ilevel = 0
        call get_laps_2dgrid(i4time,3600,i4time_nearest,ext,var
     1                ,units,comment_sst,imax,jmax,sst_k,ilevel
     1                ,istat_sst)       
        if(istat_sst .ne. 1)then
            write(6,*)' note: cannot read sst_k'
        endif

        if(l_use_39)then
          write(6,'(" getting 3.9u satellite data from lvd file")')
          ext = lvd_ext
          var = 's3a'
          ilevel = 0
          call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window      
     1                       ,i4time_nearest,lat,lon
     1                       ,subpoint_lat_clo_s3a           ! o
     1                       ,subpoint_lon_clo_s3a           ! o 
     1                       ,ext,var,units
     1                       ,comment_s3a,imax,jmax,s3a_k,ilevel
     1                       ,istat_s3a)
          if(istat_s3a .ne. 1)then
              write(6,*)' no s3a data available'
              s3a_k = r_missing_data
          endif

        else
          write(6,'(" namelist flag set to not use 3.9u (s3a) data")')
          istat_s3a = 0
          s3a_k = r_missing_data

        endif

!       obtain nesdis cloud-top pressure
        if(l_use_co2)then ! this is the "mode 2" type of co2 usage
            i4_co2_window = latency_co2

            write(6,*)' getting nesdis cloud-top pressure'
            ext = 'ctp'
            var = 'pct'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                      ,ext,var       
     1                      ,units,comment,imax,jmax,pct_pa,ilevel
     1                      ,istat_pct)       
            if(abs(istat_pct) .ne. 1)then
                write(6,*)' note: cannot read nesdis cloud-top pressure'       
            endif

!           obtain nesdis cloud-fraction
            write(6,*)' getting nesdis cloud-fraction'
            ext = 'ctp'
            var = 'lca'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                          ,ext,var
     1                          ,units,comment,imax,jmax,lca,ilevel
     1                          ,istat_lca)       
            if(abs(istat_lca) .ne. 1)then
                write(6,*)' note: cannot read nesdis cloud-fraction'
            endif

        endif ! l_use_co2

!       calculate co2-slicing cloud-top pressure

        cloud_frac_co2_a = r_missing_data
        cldtop_co2_pa_a  = r_missing_data
        lstat_co2_a = .false.
        icount = 0

        if(abs(istat_pct) .eq. 1 .and. abs(istat_lca) .eq. 1 
     1                           .and. l_use_co2                )then
            write(6,*)' extracting co2-slicing info from nesdis data'
            do j = 1,jmax
            do i = 1,imax
!               test for partial cloudiness
                if(lca(i,j) .gt. 0. .and. lca(i,j) .lt. 1.0)then ! use co2 data
                    if(pct_pa(i,j) .lt. 100000.)then
                        icount = icount + 1
                        cloud_frac_co2_a(i,j) = lca(i,j)
                        cldtop_co2_pa_a(i,j) = pct_pa(i,j) 
                        lstat_co2_a(i,j) = .true.
                    endif
                endif
            enddo ! i
            enddo ! j

            istat_co2 = 0 ! for now
        else
            istat_co2 = 0
        endif

        write(6,*)' number of utilized co2-slicing image points = '
     1           ,icount

        if(l_use_co2)then
            percent_co2_pot = float(icount) / float(imax*jmax) * 100.

            write(6,101)percent_co2_pot
101         format(' co2-slicing imagery potentially used over ',f6.2
     1            ,'% of domain')
        endif

        return
        end
