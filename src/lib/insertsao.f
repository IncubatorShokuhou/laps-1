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
        subroutine insert_sao(i4time,cldcv,cf_modelfg,t_modelfg,cld_hts   ! i
     1          ,default_clear_cover,namelist_parms                       ! i
     1          ,lat,lon,topo,t_sfc_k,wtcldcv                             ! i
     1          ,name_array,l_perimeter,ista_snd,cvr_snd,cld_snd
     1          ,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd
     1          ,ni,nj,nk                                                 ! i
     1          ,n_obs_b,lat_sta_ret,lon_sta_ret,c_stations
     1          ,wx,t,td                                                  ! o
     1          ,elev                                                     ! o
     1          ,rad,solar_ea                                             ! o
     1          ,istatus                                                  ! o
     1          ,maxstns,ix_low,ix_high,iy_low,iy_high)

!       1995 steve albers                         original version
!       11-nov-1995    steve albers     ignore saos reporting 'x'
!       12-nov-1995    steve albers     max cloud cover = 1.0, not 1.01
!       10-oct-1996    steve albers     max cloud cover = 1.0, not 1.01
!                                           (completed the change)
!       12-nov-1996    steve albers     improved logging and cleanup
!        6-dec-1996    steve albers     filter the obs string
!        1-aug-1997    ken dritz        added i_perimeter and maxstns as dummy
!                                       arguments
!        1-aug-1997    ken dritz        added maxstns, ix_low, ix_high, iy_low,
!                                       and iy_high as dummy arguments.
!        1-aug-1997    ken dritz        removed parameter statements for
!                                       ix_low, ix_high, iy_low, and iy_high.
!        1-aug-1997    ken dritz        removed include of lapsparms.for
!        6-aug-1997    steve albers     removed equivalences.

        include 'cloud.inc'

        character*150 c150_filename,directory
        character*31 ext
        character*13 filename13

!       arrays for reading in sao cloud data
        character*(*) c_stations(maxstns)
        character*5 c5_outstring

        real lat_sta_ret(maxstns)
        real lon_sta_ret(maxstns)

        logical l_sao_lso
        data l_sao_lso /.true./ ! do things the new way?
        logical l_perimeter
        character*3 lso_ext
        data lso_ext /'lso'/

        logical l_dry, l_parse

!       arrays for reading in the sao data from the lso files
        real   elev(maxstns),t(maxstns),td(maxstns)
        real   rad(maxstns),solar_ea(maxstns)
        real   ht_base_ret(maxstns,5)
c
        integer   n_cloud_layers_ret(maxstns)
        integer   obstime(maxstns)
c
        character   obstype(maxstns)*6,atype(maxstns)*6
     1             ,wx(maxstns)*8
        character   amt_ret(maxstns,5)*4

        character*8 c8_project

!       arrays for inserting the cloud data into the laps grid
        real cldcv(ni,nj,nk)
        real cf_modelfg(ni,nj,nk)
        real t_modelfg(ni,nj,nk)
        real topo(ni,nj)
        real t_sfc_k(ni,nj)
        real wtcldcv(ni,nj,nk)
        real cld_hts(nk)
        real lat(ni,nj),lon(ni,nj)
        character*1 name_array(ni,nj,nk)

!       arrays for cloud soundings
        integer ista_snd(max_cld_snd)
        real cld_snd(max_cld_snd,nk)
        real wt_snd(max_cld_snd,nk)
        real cvr_snd(max_cld_snd)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        logical l_out_of_grid

        call get_sfc_badflag(badflag,istatus)

!       initialize ista_snd
        do i = 1,max_cld_snd
            ista_snd(i) = 0
        enddo ! i

!       i4time_database = (i4time / 3600) * 3600
        i4time_database = i4time

!       this section calls a routine 'get_sao' which interfaces with the raw
!       fsl sao data in order to create an intermediate cloud layer file (.sao).
!       in the future, this task could be performed elsewhere within laps

!       construct file name for lso file
        ext = lso_ext
        call get_directory(ext,directory,len_dir) ! returns directory
        c150_filename = directory(1:len_dir)
     1                            //filename13(i4time_database,lso_ext)

        n_wait_saos = 0   ! (0,2) this is inputted, 10 min per potential wait

10      write(6,*)' calling sao ingest routine'

!       access sao data from lso files
        call read_surface_sa(i4time_database,maxstns,              ! i
     1   n_obs_b,c_stations,obstype,atype,                         ! o
     1   lat_sta_ret,lon_sta_ret,elev,wx,t,td,                     ! o
     1   n_cloud_layers_ret,amt_ret,ht_base_ret,                   ! o
     1   rad,solar_ea,obstime,istatus)                             ! o

        if(istatus .ne. 1)then
            write(6,*)' bad status returned from reading sao data'
            return
        endif

        if(.not. namelist_parms%l_use_metars)then
            write(6,*)' skipping insertion of metar/saos for analysis'       
            write(6,*)' metar/sao data used for verification only'
            return
        endif

        write(6,*)' # of obs (box) = ',n_obs_b

        write(6,*)' now we are looping to insert the stations'

        n_analyzed = 0

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)return

!       loop through the stations
        do i=1,n_obs_b

          call filter_string(obstype(i))

          if(l_parse(c8_project,'afgwc') .or. 
     1       l_parse(c8_project,'afwa')       )then
              atype(i)='u'//atype(i)(2:6)
          endif

!         determine whether we want to analyze cloud layers from this station
          if(obstype(i)(1:4) .eq. 'meso')goto126
          if(obstype(i)(1:4) .eq. 'cdot')goto126

          write(6,*)

c place station at proper laps grid point
          call latlon_to_rlapsgrid(lat_sta_ret(i),lon_sta_ret(i)
     1                                  ,lat,lon,ni,nj,ri,rj,istatus)

          ilaps = nint(ri)
          jlaps = nint(rj)

          if(  ilaps .lt. ix_low .or. ilaps .gt. ix_high
     1    .or. jlaps .lt. iy_low .or. jlaps .gt. iy_high)then
              write(6,*)' note: out of bounds ',c_stations(i)
              goto 126
          endif

          if(n_cloud_layers_ret(i) .eq. 0)then ! kick out amos stations not
                                               ! reporting clouds this time
              if(i .le. 1000)then
                  write(6,*)' no cloud layers reported - '
     1              ,'clr?/msg?/amos? - goto end of loop '       
     1              ,c_stations(i),' ',obstype(i),' ',atype(i)
              endif
              goto 126
          endif

!         what is the height limit of the cloud observation? is obstype valid?
          if(      obstype(i)(1:5) .eq. 'metar'
     1        .or. obstype(i)(1:5) .eq. 'testm' 
     1        .or. obstype(i)(1:5) .eq. 'speci' 
     1        .or. obstype(i)(1:5) .eq. 'tests' 
     1        .or. obstype(i)(1:5) .eq. 'synop' )then  ! new lso file format

              if(  atype(i)(1:1) .eq. 'a'
     1        .or. atype(i)(1:1) .eq. 'u' )then      ! use 12000' limit

                  if(atype(i)(1:1) .eq. 'a')then     ! automated station
                      i_auto = 1                       
                  else                                 ! indeterminate
                      i_auto = 0
                  endif

                  ht_defined = elev(i) + 12000./3.281

              else                                     ! non-automated
                  i_auto = -1
                  ht_defined = 99999.

              endif ! automated station

          else                                         ! non-sanctioned cld type
              write(6,*)' warning, questionable obstype having '
     1                 ,'cloud layers - reject: ',obstype(i)
     1                 ,' ',c_stations(i)

              goto 126                                 ! loop to next station

!             if(obstype(i)(5:5) .ne. ' ' .and.
!    1           obstype(i)(4:7) .ne. 'amos')then      ! automated station 
!                                                      ! (12000' limit)
!                 ht_defined = elev(i) + 12000./3.281
!             else                                     ! non-automated
!                 ht_defined = 99999.
!             endif

          endif ! sanctioned obstype for reporting cloud layers

          n_analyzed = n_analyzed + 1
          n_cld_snd = n_cld_snd + 1
          ista_snd(n_cld_snd) = i

          c5_outstring = c_stations(i)

          write(6,1,err=110)c5_outstring,n_cld_snd,lat_sta_ret(i)
     1         ,lon_sta_ret(i),n_cloud_layers_ret(i)
     1         ,ilaps,jlaps,obstype(i),atype(i),ht_defined ! ,obstime(i)
1         format(1x,a5,1x,i4,2f8.2,i3,2i4,1x,a8,1x,a6,f8.0,i5)

110       do l = 1,n_cloud_layers_ret(i)
              write(6,2,err=3)amt_ret(i,l),ht_base_ret(i,l)
2             format(1x,a4,f8.0)
3             continue
          enddo ! l

          if(  ilaps .lt. 1 .or. ilaps .gt. ni
     1    .or. jlaps .lt. 1 .or. jlaps .gt. nj)then
              l_out_of_grid = .true.
          else
              l_out_of_grid = .false.
              name_array(ilaps,jlaps,:)=c_stations(i)(1:1)
          endif

          cvr_snd(n_cld_snd) = 0.

!         initialize summation layer calculation for this ob
          call get_layer_cover(0.,cover,istatus)

          do l=1,n_cloud_layers_ret(i)

              cover = 0.

              ht_base=ht_base_ret(i,l) ! meters msl

              if(ht_base .gt. ht_defined+1.)then

                if( (.not. l_parse(amt_ret(i,l),'clr') ) .and.
     1              (.not. l_parse(amt_ret(i,l),'skc') )      )then ! clouds

                  if(.true.)then ! allow a redefinition of ht_defined
                    write(6,*)' warning, inconsistent sao data,'
     1              //' cloud base is'       
     1              //' reported to be too high for this sensor type'       

                    write(6,*)ht_base,ht_defined,' ',obstype(i)
     1                                          ,' ',atype(i)

                    write(6,*)' please check cloud layer heights in the'       
     1              //' lso file to see that they are compatable with'     
     1              //' the types of sensors used.'

                    write(6,*)' assume human augmented ob, raising '
     1                       ,'ht_defined'
                    ht_defined = ht_base

                  else ! flag qc error if ht_defined is exceeded by cloud layer
                    write(6,*)
     1                ' error, inconsistent sao data, cloud base is'      
     1              //' reported to be too high for this sensor type'       

                    write(6,*)ht_base,ht_defined,' ',obstype(i)
     1                                          ,' ',atype(i)
                    write(6,*)
     1                ' please check cloud layer heights in the lso'
     1              //' file to see that they are compatable with the'       
     1              //' types of sensors used.'

                    istatus = 0
                    return

                  endif

                else ! clr
                  write(6,*)' warning, clr sky cloud base does not'
     1            //' reflect this sensor type'
                  write(6,*)ht_base,ht_defined,' ',obstype(i)
     1                                        ,' ',atype(i)

                endif ! clouds

              endif ! ht_base > ht_defined

              if(l_parse(amt_ret(i,l),'vv'))then
                  if(l .eq. 1)then
                      write(6,*)
     1                    ' vv report, reset ht_defined to ht_base'
     1                    ,ht_defined,ht_base
                      ht_defined = ht_base
                  else
                      write(6,*)' warning: vv report with >1 layer ob'
                  endif
              endif

c clouds are now in msl
!             fill in clear for entire column for metar or up to ht_base for 
!             awos or up to ht_base for vv.
              if(l_parse(amt_ret(i,l),'clr')    .or.
     1           l_parse(amt_ret(i,l),'skc')    .or.
     1           l_parse(amt_ret(i,l),'vv')         )then

                  cover=default_clear_cover
                  if(istatus .ne. 1)goto125 ! go to next station

                  do k=1,nk
                      if(     cld_hts(k).le.ht_base
     1                  .and. cld_hts(k).le.ht_defined      )then
                          call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                      endif
                  enddo

                  write(6,*)' filled in ',amt_ret(i,l),' from bottom'
     1                     ,' of domain up to '
     1                     ,nint(min(ht_base,ht_defined)),' meters'     

                  if(.true.)then                               ! skc is now used
                      if(l_parse(amt_ret(i,l),'clr') .and. 
     1                                            i_auto .eq. -1)then       
                          write(6,*)' warning: clr reported for '
     1                             ,'non-automated station'

                      elseif(l_parse(amt_ret(i,l),'skc') .and. ! converse
     1                                            i_auto .eq. +1)then
                          write(6,*)' warning: skc reported for '
     1                             ,'automated station'

                      endif ! clr/skc test

                  endif ! .true.

                  if(n_cloud_layers_ret(i) .gt. 1)then
                      write(6,*)' warning: clr/skc/vv in ob that'
     1                         ,' has more than one layer'
                  endif
!                 go to 125 ! loop to next station

              endif


!             if this station has obscured (but not thin obscured),
!             leave the entire cloud sounding to say "missing data".
!             thin obscured simply drops through as an ignored cloud layer
              do nc = 1,4 ! parse the cloud amount string
                  if(amt_ret(i,l)(nc:nc) .eq. 'x')then ! obscured or thin obsc
                      if(nc .eq. 1)then
                          write(6,*)' obscured sky detected'
                          goto 125
                      else ! search for "thin" designation
                          if(amt_ret(i,l)(nc-1:nc-1) .ne. '-')then
                              write(6,*)' obscured sky detected'
                              goto 125
                          endif
                      endif
                  endif
              enddo

!             lcl check
              if(t(i) .ne. badflag .and. td(i) .ne. badflag)then
                  ht_lcl_agl = (t(i) - td(i)) * 250. * .3048 ! meters
                  if(td(i) .gt. 32.)then
                      lcl_thresh = 300.
                  else
                      lcl_thresh = 600.
                  endif
              else
                  ht_lcl_agl = 0. ! used for missing value
                  lcl_thresh = 0.   
              endif

              if((ht_base-elev(i)) .lt. (ht_lcl_agl - lcl_thresh))then
                  write(6,111)t(i),td(i),ht_base-elev(i),ht_lcl_agl             
111               format(' warning: ht_base_agl << ht_lcl_agl (m) ' 
     1                                                   ,2f8.1,2f9.1)       
                  goto 126                               ! loop to next station

              elseif(ht_base .lt. ht_lcl_agl)then
                  write(6,112)t(i),td(i),ht_base,ht_lcl_agl                    
112               format(' note: ht_base < ht_lcl_agl (m) ',2f8.1,2f9.1)

              else
                  write(6,113,err=114)t(i),td(i),ht_base,ht_lcl_agl 
113               format(' t,td,ht_base, ht_lcl (m) ',2f8.1,2f9.1)
114               continue

              endif

              if(l_parse(amt_ret(i,l),'few') .and.
     1           n_cloud_layers_ret(i) .eq. 1      )then
                  summation_cover=.125
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+1000.
                  do k=1,nk

                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif

              if(l_parse(amt_ret(i,l),'sct'))then
                  summation_cover=.44
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+1000.
                  do k=1,nk

                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif

              if(l_parse(amt_ret(i,l),'-bkn'))then
                  summation_cover=.4 ! .5
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+1000.
                  do k=1,nk

                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif

              if(l_parse(amt_ret(i,l),'bkn'))then
                  summation_cover=.75
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+cld_thk(ht_base) ! 1500.

                  do k=1,nk

                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif

              if(l_parse(amt_ret(i,l),'-ovc'))then
                  summation_cover=.6 ! .9
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+1000.
                  do k=1,nk

                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif

              if(l_parse(amt_ret(i,l),'ovc'))then
                  summation_cover=1.00 
                  call get_layer_cover(summation_cover,cover,istatus)      
                  if(istatus .ne. 1)goto125 ! go to next station

                  ht_top=ht_base+cld_thk(ht_base) ! 1500.

                  do k=1,nk
                      if(cld_hts(k).ge.ht_base .and. 
     1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,0,l_dry)

                          if(.not. l_dry)then
                              call spread2(
     1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
                          endif

                      else
!                         initialize the modify sounding routine
                          call modify_sounding(
     1                         cld_snd,n_cld_snd,max_cld_snd
     1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
     1                        ,ht_base,ht_top,1,l_dry)

                      endif
                  enddo
              endif


!             we should not hit this anymore
              if(l_parse(amt_ret(i,l),'x'))then
                  write(6,*)' error: insertsao - stop x'
                  if(.true.)stop

!                 summation_cover=1.00 
!                 call get_layer_cover(summation_cover,cover,istatus)      
!                 if(istatus .ne. 1)goto125 ! go to next station

!                 ht_top=ht_base+cld_thk(ht_base) ! 1500.

!                 do k=1,nk
!                     if(cld_hts(k).ge.ht_base .and. 
!    1                   cld_hts(k).le.ht_top        )then

!                         search for model d(cldcv)/dz within cloud layer
!                         call modify_sounding(
!    1                         cld_snd,n_cld_snd,max_cld_snd
!    1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
!    1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
!    1                        ,ht_base,ht_top,0,l_dry)

!                         if(.not. l_dry)then
!                             call spread2(
!    1                           cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
!    1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
!                         endif

!                     else
!                         initialize the modify sounding routine
!                         call modify_sounding(
!    1                         cld_snd,n_cld_snd,max_cld_snd
!    1                        ,cf_modelfg,t_modelfg,topo,t_sfc_k
!    1                        ,ilaps,jlaps,k,ni,nj,nk,cld_hts
!    1                        ,ht_base,ht_top,1,l_dry)

!                     endif
!                 enddo
              endif ! amt_ret

!             calculate summation cover assuming obs are in layer cover
!             cvr_snd(n_cld_snd) = 1. - ((1. - cvr_snd(n_cld_snd)) 
!    1                           * cover)

!             obtain summation cover directly from obs
              cvr_snd(n_cld_snd) = summation_cover

 115          continue

        enddo ! l (cloud layer)

!       locate the highest ceiling
        k_ceil = nk
        if(l_perimeter)then
            do k=nk,1,-1
                if(wt_snd(n_cld_snd,k) .eq. 1.00 .and.
     1           cld_snd(n_cld_snd,k) .gt. 0.5          )then
                    k_ceil = k
                    goto 1001
                endif
            enddo
        else
            do k=nk,1,-1
                if(wtcldcv(ilaps,jlaps,k) .eq. 1.00 .and.
     1           cldcv(ilaps,jlaps,k) .gt. 0.5          )then
                    k_ceil = k
                    goto 1001
                endif
            enddo
        endif


!       fill in other clear layers outside of clouds, below the ceiling,
!                        and within defined height range of sensor.
1001    cover = default_clear_cover

        ht_fill = min(ht_defined,cld_hts(k_ceil))

        write(6,*)' k_ceil/ht_ceil/ht_defined/ht_fill = ',k_ceil
     1           ,nint(cld_hts(k_ceil)),nint(ht_defined),nint(ht_fill)

        if(l_perimeter)then
          do k=1,k_ceil
            if(     wt_snd(n_cld_snd,k) .ne.  1.00
     1        .and. cld_hts(k)          .le.  ht_defined          )then
                call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
            endif
          enddo
        else
          do k=1,k_ceil
            if(     wtcldcv(ilaps,jlaps,k) .ne.  1.00
     1          .and. cld_hts(k)           .le.  ht_defined       )then
                call spread2(cld_snd,wt_snd,i_snd,j_snd,n_cld_snd
     1                          ,max_cld_snd,nk,ilaps,jlaps,k,cover,1.)
            endif
          enddo
        endif

125     continue

!       clear out name array for the station where the station has no influence
        if(  ilaps .lt. 1 .or. ilaps .gt. ni
     1  .or. jlaps .lt. 1 .or. jlaps .gt. nj)then
            l_out_of_grid = .true.
        else
            l_out_of_grid = .false.
        endif

        if(l_perimeter .and. .not. l_out_of_grid)then
          do k=1,nk
            if(wt_snd(n_cld_snd,k) .ne. 1.00)then
                name_array(ilaps,jlaps,k) = ' '
            endif
          enddo
        endif

 126    continue

        enddo ! i

        write(6,*)
        write(6,*)' num stations analyzed/cloud soundings = '
     1                          ,n_analyzed,n_cld_snd

999     continue

        istatus = 1

        return
        end

        subroutine modify_sounding(cld_snd,n_cld_snd,max_cld_snd
     1          ,cf_modelfg,t_modelfg,topo,t_sfc_k
     1          ,i_in,j_in,k,ni,nj,nk,cld_hts
     1          ,ht_base,ht_top,init,l_dry)

        real cld_snd(max_cld_snd,nk)
        real cf_modelfg(ni,nj,nk)
        real t_modelfg(ni,nj,nk)
        real topo(ni,nj)
        real t_sfc_k(ni,nj)
        real cld_hts(nk)

        logical l_wait_for_base,l_dry,l_cf,l_inversion
        save l_wait_for_base,cf_model_base,t_model_base,l_inversion,t_su
     1bcloud

        l_dry = .false.
        l_cf = .false.

!       find laps grid point nearest the sao if it is out of bounds
        i = max(min(i_in,ni),1)
        j = max(min(j_in,nj),1)

        if(init .eq. 1)then ! (below base )
                            ! reset to wait for the beginning of the layer
            l_wait_for_base = .true.
            l_inversion = .false.
            t_subcloud = t_modelfg(i,j,k)
            return

        else ! init = 0 (inside cloud layer)
            if(l_wait_for_base)then ! set reference (just within cloud base)
                l_wait_for_base = .false.
                cf_model_base = cf_modelfg(i,j,k)
                t_model_base = t_modelfg(i,j,k)

                write(6,21)t_subcloud
21              format(' modify_sounding.....          '
     1         ,'cf     t    dlt th r i   i   j kcld m-msl'
     1         /' model    t    subcloud   = ',7x,f7.2)

!               write(6,1)cf_model_base,t_modelfg(i,j,k),l_cf
!       1                               ,l_inversion,i,j,k,nint(cld_hts(k))
1               format(' model rh/t at cloud base = ',2f7.2,2l2,3i4,i6)

            endif

            if(.true.)then ! determine if cloud should be cleared out

!               set inversion strength flag
                t_dry_adiabat = t_sfc_k(i,j)
     1                     -.0098 * (cld_hts(k) - topo(i,j))
                t_inversion_strength = t_modelfg(i,j,k) - t_dry_adiabat

                if(
     1     (    (t_modelfg(i,j,k) .gt. t_model_base)
     1                                .or.
     1            (t_modelfg(i,j,k) .gt. t_subcloud .and. k .ge. 2) )
     1                         .and.
     1                (t_modelfg(i,j,k) .gt. 283.15)       ! temp check
     1                         .and.
     1                (t_inversion_strength .gt. 4.)       ! delta theta chk
     1                                              )then  ! inversion search
                    l_inversion = .true.
                    write(6,2)cf_modelfg(i,j,k),t_modelfg(i,j,k)
     1               ,t_inversion_strength,l_cf,l_inversion
     1                       ,i,j,k,nint(cld_hts(k))
2                   format(' inversion detected       = '
     1                                       ,3f7.2,2l2,3i4,i6)
                elseif(cf_modelfg(i,j,k) .lt. cf_model_base - 0.3   ! cf search
     1            .and.    cld_hts(k) - ht_base .ge. 500.)then
                    l_cf = .true.
                    write(6,3)cf_modelfg(i,j,k),t_modelfg(i,j,k)
     1               ,t_inversion_strength,l_cf,l_inversion
     1                       ,i,j,k,nint(cld_hts(k))
3                   format(' dry layer detected       = '
     1                                       ,3f7.2,2l2,3i4,i6)
                else                                        ! not newly flagged
                    write(6,4)cf_modelfg(i,j,k),t_modelfg(i,j,k)
     1               ,t_inversion_strength,l_cf,l_inversion
     1                       ,i,j,k,nint(cld_hts(k))
4                   format(' model rh/t in cloud      = ',3f7.2,2l2,3i4,
     1i6)
                endif

                if(l_cf .or. l_inversion)then
!               if(l_cf)then
                    l_dry = .true.
                endif

            endif
        endif

        return
        end

        function cld_thk(ht_base)

        if(ht_base .gt. 9000.)then
            cld_thk = 2500.
        elseif(ht_base .gt. 6000.)then
            cld_thk = 1500.
        elseif(ht_base .gt. 4000.)then
            cld_thk = 1000.
        else
            cld_thk = 600.
        endif

        return
        end


        subroutine get_layer_cover(summation_cover,layer_cover,istatus)        
!                                         i             o         o

        real layer_cover

        save summation_cover_last

!       null the changes for testing
!       layer_cover = summation_cover 
!       istatus = 1
!       return

        if(summation_cover .eq. 0.)then ! initializing a new cloud ob
            goto900
        endif

        if(summation_cover_last .ge. 1.)then
            write(6,*)' ob error, previous summation cover already ovc'       
            write(6,*)summation_cover_last
            istatus = 0
            return
        endif

        layer_cover = (summation_cover - summation_cover_last)
     1              / (1. - summation_cover_last)

        if(layer_cover .le. 0.)then
            if(layer_cover .lt. 0.)then
                write(6,*)' warning: layer_cover < 0, decrease '
     1                   ,'observed in summation cover: '
     1                   ,summation_cover_last,summation_cover
            endif

            layer_cover = .125
            write(6,*)' setting layer_cover to minimum value of '
     1               ,layer_cover 

        endif

        write(6,1)summation_cover,layer_cover
 1      format(' get_layer_cover: summation/layer ',2f8.3)

 900    summation_cover_last = summation_cover

        istatus = 1
        return
        end
