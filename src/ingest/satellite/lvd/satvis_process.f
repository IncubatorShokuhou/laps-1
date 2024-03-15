      subroutine process_vis_satellite(csatid,
     &           c_sat_type,
     &           i4time,                                           ! i
     &           ni,nj,lat,lon,
     &           n_vis_lines,n_vis_elem,
     &           r_grid_ratio,
     &           image_vis,                                        ! i
     &           r_llij_lut_ri,
     &           r_llij_lut_rj,
     &           sublat_d,sublon_d,range_m,
     &           laps_vis_raw,laps_vis_refl,laps_vis_norm,albedo,  ! o
     &           istatus)
c
c**************************************************************************
c
c       routine to collect satellite data for the laps analyses.
c
c       changes:
c        p.a. stamus       12-01-92       original (from 'get_vas_bt' routine).
c                     01-11-93       write output in laps standard format.
c                     02-01-93       add snooze call for vis data.
c                     09-16-93       add bands 3,4,5,12 data.
c       j.r. smart    03-01-94	     implement on the sun.  adapt to ispan
c                             grids. this required removing all references to
c                             goes mdals/ground station satellite receiving.
c          "          10-26-94       modified for goes-8 data. no need to compute
c                                    radiance and btemp as this conversion is included
c                                    in icnt_lut.
c          "          11-28-95       extract the visible processing part as a separate
c                                    subroutine.
c          "          10-21-96       incorporated new normalize brightness procedure to
c				     account for vis data that has been pre-normalized (ie., gvar fsl-conus)
c          "          03-07-97      new normalize brightness routine that accounts for missing data.
c          "          03-12-97      added csatid to argument list. this directs the normalization procedure.
c          "          07-06-99      added emission_angle_d as output from normalize_brightness and input
c                                   to vis_2_albedo.
c
c       notes:
c       this program processes vis satellite data from ispan database
c
c       variables:
c       image_vis        ra       i       visible (on satellite grid)
c       laps_vis         ra       o       visible (raw)
c       vis_norm         ra       o       visible (normalized)
c       albedo           ra       o       albedo (0.0 -- 1.0)
c
c
c       note: for details on the filtering and averaging methods, see the
c       vasdat2 routine or talk to s. albers.
c
c****************************************************************************
c
       implicit none
c
       integer ni,nj
c
c..... grids to put the satellite data on.
c
       real laps_vis_raw(ni,nj)
       real laps_vis_refl(ni,nj)
       real laps_vis_norm(ni,nj)
       real laps_vis_norm_natl(ni,nj)
       real albedo(ni,nj)
       real phase_angle_d(ni,nj)
       real specular_ref_angle_d(ni,nj)
       real emission_angle_d(ni,nj)
       real rland_frac(ni,nj)
c 
c..... laps lat/lon files.
c
       real lat(ni,nj)
       real lon(ni,nj)
c
       integer n_vis_lines,n_vis_elem
       integer i_dir
       integer n_missing_albedo
       integer len_dir
       integer isave,jsave
       integer ismax,jsmax
       integer ismin,jsmin

       real r_llij_lut_ri(ni,nj)
       real r_llij_lut_rj(ni,nj)
       real r_grid_ratio, bad_frac
       real image_vis(n_vis_elem,n_vis_lines)
       real albedo_max,albedo_min
       real rspacing_dum
       real r_missing_data,rmin,rmax
c
       integer istatus_a
       integer istatus_l
       integer istatus_m
       integer istatus_n
       integer istatus_v
       integer istatus_r
       integer istatus(3)
       integer i,j
c      integer imn,jmn,imx,jmx

       real sublat_d,sublon_d,range_m
c      real difference
c      real diffsum
c      real diffsum_abs
c      real meandiff
c      real meanabsdiff
c      real maxdiff
c      real mindiff
       real ave,adev,sdev,var,skew,curt

       real   visin1_g8,visin2_g8
       real   visout1_g8,visout2_g8
       real   visin1_cur,visin2_cur
       real   visout1_cur,visout2_cur

       integer   icnt
c
       integer i4time,imax,jmax
       integer isat
c
       character*255 dataroot
       character*3   var_2d
       character*3   c_sat_type
       character*6   csatid
       character*150  directory
       character*31  ext
       character*40  domain_name
       character*10  units_2d
       character*125 comment_2d
       character*9 a9time

       integer mode_norm

       include 'satellite_dims_lvd.inc'
       include 'satellite_common_lvd.inc'
c
c start
c ----------------------
       istatus(1) = 0
       istatus(2) = 0
       istatus(3) = 0

       imax = ni
       jmax = nj

       call get_r_missing_data(r_missing_data,istatus_r)
       if(istatus_r.ne.1)then
          write(6,*)'error getting r_missing_data'
          goto 999
       endif
c
c call the visible satellite data remapping subroutine.
c
       write(6,*)
       write(6,*)' processing visible satellite data - ',c_sat_type
       write(6,*)' ---------------------------------------'

       call array_range(image_vis,n_vis_elem,n_vis_lines,rmin,rmax
     1                 ,r_missing_data)
       write(6,*)'range of non-missing vis on sat grid is',rmin,rmax
       write(6,*)'vis sat grid (center) = '
     1          ,image_vis(n_vis_elem/2,n_vis_lines/2)


       call satdat2laps_vis(
     &                  r_grid_ratio,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  imax,jmax,
     &                  n_vis_lines,n_vis_elem, ! image_vis array dimensions
     &                  image_vis,              ! satellite grid
     &                  laps_vis_raw,           ! model grid
     &                  istatus_v)
       if(istatus_v .ne. 1) then
          write(*,921)istatus_v
          istatus(1) = -1
       end if
 921   format(1x,' +++ warning. bad istatus_v = ',i3, 'from
     & satdat2laps_vis+++')
c
       call check(laps_vis_raw,r_missing_data,istatus_v,imax,jmax)
       if(istatus_v .lt. 0) then
          write(6,*)'bad vis counts detected in process_vis_satellite'
          bad_frac = -float(istatus_v) / float(imax*jmax)
          write(6,916) istatus_v,bad_frac
          istatus(1) = istatus_v
 916      format(' +++ warning. visible status/frac = ',i8,f7.3)
       else
          write(6,*)
     1         'all vis data ok on model grid (process_vis_satellite)'         
       endif

       call array_range(laps_vis_raw,imax,jmax,rmin,rmax,r_missing_data)
       write(6,*)' range of non-missing vis_raw on model grid is'
     1          ,rmin,rmax

       do j=1,jmax
       do i=1,imax
          laps_vis_norm(i,j)=laps_vis_raw(i,j)
       enddo
       enddo

       write(6,*)'vis_raw mdl grid (center) = '
     1          ,laps_vis_raw(imax/2,jmax/2)
c
c.....       normalize the vis data.
c
c
c our standard is goes08 and we know how to stretch this
c to look like goes07; therefore we should always need the
c goes8-to-goes7 stretch parameters:  when processing goes08
c data we only need to stretch once with the g8 stretch parms;
c all other satellite vis data require two stretches.

       visin1_g8=0.0    !vis_cnt_range_in(1,1)
       visin2_g8=303.57 !vis_cnt_range_in(2,1)
       visout1_g8=0.0   !vis_cnt_range_out(1,1)
       visout2_g8=255.0 !vis_cnt_range_out(2,1)

c for locally produced (ground station sbn look-alike) visible we
c want to reverse the already applied (national scale) normalization 
c and then "re-normalized" for the local domain.  all satellites of
c type 'cdf' qualify for the reverse. only one satellite is selected
c for type 'cdf' by gsd's its group.
c
c currently (2007) type 'cdf' is goes12.
c

!      determine reflectance
       laps_vis_refl(:,:) = r_missing_data
       if(trim(csatid) .eq. 'coms' .or. c_sat_type .eq. 'gnp' .or.
     1      c_sat_type .eq. 'gr2'  .or. c_sat_type .eq. 'jma'      )then
         write(6,*)' scaling 1.0 raw data = 1.0 reflectance'
         where(laps_vis_raw(:,:) .ne. r_missing_data)
             laps_vis_refl(:,:) = 1.0 * laps_vis_raw(:,:) / 1.
         endwhere
       else
         write(6,*)' scaling 255 vis raw counts = 1.2 reflectance'
         where(laps_vis_raw(:,:) .ne. r_missing_data)
             laps_vis_refl(:,:) = 1.2 * laps_vis_raw(:,:) / 255.
         endwhere
       endif

       if(c_sat_type.eq.'cdf')then

          if(csatid.eq.'goes08')isat=1
          if(csatid.eq.'goes09')isat=6
          if(csatid.eq.'goes10')isat=3
          if(csatid.eq.'goes11')isat=7
          if(csatid.eq.'goes12')isat=5    

          do j=1,jmax
          do i=1,imax
             laps_vis_norm_natl(i,j)=laps_vis_raw(i,j)
          enddo
          enddo

          print*,csatid,'/',c_sat_type,' vis data already normalized'
          i_dir = -1

c override the namelist value for these data. the fsl-conus data
c have been normalized with l_national = .true.

          print*,'reverse the normalization'

          call normalize_brightness(i4time,lat,lon,laps_vis_norm
     &,imax,jmax,sublat_d,sublon_d,range_m,.true.,iskip_bilin   !note: l_national is hardwired true here.
     &,r_missing_data,6,i_dir,phase_angle_d,specular_ref_angle_d
     &,emission_angle_d,istatus_n)
c
          if(istatus_n .ne. 1) then
             write(*,*)'+++warning+++ bad status returned from
     &normalize laps vis'
          else
             write(*,*)'visible image successfully un-normalized'
          endif
c
c the standard is to save the raw (true) satellite vis counts (svs).
c
          print*,'reverse the stretch and save (raw) counts'
          print*,'stretch ',c_sat_type,' visible'
          print*,'vis-cnt-stretch-in  1/2: ',visin1_g8,visin2_g8
          print*,'vis-cnt-stretch-out 1/2: ',visout1_g8,visout2_g8
          print*
          print*,'apply goes08-to-goes07 stretch and save counts'

          do j=1,jmax
          do i=1,imax
c
c since this data looks like raw goes-7 counts (due to normalization by
c fd, now its), we must reverse the stretch.
c

             if(laps_vis_norm(i,j).ne.r_missing_data)then

               call stretch(visout1_g8,visout2_g8
     &                     ,visin1_g8,visin2_g8
     &                     ,laps_vis_norm(i,j))

               laps_vis_raw(i,j)=laps_vis_norm(i,j)                    !ok, save these. this is svs!

c              call stretch(0., 255., 0., 303.57, laps_vis_norm(i,j))  !reverse the stretch 

c for goes8 - make it look like goes7
               call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c              call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))

c laps_vis_norm is now ready for local normalization
             endif

          enddo
          enddo
 
          write(6,*)
          i_dir = 1
c
c================================================
c wfo switch
c
       elseif(c_sat_type.eq.'wfo')then 

          write(6,*)'these vis data will get normalized - wfo'

          if(csatid.eq.'goes08')then

             isat = 1
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'stretch ',c_sat_type,' visible'
             print*,'in: ',visin1_g8,visin2_g8
             print*,'out: ',visout1_g8,visout2_g8

             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes8 - make it look like goes7
                 call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))

               endif
             enddo
             enddo

          elseif(csatid.eq.'goes10')then

             isat = 3
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'two step process:'
             print*,'1:  stretch ',csatid,'to goes08'
             print*,'2:  stretch goes08 to look like goes07' 

             visin1_cur=vis_cnt_range_in(1,isat)
             visin2_cur=vis_cnt_range_in(2,isat)
             visout1_cur=vis_cnt_range_out(1,isat)
             visout2_cur=vis_cnt_range_out(2,isat)

             print*,csatid,' to goes08 stretch parameters:'
             print*,'   in:  ',visin1_cur,visin2_cur
             print*,'   out: ',visout1_cur,visout2_cur

             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c new stretch for goes10
c j. smart 2-23-99.
c                 call stretch(0., 255., 0., 193.,laps_vis_norm(i,j))
c          3-08-99 commented out and added new stretch parameters for goes10.
c          8-27-99 re-activated the 305 stretch for wfo.
c
                  call stretch(visin1_cur,visin2_cur,
     &                         visout1_cur,visout2_cur,
     &                         laps_vis_norm(i,j))

c                 call stretch(0.,305.,0.,255.,laps_vis_norm(i,j))

cc old:              call stretch(0.,295.,0.,255.,laps_vis_norm(i,j))

c for goes8 - make it look like goes7
                  call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                 call stretch(0., 303.57, 0., 255., laps_vis_norm(i,j))

               endif
             enddo
             enddo

          elseif(csatid.eq.'goes12'.or.csatid.eq.'goes11'
     1      .or. csatid.eq.'goessw'.or.csatid.eq.'goesse'
     1      .or. csatid.eq.'goesea'.or.csatid.eq.'goeswe')then

             isat = 5
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'two step process:'
             print*,'1:  stretch ',csatid,' to goes08'
             print*,'2:  stretch goes08 to look like goes07'

             visin1_cur=vis_cnt_range_in(1,isat)
             visin2_cur=vis_cnt_range_in(2,isat)
             visout1_cur=vis_cnt_range_out(1,isat)
             visout2_cur=vis_cnt_range_out(2,isat)

             print*,csatid, ' to goes08 stretch parameters:'
             print*,'   in:  ',visin1_cur,visin2_cur
             print*,'   out: ',visout1_cur,visout2_cur

             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c new stretch for goes10
c j. smart 2-23-99.
c                 call stretch(0., 255., 0., 193.,laps_vis_norm(i,j))
c          3-08-99 commented out and added new stretch parameters for goes10.
c          8-27-99 re-activated the 305 stretch for wfo.
c
                  call stretch(visin1_cur,visin2_cur,
     &                         visout1_cur,visout2_cur,
     &                         laps_vis_norm(i,j))
c for goes8 - make it look like goes7
                  call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))
               endif
             enddo
             enddo
          endif
          i_dir = 1
c
c=======================================
c gvar switch
c=======================================
       elseif(c_sat_type.eq.'gvr'.or.c_sat_type.eq.'gwc')then

          call make_fnam_lp(i4time,a9time,istatus)
          write(6,*)'gvar visible data: ',csatid,' ',a9time
     1             ,' -------------------------------'

          if(csatid.eq.'goes08')then

             isat = 1
             print*,'stretch ',csatid,' to goes07 look-a-like'
             print*,' in:  ',visin1_g8,visin2_g8
             print*,' out: ',visout1_g8,visout2_g8

             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes8 - make it look like goes7
                   call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                  call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))
                endif
             enddo
             enddo

          elseif(csatid.eq.'goes10')then

             isat = 3
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'two step process:'
             print*,'1:  stretch ',csatid,'to goes08'
             print*,'2:  stretch goes08 to look like goes07'

             visin1_cur=vis_cnt_range_in(1,isat)
             visin2_cur=vis_cnt_range_in(2,isat)
             visout1_cur=vis_cnt_range_out(1,isat)
             visout2_cur=vis_cnt_range_out(2,isat)

             print*,'goes10 to goes08 stretch parameters:'
             print*,'   in:  ',visin1_cur,visin2_cur
             print*,'   out: ',visout1_cur,visout2_cur

             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes9 - make it look like goes8; no goes09. new stretch for goes10
c j. smart 2-23-99.
c                  call stretch(0., 255., 0.,193.,laps_vis_norm(i,j))

c                  call stretch(visin1_cur,visin2_cur
c    &                     ,visout1_cur,visout2_cur
c    &                     ,laps_vis_norm(i,j))

                   call stretch(0.,305.,0.,255.,laps_vis_norm(i,j))

c for goes8 - make it look like goes7
                   call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                  call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))

                endif
             enddo
             enddo

          elseif(csatid.eq.'goes12'.or.csatid.eq.'goes11'
     1      .or. csatid.eq.'goessw'.or.csatid.eq.'goesse'
     1      .or. csatid.eq.'goesea'.or.csatid.eq.'goeswe')then

             isat = 5 
             if(csatid.eq.'goes11')isat=7
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'two step process:'
             print*,'1:  stretch ',csatid,'to goes08'
             print*,'2:  stretch goes08 to look like goes07'

             visin1_cur=vis_cnt_range_in(1,isat)
             visin2_cur=vis_cnt_range_in(2,isat)
             visout1_cur=vis_cnt_range_out(1,isat)
             visout2_cur=vis_cnt_range_out(2,isat)

             print*,csatid,' to goes08 stretch parameters:'
             print*,'   in:  ',visin1_cur,visin2_cur
             print*,'   out: ',visout1_cur,visout2_cur

             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then

                   call stretch(visin1_cur,visin2_cur
     &                     ,visout1_cur,visout2_cur
     &                     ,laps_vis_norm(i,j))

c now make it look like goes7
                   call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                  call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))

                endif
             enddo
             enddo

          endif
          i_dir = 1

       elseif(csatid.eq.'gmssat')then

          isat = 4

          if(c_sat_type.eq.'twn'.or.c_sat_type.eq.'hko')then

             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'currently 1 step process:'
             print*,'1:  stretch ',csatid,' to goes07'

             visin1_cur=vis_cnt_range_in(1,isat)
             visin2_cur=vis_cnt_range_in(2,isat)
             visout1_cur=vis_cnt_range_out(1,isat)
             visout2_cur=vis_cnt_range_out(2,isat)

             print*,'gmssat to goes08 stretch parameters:'
             print*,'   in:  ',visin1_cur,visin2_cur
             print*,'   out: ',visout1_cur,visout2_cur

c         print*,'stretch twn visible: 40.,250.,68.,220.'
c         print*,'stretch twn visible: 38.,320.,68.,220.' !11-12-02
c         print*,'stretch twn visible: 38.,350.,68.,220.' !11-15-02

             do j=1,jmax
             do i=1,imax
                if(laps_vis_norm(i,j).ne.r_missing_data)then
c j. smart 1-25-02.
c   "      2-05-03
                   call stretch(visin1_cur,visin2_cur
     &                     ,visout1_cur,visout2_cur
     &                     ,laps_vis_norm(i,j))

c               call stretch(38.,400.,68.,220.,laps_vis_norm(i,j))
                endif
             enddo
             enddo

          endif

c==========================
c this switch is obsolete
c==========================
       elseif(c_sat_type.eq.'cdf')then

        if(csatid.eq.'goes09'.or.csatid.eq.'goes11')then 

         isat = 6
         if(csatid.eq.'goes11')isat=7

         print*,'stretch ',csatid,' to goes7 look-a-like'
         print*,'two step process:'
         print*,'1:  stretch ',csatid,'to goes08'
         print*,'2:  stretch goes08 to look like goes07'

         visin1_cur=vis_cnt_range_in(1,isat)
         visin2_cur=vis_cnt_range_in(2,isat)
         visout1_cur=vis_cnt_range_out(1,isat)
         visout2_cur=vis_cnt_range_out(2,isat)

         print*,csatid, ' to goes08 stretch parameters:'
         print*,'   in:  ',visin1_cur,visin2_cur
         print*,'   out: ',visout1_cur,visout2_cur

         do j=1,jmax
         do i=1,imax
            if(laps_vis_norm(i,j).ne.r_missing_data)then
               call stretch(visin1_cur,visin2_cur
     &                     ,visout1_cur,visout2_cur
     &                     ,laps_vis_norm(i,j))

c now make it look like goes7
               call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

c                  call stretch(0.,303.57,0.,255.,laps_vis_norm(i,j))

            endif
         enddo
         enddo

        else

         print*,'------------------------------------------'
         print*,'!!warning: unknown sat id for cdf switch!!'
         print*,'!!        sat id = ',csatid,'           !!'
         print*,'!!       data not stretched             !!'
         print*,'------------------------------------------'

        endif

       elseif(csatid.eq.'noaapo')then

          write(6,*)'normalize vis data - noaa polar orbiter'

          if(c_sat_type.eq.'ncp')then 

             isat = 1
             print*,'stretch ',csatid,' to goes7 look-a-like'
             print*,'stretch ',c_sat_type,' visible'
             print*,'in: ',visin1_g8,visin2_g8
             print*,'out: ',visout1_g8,visout2_g8
             print*
             print*,'add 5 min to file time for this data'
             print*,'to approx sat movement to domain'
             print*,'===================================='
             i4time=i4time+300

             do j=1,jmax
             do i=1,imax
               if(laps_vis_norm(i,j).ne.r_missing_data)then
c for goes8 - make it look like goes7
                 call stretch(visin1_g8,visin2_g8
     &                     ,visout1_g8,visout2_g8
     &                     ,laps_vis_norm(i,j))

               endif
             enddo
             enddo

          endif ! c_sat_type

       elseif(trim(csatid) .eq. 'mtsat')then

          write(6,*)'stretch vis data prior to normalization ',csatid

          do j=1,jmax
          do i=1,imax
            if(laps_vis_norm(i,j).ne.r_missing_data)then
c make it look like goes7
              call stretch(0.,1600.,0.,256.,laps_vis_norm(i,j))

            endif
          enddo
          enddo

       elseif(trim(csatid) .eq. 'coms')then

          write(6,*)'stretch vis data prior to normalization ',csatid

          do j=1,jmax
          do i=1,imax
            if(laps_vis_norm(i,j).ne.r_missing_data)then
c make it look like goes7
              call stretch(0.,440.,10.,256.,laps_vis_norm(i,j))

            endif
          enddo
          enddo

       elseif(c_sat_type.eq.'gnp')then

          write(6,*)'stretch vis data prior to normalization ',csatid

          do j=1,jmax
          do i=1,imax
            if(laps_vis_norm(i,j).ne.r_missing_data)then
c make it look like goes7
!             call stretch(0.,4095., 0.,220.,laps_vis_norm(i,j))
!             call stretch(0.,6095.,30.,220.,laps_vis_norm(i,j))
              call stretch(0.,4688.,30.,220.,laps_vis_norm(i,j))

            endif
          enddo
          enddo

       endif ! csatid
c
c ready to normalize vis counts to local domain
c =============================================
       mode_norm = 1
       if(mode_norm .eq. 1)then

          call normalize_brightness(i4time,lat,lon,
     &             laps_vis_norm,imax,jmax,
     &             sublat_d,sublon_d,
     &             range_m,l_national,iskip_bilin,
     &             r_missing_data,6,i_dir,phase_angle_d,
     &             specular_ref_angle_d,emission_angle_d,
     &             istatus_n)
          if(istatus_n .ne. 1) then
             write(*,*)'+++warning+++ bad status returned from
     &normalize laps vis'
             istatus(2) = istatus_n
          else
             print*,'visible image normalized for local domain'
             call check(laps_vis_norm
     &              ,r_missing_data,istatus_n,imax,jmax)
             istatus(2)=istatus_n  
          endif

       elseif(mode_norm .eq. 2)then
!         call refl_to_albedo(reflectance,solalt,land_albedo    ! i
!    1                       ,cloud_albedo)                     ! o
          continue

       else

          print*,'skip normalization: ',csatid

       endif
c =============================================

! get the location of the static grid directory
       call get_directory('static',directory,len_dir)

       var_2d='ldf'
       write(6,*)'read land fraction'

       call find_domain_name(dataroot,domain_name,istatus)

       call rd_laps_static (directory,domain_name,imax,jmax
     1,1,var_2d,units_2d,comment_2d,rland_frac,rspacing_dum,istatus_l)
       if(istatus_l .ne. 1)then
           write(6,*)' error reading laps static-land frac'
           return
       endif

       istatus_a = 1
       call vis_to_albedo(i4time,csatid,
     &                    laps_vis_norm,
     &                    lat,lon,
     &                    imax,jmax,
     &                    r_missing_data,
     &                    phase_angle_d,
     &                    specular_ref_angle_d,
     &                    emission_angle_d,
     &                    rland_frac,
     &                    albedo,
     &                    albedo_min,
     &                    albedo_max,
     &                    n_missing_albedo,
     &                    istatus_a)

       if(istatus_a .ne. 1)then
          write(*,*)'+++warning.+++ err status ',istatus_a,' from
     &               vis_to_albedo'
          istatus(3) = istatus_a
       else
          write(*,*)'albedo successfully computed'
       end if
       call check(albedo,r_missing_data,istatus_a,imax,jmax)
       if(istatus_a .lt. 1) then
          print*,' +++ warning. visible albedo status = ',istatus_a
     1                                                   ,imax*jmax
          istatus(3) = istatus_a
       endif
c
       write(*,*)'successfully remapped vis data'

       call moment(albedo,imax*jmax,
     &             ave,adev,sdev,var,skew,curt,
     &             istatus_m)
       if(istatus_m.ne.0)then
          print*,'error returned from subroutine moment'
          goto 999
       endif

       icnt=0
       do j=1,jmax
       do i=1,imax
          if(albedo(i,j).ne.r_missing_data)then
           if(albedo(i,j)-ave.gt.3.0*sdev)then
              icnt=icnt+1
              isave=i
              jsave=j
           endif
           if(albedo(i,j).eq.albedo_max)then
              ismax=i
              jsmax=j
           endif
           if(albedo(i,j).eq.albedo_min)then
              ismin=i
              jsmin=j
           endif
          endif
       enddo
       enddo

       i=isave
       j=jsave
       print*,'  ================='
       print*,'  albedo statistics'
       print*,'  ================='
       if(icnt.gt.0.and.ismax.gt.0.and.jsmax.gt.0.and.
     &ismin.gt.0.and.jsmin.gt.0)then
          print*,'  n > 3 standard dev = ',icnt
          print*,'  last i/j > 3 stand dev = ',i,j
          print*,'  laps_vis_norm(i,j)= ',laps_vis_norm(i,j)
          print*,'  average albedo = ', ave
          print*,'  standard dev = ',sdev
          print*,'  i/j/count of max albedo ',ismax,jsmax,
     &laps_vis_norm(ismax,jsmax)
          print*,'  i/j/count of min albedo ',ismin,jsmin,
     &laps_vis_norm(ismin,jsmin)
       else
          print*,'all albedo < 3 standard dev of ave'
          print*,'no max/min output in satvis_process'
       endif

       goto 999
898    write(6,*)'error getting r_msng_sat_flag'
c     
 999   return
       end
