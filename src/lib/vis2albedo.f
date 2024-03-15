

        subroutine vis_to_albedo(i4time,csatid,
     &                           r_norm_vis_cnts_in,
     &                           lat,lon,
     &                           imax,jmax,
     &                           r_missing_data,
     &                           phase_angle_d,
     &                           specular_ref_angle_d,
     &                           emission_angle_d,
     &                           rland_frac,
     &                           albedo_out,
     &                           albedo_min,
     &                           albedo_max,
     &                           n_missing_albedo,
     &                           istatus)
c
c program computes albedo given visible image array normalized for brightness.
c
c     j. smart       17-mar-1994           original version taken from vis_to_
c                                          albedo written by s. albers.
c     s. albers       2-mar-1995           set threshold solar alt to 15.
c     s. albers      22-nov-1995           extra phase angle constraint
c     s. albers      22-aug-1996           extra specular ref angle constraint
c     j. smart        8-apr-1998           added csatid to argument list, added to
c                                          formatted output.
c       "             6-jul-1999           added emission_angle_d to argument list
c
c***parameter and variables list
c
        use mem_namelist, only: solalt_thr_vis

        implicit none

        real          cld_cnts,cld_albedo,frac,term1,term2
        real          cloud_frac_vis
        real          albedo_to_cloudfrac,cloudfrac_to_albedo
        real          rlnd_cnts,rlnd_albedo
        parameter       (cld_cnts=220.,
     &                   rlnd_cnts=68.,
     &                   cld_albedo=0.85,
     &                   rlnd_albedo=0.15)

        integer         imax,jmax
        integer         i4time

        real          r_norm_vis_cnts_in(imax,jmax)
        real          lat(imax,jmax)
        real          lon(imax,jmax)
        real          phase_angle_d(imax,jmax)
        real          specular_ref_angle_d(imax,jmax)
        real          rland_frac(imax,jmax)
        real          solar_alt_d
        real          albedo
        real          albedo_out(imax,jmax)
        real          albedo_min,albedo_max
        real          r_norm_vis_cnts_mn,r_norm_vis_cnts_mx 
        real          r_missing_data
        real          jline, iline, jdiff, idiff
        real          emission_angle_d(imax,jmax)
        integer         istatus, n_missing_albedo
        integer         i,j

        real arg
        character*(*)   csatid
c
c     ------------------------- begin ---------------------------------

        write(6,*)' solalt_thr_vis (from namelist) =',solalt_thr_vis
!       solalt_thr_vis = 15.
!       write(6,*)' solalt_thr_vis (from hardwire) =',solalt_thr_vis

        albedo_min=1.0
        albedo_max=0.0

        r_norm_vis_cnts_mx = 0.
        r_norm_vis_cnts_mn = 256.

        write(6,*)' subroutine vis2albedo:'
c
c       write(6,28)
c28      format(1x,' i   j   n vis cnts   solalt deg',/,40('-'))
        do j = 1,jmax
           jline = float(j)/10.
           jdiff = jline - int(jline)
           do i = 1,imax
              iline = float(i)/10.
              idiff = iline - int(iline)

              if(r_norm_vis_cnts_in(i,j).ne.r_missing_data)then

                 call solalt(lat(i,j),lon(i,j),i4time,solar_alt_d)

c             if(idiff.eq.0.00 .and. jdiff.eq.0.00)then
c                write(6,29)i,j,r_norm_vis_cnts_in(i,j),solar_alt_d
c29               format(1x,2i3,2x,2f8.2)
c             end if

!                test for favorable geometry
                 if(      solar_alt_d .gt. solalt_thr_vis
     1                            .and.
     1          (solar_alt_d .gt. 23. .or. phase_angle_d(i,j) .gt. 20.)
     1                            .and.
     1          (rland_frac(i,j) .gt. 0.5 
     1                         .or. specular_ref_angle_d(i,j) .gt. 10.)
     1                            .and.
     1                    emission_angle_d(i,j) .gt. 15.       
     1                                                            )then       

                   arg = (r_norm_vis_cnts_in(i,j)- rlnd_cnts) /
     &                 (cld_cnts - rlnd_cnts)
         
                   albedo = rlnd_albedo + arg *
     &                   (cld_albedo - rlnd_albedo)
                   albedo_out(i,j)=min(max(albedo,-0.5),+1.5) ! reasonable

                   if(solar_alt_d .lt. 20.)then ! enabled for now
!                    fudge the albedo at low solar elevation angles < 20 deg
                     frac = (20. - solar_alt_d) / 10.
                     term1 = .13 * frac
                     term2 = 1. + term1

                     cloud_frac_vis = 
     1                           albedo_to_cloudfrac(albedo_out(i,j))
                     cloud_frac_vis = (cloud_frac_vis + term1) * term2
                     albedo_out(i,j) = 
     1                           cloudfrac_to_albedo(cloud_frac_vis)
                   endif

!                  additional stretch
!                  call stretch2(0.0,1.0,.09,1.0 ,albedo_out(i,j))
                   call stretch2(0.0,1.0,.04,1.15,albedo_out(i,j))      
c                                                               excesses
c accumulate extrema
c
                   r_norm_vis_cnts_mn = min(r_norm_vis_cnts_in(i,j)
     1                                     ,r_norm_vis_cnts_mn)
                   r_norm_vis_cnts_mx = max(r_norm_vis_cnts_in(i,j)
     1                                     ,r_norm_vis_cnts_mx)
                   albedo_min = min(albedo_out(i,j),albedo_min)
                   albedo_max = max(albedo_out(i,j),albedo_max)
   
                 else              ! albedo .eq. missing_data

                   albedo_out(i,j) = r_missing_data
                   n_missing_albedo = n_missing_albedo + 1

                 endif             ! qc based on geometry

              else

                 albedo_out(i,j) = r_missing_data
                 n_missing_albedo = n_missing_albedo + 1

              endif                ! r_norm_vis_cnts_in = r_missing_data

           end do
         end do

         write(6,*)' n_missing_albedo = ',n_missing_albedo
         write(6,105)csatid,r_norm_vis_cnts_mn,r_norm_vis_cnts_mx
         write(6,106)csatid,albedo_min,albedo_max

 105     format(1x,a6,'  normalized counts range: ',2f10.2)
 106     format(1x,a6,'  albedo range:            ',2f10.2)

        return
        end

        function albedo_to_cloudfrac(albedo)

!       we might refer to a new version of this to return several quantities:
!       cloud fraction/opacity 
!       cloud albedo (assuming a black terrain surface), related to back scatter
!       cloud optical depth

!       such a routine, 'albedo_to_clouds', is in 
!           'src/lib/modules/module_cloud_rad.f90'

        clear_albedo = .2097063
        cloud_albedo = .4485300
!       cloud_albedo = .40

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo,0.,1.,arg)

        albedo_to_cloudfrac = arg

        return
        end

        function cloudfrac_to_albedo(cloud_frac_vis)

        clear_albedo = .2097063
        cloud_albedo = .4485300
!       cloud_albedo = .40

        arg = cloud_frac_vis

        call stretch2(0.,1.,clear_albedo,cloud_albedo,arg)

!       cloudfrac_to_albedo = (cloud_frac_vis + .87808) / 4.18719

        cloudfrac_to_albedo = arg

        return
        end

c-------------------------------------------------------------------------------
        subroutine stretch2(il,ih,jl,jh,rarg)

        implicit        none

        real          a,b,il,ih,jl,jh,rarg

        a = (jh - jl) / (ih - il)
        b =  jl - il * a

        rarg = a * rarg + b

        return
        end

