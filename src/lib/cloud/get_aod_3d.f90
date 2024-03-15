
       subroutine get_aod_3d(pres_3d,heights_3d,topo_2d,ni,nj,nk &
                            ,aod,aod_ref,i_aero_synplume,i_aero_1d,aod_3d)

       use mem_namelist, only: redp_lvl,aero_scaleht,grid_spacing_m,aod_ha,ht_ha,alpha_ha
       use mem_allsky, only: mode_aero_cld,nc

       include 'rad.inc'

       real aod_3d(ni,nj,nk) ! aerosol extinction coefficient
       real pres_3d(ni,nj,nk)
       real heights_3d(ni,nj,nk)
       real topo_2d(ni,nj)

       scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0-1
       fraclim(x,x1,x2) = min(max((x-x1)/(x2-x1),0.),1.)
       scurvel(x,x1,x2) = scurve(fraclim(x,x1,x2)) ! scurve over bounds

       gas_scale_height = 8000.

       write(6,*)'subroutine get_aod_3d: aod/redp_lvl/aero_scaleht = ' &
                                        ,aod,redp_lvl,aero_scaleht

       write(6,11)mode_aero_cld,i_aero_1d,i_aero_synplume
11     format('  mode_aero_cld/i_aero_1d/i_aero_synplume = ',3i5)

       write(6,*)'aod_ha / alpha_ha new method ',aod_ha,alpha_ha

       alpha_ha = aod_ha / ((h2_ha-h1_ha)+0.5*(h3_ha-h2_ha))
       write(6,*)'aod_ha / alpha_ha old method ',aod_ha,alpha_ha

!      transitional assignments if needed
!      h1_ha = ht_ha(2)
!      h2_ha = ht_ha(3)

       if(i_aero_1d .eq. 1)then
         write(6,*)' set aerosols from 1d parameters'
         do k = 1,nk
           do i = 1,ni
           do j = 1,nj
!            sum_aod = 0.
!            pratio = pres_3d(i,j,k) / pref
             htarg = heights_3d(i,j,k)
             h_agl = htarg - redp_lvl 
             alpha_low = (aod/aero_scaleht) * exp(-h_agl/aero_scaleht)
             if(htarg .gt. h1_ha .and. & 
                htarg .le. h2_ha .and. aod_ha .gt. 0.)then
                 alpha_high = alpha_ha
             else
                 alpha_high = 0.
             endif
             aod_3d(i,j,k) = alpha_low + alpha_high
           enddo ! j
           enddo ! i
         enddo ! k

         i = ni/2
         j = nj/2  
         sum_aod = 0.
         do k = 1,nk
             h_agl = heights_3d(i,j,k) - redp_lvl 
             if(h_agl .gt. 0.)then
                 ave_aod = 0.5 * (aod_3d(i,j,k)+aod_3d(i,j,k-1))
                 sum_aod = sum_aod + ave_aod       & 
                         * (heights_3d(i,j,k) - heights_3d(i,j,k-1))
             endif
             write(6,101)k,h_agl,aod_3d(i,j,k),sum_aod
101          format('k,agl,aod_3d,sum',i5,f9.1,e15.8,f9.5)
         enddo ! k
         i4_elapsed = ishow_timer()
       else
         write(6,*)' skip aerosols from 1d parameters'
       endif ! mode_aero_cld

       if(i_aero_synplume .eq. 1)then ! synthetic aerosol plume
         ext_syn = 0.001
         iplume = ni/2 + 1
         jplume = nj/2 + 1
         iwp = 1
         write(6,*)' adding synthetic aerosol plume of',ext_syn,iplume,jplume
         aod_3d(iplume-iwp:iplume+iwp,jplume-iwp:jplume+iwp,1:nk-8) = ext_syn
       elseif(i_aero_synplume .eq. 20)then ! synthetic aerosol gradient (co)
         ri_full = float(ni/2) + 0000. / grid_spacing_m
         ri_none = float(ni/2) + 4000. / grid_spacing_m
         write(6,*)' adding synthetic aerosol i gradient ',ni/2,ri_none,ri_full
         do i = 1,ni
           aero_scale = max(scurvel(float(i),ri_none,ri_full),.0001)
           do k = 1,nk
           do j = 1,nj
               if(heights_3d(i,j,k) .gt. 2000.)then
                   aod_3d(i,j,k) = aod_3d(i,j,k) * 0.15
               else
                   aod_3d(i,j,k) = aod_3d(i,j,k) * (0.15 + 0.85 * aero_scale)
               endif
           enddo ! j
           enddo ! k
           if(i .eq. (i/10) * 10)then
              write(6,*)' i,aero_scale ',i,aero_scale,aod_3d(i,nj/2,:)
           endif
         enddo ! i
         aod = 0.
         aod_ref = 0.
       elseif(i_aero_synplume .eq. 2)then ! synthetic aerosol gradient (gibraltar)
         rj_full = float(nj) - 400000. / grid_spacing_m
         rj_none = float(nj) - 300000. / grid_spacing_m
         write(6,*)' adding synthetic aerosol j gradient ',ni/2,rj_none,rj_full
         do j = 1,nj
           aero_scale = max(scurvel(float(j),rj_none,rj_full),.0001)
           do k = 1,nk
           do i = 1,ni
               if(heights_3d(i,j,k) .gt. 2000.)then
                   aod_3d(i,j,k) = aod_3d(i,j,k) * 0.15
               else
                   aod_3d(i,j,k) = aod_3d(i,j,k) * (0.15 + 0.85 * aero_scale)
               endif
           enddo ! j
           enddo ! k
           if(j .eq. (j/50) * 50)then
              write(6,*)' j,aero_scale ',j,aero_scale,aod_3d(ni/2,j,:)
           endif
         enddo ! j
         aod = 0.
         aod_ref = 0.
       else
         write(6,*)' skip synthetic aerosol plume'
       endif

       return
       end
