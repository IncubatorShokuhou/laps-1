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

        subroutine helicity_laps(uanl,vanl,ustorm,vstorm
     1                          ,heights_3d,topo
     1                          ,imax,jmax,kmax,helicity,istatus)   

cdoc    calculate storm relative helicity over a 2-d grid    

        real ustorm(imax,jmax),vstorm(imax,jmax)
        real usfc(imax,jmax),vsfc(imax,jmax)
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)
        real helicity(imax,jmax)
        real topo(imax,jmax)

!       1998 steve albers - overhauled

        icount_write = 0

        write(6,*)' computing helicity for 0-3 km agl'

        do j = 1,jmax
        do i = 1,imax

            area_sum = 0.

!           layer is 0-3 km agl, denoted from "sfc" to "top"
            rksfc = height_to_zcoord2(topo(i,j)      ,heights_3d
     1                               ,imax,jmax,kmax,i,j,istatus)
            if(istatus .ne. 1)return

            rktop = height_to_zcoord2(topo(i,j)+3000.,heights_3d
     1                               ,imax,jmax,kmax,i,j,istatus)
            if(istatus .ne. 1)return

!           get storm relative wind at the sfc, using interpolated sfc wind
            ksfc  = int(rksfc)
            frack = rksfc - ksfc
            u_sfc = uanl(i,j,ksfc)   * (1.-frack) 
     1            + uanl(i,j,ksfc+1) * frack       
            v_sfc = vanl(i,j,ksfc)   * (1.-frack) 
     1            + vanl(i,j,ksfc+1) * frack

            u_rel_l = u_sfc - ustorm(i,j)
            v_rel_l = v_sfc - vstorm(i,j)

            klow  = int(rksfc) + 1          ! 1st level above the sfc
            khigh = int(rktop) + 1          ! 1st level above top of layer

            if(i .eq. i/10*10 .and. j .eq. 42)then
                iwrite = 1
                write(6,*)
                write(6,*)'i/j = ',i,j
     1                   ,'u/vstorm = ',ustorm(i,j),vstorm(i,j)
            else
                iwrite = 0
            endif

            do k = klow,khigh

                if(k .lt. khigh)then
!                   get storm relative wind at this level
                    u_rel_u = uanl(i,j,k) - ustorm(i,j)
                    v_rel_u = vanl(i,j,k) - vstorm(i,j)

                else ! k = khigh, use 3km agl values instead of this laps level
                    frack = rktop - int(rktop)
                    u_anl = uanl(i,j,k-1) * (1.-frack) 
     1                    + uanl(i,j,k)   * frack
                    v_anl = vanl(i,j,k-1) * (1.-frack) 
     1                    + vanl(i,j,k)   * frack

                    u_rel_u = u_anl - ustorm(i,j)
                    v_rel_u = v_anl - vstorm(i,j)

                endif


!               cross product of wind vectors on top and bottom of layer
                xprod = u_rel_l * v_rel_u - v_rel_l * u_rel_u       

!               incremental area of hodograph
                area = .5 * xprod

!               total area of hodograph
                area_sum = area_sum + area

                if(iwrite .eq. 1)then
                    write(6,51)k,u_rel_l,v_rel_l,u_rel_u,v_rel_u
     1                        ,area,area_sum
 51                 format(i6,6f11.2)
                endif

                u_rel_l = u_rel_u
                v_rel_l = v_rel_u

            enddo ! k

            helicity(i,j) = area_sum * (-2.) 

            if(iwrite .eq. 1)then
                write(6,101)i,j,klow,khigh,helicity(i,j),ustorm(i,j)
     1                    ,(uanl(i,j,k),k=1,min(kmax,21))        
101             format(4i4,f10.5,f8.3/7e11.3/7e11.3/7e11.3)
            endif

        enddo ! j
        enddo ! i

        istatus = 1

        return
        end
