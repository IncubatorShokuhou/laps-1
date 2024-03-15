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

        subroutine get_heights_hydrostatic
     1              (temp_3d_k          ! input (3d temp k)
     1              ,pres_sfc_pa        ! input (surface pressure pa)
     1              ,pres_3d            ! input (3d pressure pa)
     1              ,sh_3d              ! input (3d specific humidity)
     1              ,topo               ! input (terrain m)
     1              ,ni,nj,nk           ! input (dimensions)
     1              ,heights_3d         ! output (m)
     1              ,istatus)           ! output

!           1991        steve albers
!       nov 1991        steve albers    change to sfc p input + terrain
!           1996        steve albers    expanded 'esat_lut' array from -100c
!                                       to the colder value of -120
!
!       for each horizontal gridpoint, the 'pres_sfc_pa' at the lowest level 
!       must be greater than or equal to the reference (surface) pressure.

        real esat_lut(-120:+100)

        real temp_3d_k(ni,nj,nk)
        real pres_sfc_pa(ni,nj)
        real topo(ni,nj)
        real heights_3d(ni,nj,nk)
        real sh_3d(ni,nj,nk)

        real heights_below(ni,nj)
        real a_below(ni,nj)
        real pres_sfc_mb(ni,nj)
        real z_correction(ni,nj)

        real pres_3d(ni,nj,nk)                    

        real pres_3d_mb(ni,nj,nk)               ! local
        real p_1d_mb(nk)                        ! local
        real alog_array(ni,nj,nk)               ! local

        real make_td, k_to_c

        include 'constants.inc'
        real c1
        parameter (c1 = ep_1) 

        real c2
        parameter (c2 = r_d / (2. * grav) )

cdoc    generate heights of grid points using laps data and hydrostatic equation
        write(6,*)' generating 3d height grid'
        write(6,*)' integrating hydrostatic equation over the laps grid'
c       write(6,*)' initialize and calculate first level'

        do it = -120,+100
            esat_lut(it) = esat(float(it))
        enddo

        call get_r_missing_data(r_missing_data,istatus)
        z_correction = r_missing_data

        do k = 1,nk
            do j = 1,nj
            do i = 1,ni
!               p_1d_mb(k) = zcoord_of_level(k)/100.
                pres_3d_mb(i,j,k) = pres_3d(i,j,k) / 100.
                if(k .gt. 1)then
                    alog_array(i,j,k-1) = 
     1                  alog(pres_3d_mb(i,j,k-1)/pres_3d_mb(i,j,k))
                endif
            enddo ! i
            enddo ! j
        enddo


        do j = 1,nj
        do i = 1,ni
            pres_sfc_mb(i,j) = pres_sfc_pa(i,j) / 100.
            heights_below(i,j) = 0.
            t_2d_c_below = k_to_c(temp_3d_k(i,j,1))
            sh_below = sh_3d(i,j,1)
            w_below = sh_below / (1. - sh_below)
            a_below(i,j)= temp_3d_k(i,j,1) * (1. + c1 * w_below)
        enddo ! i
        enddo ! j

        do k = 2,nk ! integrate upward one level at a time.

ccccccccccccccccccc                         istat = lib$show_timer(,,,)

            do j = 1,nj
            do i = 1,ni
                sh_value = sh_3d(i,j,k)
                w_value = sh_value / (1. - sh_value)

                a1= temp_3d_k(i,j,k) * (1. + c1 * w_value)

                z_add = c2 * (a1+a_below(i,j)) * (alog_array(i,j,k-1))       

                heights_3d(i,j,k) = heights_below(i,j) + z_add

!               test whether this layer contains the surface
                if(pres_3d_mb(i,j,k  ) .le. pres_sfc_mb(i,j) .and.
     1             pres_3d_mb(i,j,k-1) .ge. pres_sfc_mb(i,j)      )then       

                    alog_factor = alog( pres_3d_mb(i,j,k-1)
     1                                 /pres_sfc_mb(i,j)   )
                    z_to_sfc = z_add * (  alog_factor 
     1                                  / alog_array(i,j,k-1) )

                    z_below = topo(i,j) - z_to_sfc

!                   z_below = topo(i,j) -
!       1           z_thk(pres_sfc_mb(i,j),pres_3d_mb(i,j,k-1),t_1d_c,td_1d_c
!       1                                               ,alog_array,esat_lut,2)

                    z_correction(i,j) = z_below - heights_below(i,j)
                endif

                heights_below(i,j) = heights_3d(i,j,k)
                a_below(i,j) = a1

            enddo ! i
            enddo ! j

c           write(6,*)' completed level',k

        enddo ! k

ccccccccccccccccccc                         istat = lib$show_timer(,,,)

        ndiffp = 0
        diffp_max = -abs(r_missing_data)
        diffp_min = +abs(r_missing_data)
        iwrite = 0

      ! apply constant of integration to the heights
        do j = 1,nj
        do i = 1,ni
            z_correction_orig = z_correction(i,j)
            if(z_correction(i,j) .eq. r_missing_data)then
                z_correction(i,j) = 0. ! level 1 becomes the reference height
                ndiffp = ndiffp + 1
            endif

            heights_3d(i,j,1) = z_correction(i,j)

!           obtain diffp information
            diffp = pres_sfc_mb(i,j) - pres_3d_mb(i,j,1) 
            diffp_max = max(diffp,diffp_max)
            diffp_min = min(diffp,diffp_min)
            iwrite = iwrite + 1
            if(iwrite .le. 10)then
                write(6,*)i,j,pres_sfc_mb(i,j),pres_3d_mb(i,j,1)
     1               ,diffp,diffp_max,diffp_min,ndiffp,z_correction_orig       
            endif
        enddo ! i
        enddo ! j

        do k = 2,nk
            do j = 1,nj
            do i = 1,ni
                heights_3d(i,j,k) = heights_3d(i,j,k) + z_correction(i,j
     1)
            enddo ! j
            enddo ! i
        enddo ! k

c       write(6,*)' added in constant of integration to the heights'
ccccccccccccccccccc                         istat = lib$show_timer(,,,)

        write(6,*)' range of diffp ( psfc - p[1] ) = ',diffp_min
     1                                                ,diffp_max,'mb'

!       return with diagnostic information
        ngrid = ni*nj
        if(ndiffp .eq. ngrid)then
            write(6,*)' warning: surface extends outside '
     1                ,'3d grid at all grid points'    
            write(6,*)' output heights are relative to level 1'
            istatus = -1
            return
        elseif(ndiffp .gt. 0)then
            write(6,*)' error: surface extends outside '
     1                ,'3d grid at some grid points',ndiffp    
            istatus = 0
            return
        else
            write(6,*)' success in get_heights_hydrostatic'
            istatus = 1
            return
        endif

        end


        function z_thk(pt,p,t,td,alog_array,esat_lut,n)
c
cdoc    this function returns the thickness of a layer bounded by pressure
cdoc    p(1) at the bottom and pressure pt at the top.
c
c       baker,schlatter 17-may-1982     original version
c       albers                 1990     restructured
c       albers                 1996     expanded 'esat_lut' array from
c                                       -100c to the colder value of -120
c
c   on input:
c       p = pressure (mb).  note that p(i).gt.p(i+1).
c       t = temperature (celsius)
c       td = dew point (celsius)
c       n = number of levels in the sounding and the dimension of
c           p, t and td
c   on output:
c       z = geometric thickness of the layer (m)
c
c   the algorithm involves numerical integration of the hydrostatic
c   equation from p(1) to pt. it is described on p.15 of stipanuk
c   (1973).
c
        dimension t(n),p(n),td(n),tk(100),alog_array(n)
c       c1 = .001*(1./eps-1.) where eps = .62197 is the ratio of the
c                             molecular weight of water to that of
c                             dry air. the factor 1000. converts the
c                             mixing ratio w from g/kg to a dimension-
c                             less ratio.
c       c2 = r/(2.*g) where r is the gas constant for dry air
c                     (287 kg/joule/deg k) and g is the acceleration
c                     due to the earth's gravity (9.8 m/s**2). the
c                     factor of 2 is used in averaging two virtual
c                     temperatures.
        real esat_lut(-120:+100)

        include 'constants.inc' 
        real c1,c2
        parameter (c1 = .001 * ep_1)
        parameter (c2 = r_d / (2. * grav) )

        do 5 i= 1,n
           tk(i)= c_to_k(t(i))
    5   continue

        do i = 1,n-1
            j = i + 1

            if (p(j) .gt. pt)then             ! we are still in complete layers
                a1= tk(j)*(c2 + w_laps(td(j),p(j),esat_lut))
                a2= tk(i)*(c2 + w_laps(td(i),p(i),esat_lut))
                z_thk = z_thk + (a1+a2)*(alog_array(i))

            elseif(p(j) .eq. pt)then          ! finish up on complete layer
                a1= tk(j)*(c2 + w_laps(td(j),p(j),esat_lut))
                a2= tk(i)*(c2 + w_laps(td(i),p(i),esat_lut))
                z_thk = z_thk + (a1+a2)*(alog_array(i))
                return

            else ! p(j) .lt. pt               ! finish up on partial layer
                a1= tk(j)*(c2 + w_laps(td(j),p(j),esat_lut))
                a2= tk(i)*(c2 + w_laps(td(i),p(i),esat_lut))
                z_thk = z_thk + (a1+a2)*(alog(p(i)/pt))
                return

            endif

        enddo ! i


        write(6,*)' error, pt out of bounds in z_thk'
        z_thk = -1.0

        return
        end


        function w_laps(t,p,esat_lut)

cdoc    convert t(c) and p to w. this is a fast approximate routine.
cdoc    output w_laps is dimensionless mixing ratio

!       this function really only works when t > -50c but the efficiency will
!       outweigh the error when t < -50c in this application

        real esat_lut(-120:+100)

        include 'constants.inc'
        real c1,c2,const
        parameter (c1 = .001 * ep_1)
        parameter (c2 = r_d / (2. * grav) )

        parameter (const = 622. * c1 * c2)

c       es(x)=6.1078+x*(.443652+x*(.014289+x*(2.65065e-4+x*
c    1 (3.03124e-6+x*(2.034081e-8+x*(6.13682e-11))))))
c
c
        x= esat_lut(nint(t))
        w_laps= const*x/(p-x)
        return
        end


