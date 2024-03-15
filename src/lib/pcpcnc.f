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


        subroutine cpt_pcp_cnc(ref_3d,temp_3d
     1                                  ,rh_3d_pct    ! input
     1                                  ,cldpcp_type_3d  ! input
     1                                  ,ni,nj,nk     ! input
     1                                  ,c_z2m        ! input
     1                                  ,pres_3d      ! input
     1                                  ,pcp_cnc_3d   ! output (kg/m^3)
     1                                  ,rai_cnc_3d   ! output (rain)
     1                                  ,sno_cnc_3d   ! output (snow)
     1                                  ,pic_cnc_3d)  ! output (precip ice)


        include 'constants.inc'

        real temp_3d(ni,nj,nk)
        real rh_3d_pct(ni,nj,nk)
        real ref_3d(ni,nj,nk)
        integer cldpcp_type_3d(ni,nj,nk)

        character*20 c_z2m

        real pres_3d(ni,nj,nk)    ! pa
        real pcp_cnc_3d(ni,nj,nk) ! kg/m^3
        real rai_cnc_3d(ni,nj,nk) ! kg/m^3
        real sno_cnc_3d(ni,nj,nk) ! kg/m^3
        real pic_cnc_3d(ni,nj,nk) ! kg/m^3

        real rate_3d(ni,nj,nk)

        npcp = 0

!       get 3d precip rates one layer at a time with zr routine
        do k = 1,nk
            call zr(ref_3d(1,1,k),ni,nj,rate_3d(1,1,k))
        enddo ! k

        write(6,*)' c_z2m = ',c_z2m

        if(c_z2m .eq. 'albers')then
          do j = 1,nj
          do i = 1,ni
          do k = 1,nk
            pressure = pressure_of_level(k)

            ipcp_type = cldpcp_type_3d(i,j,k) / 16   ! pull out precip type

            if(ipcp_type .ne. 0)then
!                call cpt_fall_velocity(ipcp_type,pressure,temp_3d(i,j,k)
!     1                                                  ,fall_velocity)
! adan add
                call cpt_fall_velocity(ipcp_type,pressure,temp_3d(i,j,k)
     1                                 ,ref_3d(i,j,k),fall_velocity)
                call cpt_concentration(rate_3d(i,j,k),fall_velocity
     1                                          ,pcp_cnc_3d(i,j,k))

                npcp = npcp + 1

                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    rai_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    sno_cnc_3d(i,j,k) = 0.
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 2)then                       ! snow
                    rai_cnc_3d(i,j,k) = 0.
                    sno_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! ip or hail
                    rai_cnc_3d(i,j,k) = 0.
                    sno_cnc_3d(i,j,k) = 0.
                    pic_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)

                endif ! ipcp_type

            else  ! ipcp_type = 0
                pcp_cnc_3d(i,j,k) = 0.
                rai_cnc_3d(i,j,k) = 0.
                sno_cnc_3d(i,j,k) = 0.
                pic_cnc_3d(i,j,k) = 0.

            endif ! ipcp_type .ne. 0

          enddo ! k
          enddo ! i
          enddo ! j

        elseif(c_z2m .eq. 'rams')then                           
          do j = 1,nj
          do i = 1,ni
          do k = 1,nk

            ipcp_type = cldpcp_type_3d(i,j,k) / 16   ! pull out precip type

!           compute the basic reflectivity 
!           w=264083.11*(rainmr(i,j,k)  &
!               +0.2*(icemr(i,j,k)+snowmr(i,j,k))  &
!               +2.0*graupelmr(i,j,k))
!           w=max(1.,w)
!           refl(i,j,k)=17.8*alog10(w)

!           consider mixed graupel/rain/snow between +5 and -15c?

            if(ipcp_type .ne. 0)then
                npcp = npcp + 1
                rho = pres_3d(i,j,k) / (r_d * temp_3d(i,j,k))      ! kg_m^3
                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    rainfrac = 1.0
                    snowfrac = 0.0
                    icefrac = 0.0

                    w = 10.**(ref_3d(i,j,k)/17.8)
                    rainmr = w / 264083.11                         
                    pcp_cnc_3d(i,j,k) = rainmr * rho               ! kg_m^3
                    rai_cnc_3d(i,j,k) = rainmr * rho               ! kg_m^3

                elseif(ipcp_type .eq. 2)then                       ! snow (or mixed)

!                   calculate wet bulb temperature
                    zero_c = 273.15
                    t_c         = temp_3d(i,j,k) - zero_c
                    td_c        = dwpt(t_c,rh_3d_pct(i,j,k))
                    pressure_mb = pres_3d(i,j,k) / 100.
                    t_wb_c = twet_fast(t_c,td_c,pressure_mb)

                    rainfrac = max(min((t_wb_c / 1.3),1.),0.)
                    snowfrac = 1.0 - rainfrac

                    rainfrac = 0.0
                    snowfrac = 1.0
                    icefrac = 0.0

                    w = 10.**(ref_3d(i,j,k)/17.8)
                    snowmr = w / (264083.11*0.2)                   
                    pcp_cnc_3d(i,j,k) = snowmr * rho               ! kg_m^3
                    sno_cnc_3d(i,j,k) = snowmr * rho               ! kg_m^3

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! ip or hail
                    rainfrac = 0.0
                    snowfrac = 0.0
                    icefrac = 1.0

                    w = 10.**(ref_3d(i,j,k)/17.8)
                    graupelmr = w / (264083.11*2.0)                
                    pcp_cnc_3d(i,j,k) = graupelmr * rho            ! kg_m^3
                    pic_cnc_3d(i,j,k) = graupelmr * rho            ! kg_m^3

                endif

                refarg = 1.0*rainfrac + 0.2*snowfrac + 2.0*icefrac
                w = 10.**(ref_3d(i,j,k)/17.8)
                totalmr = w / (264083.11*refarg)

                rainmr    = totalmr * rainfrac
                snowmr    = totalmr * snowfrac
                graupelmr = totalmr * icefrac

                pcp_cnc_3d(i,j,k) = totalmr * rho              ! kg_m^3
                rai_cnc_3d(i,j,k) = rainmr * rho               ! kg_m^3
                sno_cnc_3d(i,j,k) = snowmr * rho               ! kg_m^3
                pic_cnc_3d(i,j,k) = graupelmr * rho            ! kg_m^3

            else  ! ipcp_type = 0
                pcp_cnc_3d(i,j,k) = 0.
                rai_cnc_3d(i,j,k) = 0.
                sno_cnc_3d(i,j,k) = 0.
                pic_cnc_3d(i,j,k) = 0.

            endif

          enddo ! k
          enddo ! i
          enddo ! j

        elseif(c_z2m .eq. 'kessler')then 
          do j = 1,nj
          do i = 1,ni
          do k = 1,nk

            ipcp_type = cldpcp_type_3d(i,j,k) / 16   ! pull out precip type

!           compute the basic reflectivity 
!           refl(i,j,k) =17300.0 * &
!                 (rho(i,j,k) * 1000.0 * &
!                  max(0.0,rainmr(i,j,k)))**svnfrth

!           add the ice component
!           refl(i,j,k)=refl(i,j,k) + &
!                 38000.0*(rho(i,j,k) * 1000.0 * &
!                 max(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

            if(ipcp_type .ne. 0)then
                npcp = npcp + 1
                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    a = 17300.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/a)**(4./7.)                         ! g/m^3
                    pcp_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3
                    rai_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3

                elseif(ipcp_type .eq. 2)then                       ! snow
                    b = 38000.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/b)**(1./2.2)                        ! g/m^3
                    pcp_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3
                    sno_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! ip or hail
                    b = 38000.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/b)**(1./2.2)                        ! g/m^3
                    pcp_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3
                    pic_cnc_3d(i,j,k) = rho_q / 1000.              ! kg_m^3

                endif

            else  ! ipcp_type = 0
                pcp_cnc_3d(i,j,k) = 0.
                rai_cnc_3d(i,j,k) = 0.
                sno_cnc_3d(i,j,k) = 0.
                pic_cnc_3d(i,j,k) = 0.

            endif

          enddo ! k
          enddo ! i
          enddo ! j

        else
          write(6,*)' error - unknown value of c_z2m in pcpcnc ',c_z2m

        endif

        write(6,*)' npcp/max = ',npcp,maxval(pcp_cnc_3d)

        return
        end

        subroutine cpt_fall_velocity(ipcp_type,p,t,dbz,fall_velocity)

        dbz_eff = max(dbz,1.0)
        vvmax = 4.32*dbz_eff**0.0714286      ! adan add

        if(ipcp_type .eq. 1)then ! rain
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! adan add
        elseif(ipcp_type .eq. 2)then ! snow
            fall_velocity = 1.
!            fall_velocity = vvmax*0.2      ! adan change
        elseif(ipcp_type .eq. 3)then ! freezing rain
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! adan add
        elseif(ipcp_type .eq. 4)then ! sleet
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! adan add
        elseif(ipcp_type .eq. 5)then ! hail
            fall_velocity = 10.0
        endif

        if(p .ne. 101300. .or. t .ne. 273.15)then
            density_norm = (p / 101300.) * (273.15 / t)
            sqrt_density_norm = sqrt(density_norm)
        else
            sqrt_density_norm = 1.
        endif

        rate = rate * sqrt_density_norm

        return
        end


        subroutine cpt_concentration(rate,fall_velocity,concentration)

!       this gets m**3 of water per m**3 of volume (non-dimensional)

c lsw   stop if fall_velocity = 0
        if (fall_velocity.eq.0) then
          write(6,*) 'fall_velocity=',fall_velocity
          write(6,*)' error src/lib/pcpcnc.f, cpt_concentration, stop'
          stop
        else
          concentration = rate / fall_velocity

          rho = 1e3   ! mass in kilograms of a cubic meter of water

          concentration = concentration * rho
        endif

        return

        end
