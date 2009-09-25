cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis


        subroutine cpt_pcp_cnc(ref_3d,temp_3d,cldpcp_type_3d  ! Input
     1                                  ,ni,nj,nk     ! Input
     1                                  ,c_z2m        ! Input
     1                                  ,pcp_cnc_3d   ! Output (kg/m^3)
     1                                  ,rai_cnc_3d   ! Output (rain)
     1                                  ,sno_cnc_3d   ! Output (snow)
     1                                  ,pic_cnc_3d)  ! Output (precip ice)

        real temp_3d(ni,nj,nk)
        real ref_3d(ni,nj,nk)
        integer cldpcp_type_3d(ni,nj,nk)

        character*20 c_z2m

        real pcp_cnc_3d(ni,nj,nk)
        real rai_cnc_3d(ni,nj,nk)
        real sno_cnc_3d(ni,nj,nk)
        real pic_cnc_3d(ni,nj,nk)

        real rate_3d(ni,nj,nk)

!       Get 3D precip rates one layer at a time with ZR routine
        do k = 1,nk
            call zr(ref_3d(1,1,k),ni,nj,rate_3d(1,1,k))
        enddo ! k

        if(c_z2m .eq. 'albers')then
          do j = 1,nj
          do i = 1,ni
          do k = 1,nk
            pressure = pressure_of_level(k)
            ipcp_type = cldpcp_type_3d(i,j,k) / 16   ! Pull out precip type
            if(ipcp_type .ne. 0)then
!                call cpt_fall_velocity(ipcp_type,pressure,temp_3d(i,j,k)
!     1                                                  ,fall_velocity)
! Adan add
                call cpt_fall_velocity(ipcp_type,pressure,temp_3d(i,j,k)
     1                                 ,ref_3d(i,j,k),fall_velocity)
                call cpt_concentration(rate_3d(i,j,k),fall_velocity
     1                                          ,pcp_cnc_3d(i,j,k))

                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    rai_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    sno_cnc_3d(i,j,k) = 0.
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 2)then                       ! snow
                    rai_cnc_3d(i,j,k) = 0.
                    sno_cnc_3d(i,j,k) = pcp_cnc_3d(i,j,k)
                    pic_cnc_3d(i,j,k) = 0.

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! IP or Hail
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

        elseif(c_z2m .eq. 'kessler')then 
          do j = 1,nj
          do i = 1,ni
          do k = 1,nk

!           Compute the basic reflectivity 
!           refl(i,j,k) =17300.0 * &
!                 (rho(i,j,k) * 1000.0 * &
!                  MAX(0.0,rainmr(i,j,k)))**svnfrth

!           Add the ice component
!           refl(i,j,k)=refl(i,j,k) + &
!                 38000.0*(rho(i,j,k) * 1000.0 * &
!                 MAX(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

            if(ipcp_type .ne. 0)then
                if(ipcp_type .eq. 1 .or. ipcp_type .eq. 3)then     ! rain or zr
                    a = 17300.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/a)**(4./7.)
                    pcp_cnc_3d(i,j,k) = rho_q
                    rai_cnc_3d(i,j,k) = rho_q

                elseif(ipcp_type .eq. 2)then                       ! snow
                    b = 38000.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/b)**(1./2.2)
                    pcp_cnc_3d(i,j,k) = rho_q
                    sno_cnc_3d(i,j,k) = rho_q

                elseif(ipcp_type .eq. 4 .or. ipcp_type .eq. 5)then ! IP or Hail
                    b = 38000.
                    z = 10.**(ref_3d(i,j,k)/10.)
                    rho_q = (z/b)**(1./2.2)
                    pcp_cnc_3d(i,j,k) = rho_q
                    pic_cnc_3d(i,j,k) = rho_q

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

        endif

        end

        subroutine cpt_fall_velocity(ipcp_type,p,t,dbz,fall_velocity)

        dbz_eff = max(dbz,1.0)
        vvmax = 4.32*dbz_eff**0.0714286      ! Adan add

        if(ipcp_type .eq. 1)then ! Rain
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! Adan add
        elseif(ipcp_type .eq. 2)then ! Snow
            fall_velocity = 1.
!            fall_velocity = vvmax*0.2      ! Adan change
        elseif(ipcp_type .eq. 3)then ! Freezing Rain
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! Adan add
        elseif(ipcp_type .eq. 4)then ! Sleet
!            fall_velocity = 5.0
            fall_velocity = vvmax    ! Adan add
        elseif(ipcp_type .eq. 5)then ! Hail
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

!       This gets m**3 of water per m**3 of volume (non-dimensional)

C LSW   STOP if fall_velocity = 0
        if (fall_velocity.eq.0) then
          write(6,*) 'fall_velocity=',fall_velocity
          write(6,*)' Error src/lib/pcpcnc.f, cpt_concentration, STOP'
          stop
        else
          concentration = rate / fall_velocity

          rho = 1e3   ! mass in kilograms of a cubic meter of water

          concentration = concentration * rho
        endif

        return

        end
