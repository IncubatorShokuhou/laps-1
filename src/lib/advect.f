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

        subroutine advect(u,v,array_in,array_buf,grid_spacing_m
     1          ,imax,jmax,array_out,time,frac,lon,r_missing_data)

!       1992                            steve albers

        use mem_namelist, only: iverbose

        integer imax,jmax            ! input array dimensions
        real u(imax,jmax),v(imax,jmax) ! input wind field (wrt true north)
        real lon(imax,jmax)          ! input longitude field
        real array_in(imax,jmax)     ! input field to be advected
        real array_out(imax,jmax)    ! output field to be advected

        real array_buf(imax,jmax)    ! dummy array
        real ugrid(imax,jmax)
        real vgrid(imax,jmax)

        real time                    ! input seconds for advection
        real frac                    ! input scaling factor (normally 1.0)

        write(6,*)' calculating advected field sec/ratio:',time,frac

        seconds = time * frac

        do j = 1,jmax
        do i = 1,imax
            array_buf(i,j) = r_missing_data
        enddo
        enddo

!       rotate the winds to grid north
!       call uvtrue_to_uvgrid_2d for better efficiency
        call uvtrue_to_uvgrid_2d(u,v,ugrid,vgrid,lon,imax,jmax)

        do j = 1,jmax
        do i = 1,imax

            delta_u_m = ugrid(i,j) * seconds / grid_spacing_m
            delta_v_m = vgrid(i,j) * seconds / grid_spacing_m

!           write(6,*)delta_u_m,delta_v_m

            idelta_u_m = nint(delta_u_m)
            jdelta_v_m = nint(delta_v_m)

            inew = i + idelta_u_m
            jnew = j + jdelta_v_m

            if(   inew .ge. 1 .and. inew .le. imax
     1     .and.  jnew .ge. 1 .and. jnew .le. jmax ) then

!               array_buf(inew,jnew) = array_in(i,j)

!               preserve the maxima
                if(array_buf(inew,jnew) .eq. r_missing_data)then
                    array_buf(inew,jnew) = array_in(i,j)
                else
                    array_buf(inew,jnew) =
     1                  max(array_buf(inew,jnew),array_in(i,j))
                endif

!               if(array_in(i,j) .ge. 40.)then
!                   write(6,102)i,j,inew,jnew,idelta_u_m,jdelta_v_m,
!    1                          delta_u_m,delta_v_m,array_buf(inew,jnew)
102                 format(6i4,2f10.5,f8.1)
!               endif

            endif

        enddo
        enddo

!       fill in the cracks
        nfill1 = 0
        nfill2 = 0

        do j = 1,jmax
          jjmin = max(j-1,1)
          jjmax = min(j+1,jmax)

          do i = 1,imax
            if(array_buf(i,j) .eq. r_missing_data)then
                isum = 0
                ref_sum = 0.

                iimin = max(i-1,1)
                iimax = min(i+1,imax)

                do jj=jjmin,jjmax
                do ii=iimin,iimax
                    if(array_buf(ii,jj) .ne. r_missing_data)then
                        isum = isum + 1
                        ref_sum = ref_sum + array_buf(ii,jj)
                    endif
                enddo
                enddo

                if(isum .gt. 0)then
                    array_out(i,j) = ref_sum / float(isum)
                    nfill1 = nfill1 + 1
                else
                    array_out(i,j) = array_in(i,j) ! 0.
                    nfill2 = nfill2 + 1
                endif

c               write(6,*)' upgrade',array_buf(i,j),array_out(i,j)

            else ! array_buf(i,j) .ne. r_missing_data
                array_out(i,j) = array_buf(i,j)

            endif

          enddo ! i
        enddo ! j

        if(iverbose .ge. 2)then
            write(6,*)' nfill1,nfill2 ',nfill1,nfill2
        endif

        return
        end

        subroutine cpt_advection(field,u,v,dx,dy,ni,nj,advection)

        real field(ni,nj)
        real u(ni,nj)
        real v(ni,nj)
        real dx(ni,nj)
        real dy(ni,nj)
        real advection(ni,nj)

	do j=2,nj-1
	do i=2,ni-1
	  dth1 = (field(i,j) - field(i-1,j)) / dx(i,j)
	  dth2 = (field(i,j-1) - field(i-1,j-1)) / dx(i,j)
	  dtdx = (u(i,j-1) + u(i-1,j-1)) * (dth1 + dth2) * .25
	  dth3 = (field(i,j) - field(i,j-1)) / dy(i,j)
	  dth4 = (field(i-1,j) - field(i-1,j-1)) / dy(i,j)
	  dtdy = (v(i-1,j) + v(i-1,j-1)) * (dth3 + dth4) * .25
	  advection(i,j) = - dtdx - dtdy                       ! field / sec

        enddo ! i
        enddo ! j 

	call bounds(advection,ni,nj)

        return
        end
