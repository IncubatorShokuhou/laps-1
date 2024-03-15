        subroutine write_bal_laps(i4time, phi, u, v, temp, om, rh
        ., sh, imax, jmax, kmax, p, istatus)
        c
        implicit none
        c
        integer i4time,
1       imax, jmax, kmax,
1       kmax3, kmax2,
1       ip(kmax),
1       lvl(kmax*3),
1       i, j, k,
1       istatus
        c
        integer lend

        real phi(imax, jmax, kmax),
1       u(imax, jmax, kmax),
1       v(imax, jmax, kmax),
1       om(imax, jmax, kmax),
1       temp(imax, jmax, kmax),
1       bal(imax, jmax, kmax*3),
1       rh(imax, jmax, kmax),
1       sh(imax, jmax, kmax),
1       p(kmax)
        c
        character*150 dir
        character*150 directory
        character*31 ext
        character*3 var(kmax*3)
        character*4 lvl_coord(kmax*3)
        character*9 fname9
        character*10 units(kmax*3)
        character*125 comment(kmax*3)
        c
        c - ------------------------------------------------------------------------------
        c
        c first is the wind(lw3)
        c

        ext = 'lw3'
        call get_directory(ext, directory, lend)
        lend = lend - 4  ! get_directory insures that a "/" is at the end of the directory string.

        dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
        c
        do k = 1, kmax
           ip(k) = int(p(k)/100.)
           var(k) = 'u3 '
           lvl(k) = ip(k)
           lvl_coord(k) = 'mb  '
           units(k) = 'm/s   '
           comment(k) = 'non-linear balanced u-component wind.'
           do j = 1, jmax
           do i = 1, imax
              bal(i, j, k) = u(i, j, k)
           end do
           end do
        end do
        c
        do k = kmax + 1, 2*kmax
           var(k) = 'v3 '
           lvl(k) = ip(k - kmax)
           lvl_coord(k) = 'mb  '
           units(k) = 'm/s   '
           comment(k) = 'non-linear balanced v-component wind.'
           do j = 1, jmax
           do i = 1, imax
              bal(i, j, k) = v(i, j, k - kmax)
           end do
           end do
        end do
        c
        do k = 2*kmax + 1, 3*kmax
           var(k) = 'om '
           lvl(k) = ip(k - 2*kmax)
           lvl_coord(k) = 'mb  '
           units(k) = 'pa/s  '
           comment(k) = 'non-linear balanced omega.           '
           do j = 1, jmax
           do i = 1, imax
              bal(i, j, k) = om(i, j, k - 2*kmax)
           end do
           end do
        end do

        call make_fnam_lp(i4time, fname9, istatus)
        write (6, *) ' writing grids ', ext(1:3), ' ', fname9

        kmax3 = kmax*3
        call write_laps_data(i4time, dir, ext, imax, jmax
        +, kmax3, kmax3, var, lvl, lvl_coord, units, comment, bal, istatus)
        c
        c - ---------------------
        c now lt1

        ext = 'lt1'
        dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
        do k = 1, kmax
           var(k) = 'ht '
           lvl(k) = ip(k)
           lvl_coord(k) = 'mb  '
           units(k) = 'meters'
           comment(k) = 'non-linear balanced height.'
           do j = 1, jmax
           do i = 1, imax
              bal(i, j, k) = phi(i, j, k)
           end do
           end do
        end do
        c
        do k = kmax + 1, 2*kmax
           var(k) = 't3'
           lvl(k) = ip(k - kmax)
           lvl_coord(k) = 'mb  '
           units(k) = 'kelvin'
           comment(k) = 'non-linear balanced temp.'
           do j = 1, jmax
           do i = 1, imax
              bal(i, j, k) = temp(i, j, k - kmax)
           end do
           end do
        end do
        c
        kmax2 = 2*kmax
        c
        write (6, *) ' writing grids ', ext(1:3), ' ', fname9

        call write_laps_data(i4time, dir, ext, imax, jmax, kmax2, kmax2, var, lvl
1       , lvl_coord, units, comment, bal, istatus)
        c
        if (istatus .ne. 1) then
           print *, 'error writing balanced data.'
           istatus = 0
           return
        end if
        c
        c next, the rh.
        c
        ext = 'lh3'
        dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
        do k = 1, kmax
           var(k) = 'rhl'
           lvl(k) = ip(k)
           lvl_coord(k) = 'mb  '
           units(k) = 'percent'
           comment(k) = 'balanced rh'
        end do

        write (6, *) ' writing grids ', ext(1:3), ' ', fname9

        call write_laps_data(i4time, dir, ext, imax, jmax
        +, kmax, kmax, var, lvl, lvl_coord, units, comment, rh, istatus)

        if (istatus .ne. 1) then
           print *, 'error writing balanced rh data.'
           istatus = 0
           return
        end if
        c
        c finally, the specific humidity
        c
        ext = 'lq3'
        dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
        do k = 1, kmax
           var(k) = 'sh '
           lvl(k) = ip(k)
           lvl_coord(k) = 'mb  '
           units(k) = 'kg/kg'
           comment(k) = 'balanced specific humidity.'
        end do

        write (6, *) ' writing grids ', ext(1:3), ' ', fname9

        call write_laps_data(i4time, dir, ext, imax, jmax
        +, kmax, kmax, var, lvl, lvl_coord, units, comment, sh, istatus)

        if (istatus .ne. 1) then
           print *, 'error writing balanced sh data.'
           istatus = 0
           return
        end if
        c
        istatus = 1
        return
        c
     end
!c
!============
! hongli jiang add to write sigma_ht output. 11/2/2011
!========
     subroutine write_bal_laps_ht(i4time, p3, u, v, temp, w, rh
     ., sh, imax, jmax, kmax, p, istatus)
     c
     implicit none
     c
     integer i4time,
1    imax, jmax, kmax,
1    kmax3, kmax2,
1    ip(kmax),
1    lvl(kmax*3),
1    i, j, k,
1    istatus
     c
     integer lend

     real p3(imax, jmax, kmax),
1    u(imax, jmax, kmax),
1    v(imax, jmax, kmax),
1    w(imax, jmax, kmax),
1    temp(imax, jmax, kmax),
1    bal(imax, jmax, kmax*3),
1    rh(imax, jmax, kmax),
1    sh(imax, jmax, kmax),
1    p(kmax)
     c
     character*150 dir
     character*150 directory
     character*31 ext
     character*3 var(kmax*3)
     character*4 lvl_coord(kmax*3)
     character*9 fname9
     character*10 units(kmax*3)
     character*125 comment(kmax*3)
     c
     c - ------------------------------------------------------------------------------
     c
     c first is the wind(lw3)
     c

     ext = 'lw3'
     call get_directory(ext, directory, lend)
     lend = lend - 4  ! get_directory insures that a "/" is at the end of the directory string.

     dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
     c
     do k = 1, kmax
!hj                ip(k)=int(p(k)/100.)
        ip(k) = int(p(k))
        var(k) = 'u3 '
        lvl(k) = ip(k)
        lvl_coord(k) = 'mb  '
        units(k) = 'm/s   '
        comment(k) = 'non-linear balanced u-component wind.'
        do j = 1, jmax
        do i = 1, imax
           bal(i, j, k) = u(i, j, k)
        end do
        end do
     end do
     c
     do k = kmax + 1, 2*kmax
        var(k) = 'v3 '
        lvl(k) = ip(k - kmax)
        lvl_coord(k) = 'mb  '
        units(k) = 'm/s   '
        comment(k) = 'non-linear balanced v-component wind.'
        do j = 1, jmax
        do i = 1, imax
           bal(i, j, k) = v(i, j, k - kmax)
        end do
        end do
     end do
     c
     do k = 2*kmax + 1, 3*kmax
        var(k) = 'w3 '
        lvl(k) = ip(k - 2*kmax)
        lvl_coord(k) = 'mb  '
        units(k) = 'pa/s  '
        comment(k) = 'non-linear balanced w.           '
        do j = 1, jmax
        do i = 1, imax
           bal(i, j, k) = w(i, j, k - 2*kmax)
        end do
        end do
     end do

     call make_fnam_lp(i4time, fname9, istatus)
     write (6, *) ' writing grids ', ext(1:3), ' ', fname9

     kmax3 = kmax*3
     call write_laps_data(i4time, dir, ext, imax, jmax
     +, kmax3, kmax3, var, lvl, lvl_coord, units, comment, bal, istatus)
     c
     c - ---------------------
     c now lt1

     ext = 'lt1'
     dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
     do k = 1, kmax
        var(k) = 'p3 '
        lvl(k) = ip(k)
        lvl_coord(k) = 'mb  '
        units(k) = 'meters'
        comment(k) = 'non-linear balanced height.'
        do j = 1, jmax
        do i = 1, imax
           bal(i, j, k) = p3(i, j, k)
        end do
        end do
     end do
     c
     do k = kmax + 1, 2*kmax
        var(k) = 't3'
        lvl(k) = ip(k - kmax)
        lvl_coord(k) = 'mb  '
        units(k) = 'kelvin'
        comment(k) = 'non-linear balanced temp.'
        do j = 1, jmax
        do i = 1, imax
           bal(i, j, k) = temp(i, j, k - kmax)
        end do
        end do
     end do
     c
     kmax2 = 2*kmax
     c
     write (6, *) ' writing grids ', ext(1:3), ' ', fname9

     call write_laps_data(i4time, dir, ext, imax, jmax, kmax2, kmax2, var, lvl
1    , lvl_coord, units, comment, bal, istatus)
     c
     if (istatus .ne. 1) then
        print *, 'error writing balanced data.'
        istatus = 0
        return
     end if
     c
     c next, the rh.
     c
     ext = 'lh3'
     dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
     do k = 1, kmax
        var(k) = 'rhl'
        lvl(k) = ip(k)
        lvl_coord(k) = 'mb  '
        units(k) = 'percent'
        comment(k) = 'cloud liquid balanced rh.'
     end do

     write (6, *) ' writing grids ', ext(1:3), ' ', fname9

     call write_laps_data(i4time, dir, ext, imax, jmax
     +, kmax, kmax, var, lvl, lvl_coord, units, comment, rh, istatus)

     if (istatus .ne. 1) then
        print *, 'error writing balanced rh data.'
        istatus = 0
        return
     end if
     c
     c finally, the specific humidity
     c
     ext = 'lq3'
     dir = directory(1:lend)//'balance/'//ext(1:3)//'/'
     do k = 1, kmax
        var(k) = 'sh '
        lvl(k) = ip(k)
        lvl_coord(k) = 'mb  '
        units(k) = 'kg/kg'
        comment(k) = 'cloud liquid balanced specific humidity.'
     end do

     write (6, *) ' writing grids ', ext(1:3), ' ', fname9

     call write_laps_data(i4time, dir, ext, imax, jmax
     +, kmax, kmax, var, lvl, lvl_coord, units, comment, sh, istatus)

     if (istatus .ne. 1) then
        print *, 'error writing balanced sh data.'
        istatus = 0
        return
     end if
     c
     istatus = 1
     return
     c
  end
