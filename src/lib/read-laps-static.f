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
c read time-invariant laps data.  i4time is set to zero which results
c in a "file time" of jan 1, 1960 or 1970.
c
        subroutine read_laps_static (dir, ext, idim, jdim, kmax, kdim,
     1                          var, units, comment, grid, status)
c
        implicit none
c
        character*(*) dir, ext
        integer idim, jdim, kmax, kdim
        character*(*) var(1), units(1), comment(1)
        real grid(idim,jdim,kdim)
        integer status
c
        integer nvarsmax, i
        parameter (nvarsmax=20)         ! may need to increase this someday
        character*4 lvl_coord(nvarsmax)
        integer lvl(nvarsmax)

        write(6,*)
     1 ' warning: this routine reads static data in the old format.'
        write(6,*)
     1 ' the format is non-netcdf and is not portable.'
        write(6,*)
     1 ' it is recommended to call rd_laps_static instead to read'
        write(6,*)
     1 ' the new, portable netcdf files.'

        do i=1,nvarsmax
                lvl(i) = 0
        enddo

        call read_laps_data (0, dir, ext, idim, jdim, kmax, kdim,
     1  var, lvl, lvl_coord, units, comment, grid, status)

        return
        end
