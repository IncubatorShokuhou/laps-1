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
c
c
        subroutine add(a,b,result,imax,jmax)
c
cdoc    routine to add array 'b' to array 'a' and put the
cdoc    result into array 'result' .
c
        real a(imax,jmax), b(imax,jmax), result(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          result(i,j) = a(i,j) + b(i,j)
        enddo !i
        enddo !j
c
        return
        end
c
c
        subroutine add_miss(a,b,result,imax,jmax)
c
cdoc    routine to add array 'b' to array 'a' and put the
cdoc    result into array 'result'. this takes account of 'r_missing_data'.
c
        real a(imax,jmax), b(imax,jmax), result(imax,jmax)
c
        call get_r_missing_data(r_missing_data,istatus)

        do j=1,jmax
        do i=1,imax
          if(a(i,j) .ne. r_missing_data .and. 
     1       b(i,j) .ne. r_missing_data       )then
              result(i,j) = a(i,j) + b(i,j)
          else
              result(i,j) = r_missing_data
          endif
        enddo !i
        enddo !j
c
        return
        end
