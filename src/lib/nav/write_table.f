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
	subroutine write_table (table_path,
     &             nx,ny,lat,lon,ri,rj,istatus)

	integer nx,ny
	real lat (nx,ny)
	real lon (nx,ny)
	real ri (nx,ny)
	real rj (nx,ny)
	integer istatus
        character*(*) table_path

	istatus = 0
        n=index(table_path,' ')

	open (unit=12,file=table_path(1:n-1),
     &form='unformatted',status='unknown',err=23)

	write(12,err=24) lat
	write(12,err=24) lon
	write(12,err=24) ri
	write(12,err=24) rj

	close (12)

	istatus = 1

        goto 100

23      write(6,*)'error openning ',table_path(1:n-1)
        goto 100

24      write(6,*)'error writing to ',table_path(1:n-1)

100     return
        end
