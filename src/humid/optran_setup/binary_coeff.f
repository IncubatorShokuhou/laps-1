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
        implicit none
	include 'constants_optran.inc'

	real*8 bw(nwet+1,nw,nchan)
	real*8 bd(ndry+1,nw,nchan)
	real*8 bo(nozo+1,nw,nchan)

	integer i,j,k,ichan

	open(12,file='coef.dat',
     &		form='formatted')
	open(42,file='wet_coeff.dat',
     &          form='unformatted')
	open(43,file='dry_coeff.dat',
     &          form='unformatted')
	open(44,file='ozo_coeff.dat',
     &          form='unformatted')

	do ichan = 1 , nchan
1200     format(i4,6e20.12)

	 do i = 1 , nw
	  read(12,1200) k,(bd(j,i,ichan),j=1,ndry+1)
	 enddo
	 do i = 1 , nw
	  read(12,1200) k,(bw(j,i,ichan),j=1,nwet+1)
  	 enddo
	 do i = 1 , nw
	  read(12,1200) k,(bo(j,i,ichan),j=1,nozo+1)
  	 enddo

	enddo

	write(42) bw
	write(43) bd
	write(44) bo

	close(42)
	close(43)
	close(44)


	stop
	end
