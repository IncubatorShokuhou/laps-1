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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine mak_cld_grid (i4time_d,i4time_r,cg,ii,jj,kk,
     1                            lct,istatus)

c       $log: mak_cld_grid.for,v $
c revision 1.6  1995/09/13  21:36:15  birk
c added disclaimer to files
c
c revision 1.5  1994/09/09  15:53:09  birk
c  corrected mixed mode problem in pressure level assignment
c
c revision 1.4  1994/07/22  22:03:39  birk
c made code compatible with lapsparms.inc
c
c revision 1.3  1994/06/17  22:35:32  birk
c modification to peg cloud amount to upper limit of 1.0
c
c revision 1.2  1994/06/15  22:22:54  birk
c improved initialization
c
c revision 1.1  1994/04/25  15:11:00  birk
c initial revision
c

        implicit none


        include 'laps_cloud.inc'

	integer cloud_levels
	parameter (cloud_levels = kcloud)

	character*31 cloud_ext
	data cloud_ext/'lc3'/

c parameter variables

      integer ii,jj,kk
      integer lct

      integer i4time_d     !i4time of the desired grid
      integer i4time_r       !i4time of returned grid
      integer istatus
      real cg(ii,jj,kk)    !cloud grid (1=clear, 0=cloudy)

c variables dependent on dynamic allocation

      integer lvl(kk)
      real lma_ht(ii,jj,kk)
      real heights(ii,jj,kk)
      real clouds_3d(ii,jj,cloud_levels)
      real cld_pres(cloud_levels)

c internal variables

      integer i,j,k,l
      real cld_amt
      character*3 desired_field


c       threshold is dtermined by variable thresh
      save thresh
      real thresh
      data thresh /.65/
      real pressure_of_level !function call

c       initialize

        do k = 1,kk
        lvl(k) =  nint( pressure_of_level(k)  * .01 )
        enddo

        do k = 1,kk
        do j = 1,jj
        do i = 1,ii

        cg(i,j,k) = 0.  ! clear

        enddo
        enddo
        enddo

        desired_field = 'ht  '

        call get_maps_df (i4time_d,desired_field,lma_ht,
     1              ii,jj,kk,lct,istatus)

        if (istatus .ne. 1 ) return

        do k = 1,kk

        call bilinear_interp(ii,jj,lma_ht(1,1,k),ii,
     1   jj,heights(1,1,k),ii,jj)

        enddo


        call get_clouds_3dgrid(i4time_d,i4time_r,
     1  ii,jj,cloud_levels,cloud_ext,
     1  clouds_3d,cld_hts,cld_pres,istatus)


        if(istatus .ne. 1 ) return

        do k = 1,kk
        do j = 1,jj
        do i = 1,ii

        do l = 1,cloud_levels-1

        if (heights(i,j,k).gt.cld_hts(l)
     1  .and. heights(i,j,k).lt.cld_hts(l+1)) then

        call interp (heights(i,j,k),cld_hts(l),cld_hts(l+1),
     1  clouds_3d(i,j,l),clouds_3d(i,j,l+1),cld_amt)


c   ... added mod to limit cloud fraction to upper limit of 1.0

        if(cld_amt .gt. 1.) then
                cld_amt = 1
        endif

        cg(i,j,k) = cld_amt

        endif

        enddo

        enddo
        enddo
        enddo


        return
        end
