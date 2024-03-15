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
       subroutine satdat2laps_ir(imax,jmax,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  image_ir,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  line_dim,elem_dim, ! input array dimensions
     &                  sa,sc,st,
     &                  istatus)

c.....  this routine can be called for any goes ir band
c.....  this is mcalbers version returing warmest and coldest pixels
c.....  the warmest go into st, the "extremum" to sc, the mean to sa
c       sc is the result of an edge enhancing filter that tries to filter
c       out inbetween pixels in favor of either the warmest or the coldest.
c
c       j. smart          jul 1995          orginal subroutine for goes 8
c.....          changes:  02-oct-1990       set up for prodgen.
c.....          changes:     sep-1993       add average (sa)
c
	implicit none

	integer max_elem,max_line
        integer imax,jmax
	parameter (max_elem = 15)
	parameter (max_line = 15)
	integer line_dim,elem_dim
        real    r_grid_ratio
        real    rgrid_ratio_thresh

	real image_ir(elem_dim,line_dim)
        real st(imax,jmax)
        real sc(imax,jmax)
	real sa(imax,jmax)
        real r_llij_lut_ri(imax,jmax)
        real r_llij_lut_rj(imax,jmax)

	real line_mx,line_mn,elem_mx,elem_mn
	real t_array(max_elem*max_line)
	real wm, wc, btemp, tmean
        real frac
c       real fraci,fracj
        real pixsum 
        real r_missing_data,rmin,rmax
        real pcnt_msng_thresh
        real result

        integer i,j,ii,jj,iir,jir,max_wi,max_wj
        integer istart,jstart
        integer iend,jend
	integer npix, nwarm
        integer maxpix
        integer ipix
        integer istatus
        integer qcstatus
        integer fcount
        integer insufdata
        integer icnt

        logical lforce_switch

        call zero(st,imax,jmax)
        istatus = -1
        qcstatus=0
        fcount=0
c
        write(6,*)' subroutine satdat2laps_ir'
        write(6,*)'   i   j   warmpix  coldpix  npix nwarm  cldtemp'
c
c the "10" loop represents input image resolution < output grid resolution such
c that there are enough pixels from the input image to get a representative
c mean value for the remapped output grid value. the "10" loop is also used
c when we have significant missing data in the input image since the "gdtost"
c routine can have deleterious effects.
c
        pcnt_msng_thresh=0.05
        rgrid_ratio_thresh=1.5
        lforce_switch=.false.
        icnt=0
        do j=1,line_dim
        do i=1,elem_dim
           if(image_ir(i,j).eq.r_missing_data)then
              icnt=icnt+1
           endif
        enddo
        enddo

        call array_range(image_ir,elem_dim,line_dim,rmin,rmax
     1                  ,r_missing_data)
        write(6,*)' image_ir input (non-missing) range is ',rmin,rmax

        if(icnt.gt.(pcnt_msng_thresh*elem_dim*line_dim))then
           lforce_switch=.true.
           print*,'more than 5% of data missing: '
     &           ,float(icnt)/float(imax*jmax)
           print*,'thus force grid point averaging in satir2laps'
        else
           print*,'less than 5% of data missing: '
     &           ,float(icnt)/float(imax*jmax)
        endif
 
        insufdata=0
        if(r_grid_ratio .le.0.0)then
            write(6,*)'warning: r_grid_ratio out of range',r_grid_ratio
            goto 1000
        endif

        if(r_grid_ratio.lt.rgrid_ratio_thresh.or.lforce_switch)then

          if(r_grid_ratio.lt.rgrid_ratio_thresh)then
            write(6,*)'grid ratio < thresh',r_grid_ratio
     1                                     ,rgrid_ratio_thresh
            write(6,*)'use pixel avg to get ir tb (coarse laps grid)'
          else
            write(6,*)'use pixel avg to get ir tb (due to force switch)'
          endif

          max_wi=10.
          max_wj=10.

          do 10 j=1,jmax
          do 10 i=1,imax

           if(st(i,j).ne.0.) go to 10
c
c line/elem are floating point i/j positions in ispan grid for input lat/lon
c also, use the lat/lon to real i/j look up table (r_llij_lut) to map out points
c needed for satellite pixels.
c****************************************************************************

           if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &        r_llij_lut_rj(i,j).ne.r_missing_data)then

             elem_mx = r_llij_lut_ri(i,j) + ((1./r_grid_ratio) * 0.5)
             elem_mn = r_llij_lut_ri(i,j) - ((1./r_grid_ratio) * 0.5)
             line_mx = r_llij_lut_rj(i,j) + ((1./r_grid_ratio) * 0.5)
             line_mn = r_llij_lut_rj(i,j) - ((1./r_grid_ratio) * 0.5)
             jstart = nint(line_mn) !+0.5)
             jend   = nint(line_mx)
             istart = nint(elem_mn) !+0.5)
             iend   = nint(elem_mx)

             if(istart.le.0 .or. jstart.le.0 .or.
     &iend.gt.elem_dim .or. jend.gt.line_dim)then
                insufdata=insufdata+1
c            write(*,*)'insufficient data for lat/lon sector'
c               write(*,1020)i,j
c1020	        format(1x,'laps grid (i,j) = ',i3,1x,i3)
c               write(6,1021)elem_mx,elem_mn,line_mx,line_mn
c1021            format(1x,'elem mx/mn  line mx/mn ',4f7.1)
             else
c
c **** find pixels around grid point
c
                wm=-1.e15
                wc=1.e15
                npix = 0
                pixsum = 0.
                maxpix = 0

                do 3 jj=jstart,jend
                   do 3 ii=istart,iend

                      if(image_ir(ii,jj) .ne.
     &r_missing_data)then
                         npix=npix+1
                         t_array(npix) = image_ir(ii,jj)

                      endif

    3           continue  
c
                if(npix.gt.1)then

c
c...  this section finds the warmest pixel, coldest pixel, and mean pixel temp.
c
                   do ii=1,npix

                      btemp  = t_array(ii)
                      if(btemp .le. 0.)then
                         write(6,*)
     1                  ' error: btemp too low in satdat2laps_ir ',btemp
                         return
                      endif

                      pixsum = pixsum + btemp

                      if(btemp.gt.wm) then
                         wm=btemp
                      end if

                      if(btemp.lt.wc) then
                         wc=btemp
                      end if

                   enddo

                elseif(npix.eq.1)then

                   btemp = t_array(npix)

                else   

!                  btemp = image_ir(ii,jj)
                   btemp=r_missing_data
                   fcount=fcount+1

                endif

                if(npix .ge. maxpix)maxpix=npix

                if(npix .gt. 1)then

                   sa(i,j) = pixsum / float(npix)

c.....  operate on t_array to find cloud top temp. this chooses the warmest
c.....  or coldest pixel depending on whether most pixels are closer to the
c.....  warmest or coldest. the net result is an edge sharpening filter.

                   tmean = (wc + wm)/2.

                   nwarm = 0
                   do ipix = 1,npix
                      if(t_array(ipix) .ge. tmean)then
                         nwarm = nwarm + 1
                      endif
                   enddo

c              write(6,1112) wm,wc
 1112          format(1x,11f7.0)

                   frac = float(nwarm) / float(npix)
                   if(frac .gt. 0.5)then
                      sc(i,j)=wm
                   else
                      sc(i,j)=wc
                   endif

                   st(i,j)=wm

                else

                   sa(i,j)=btemp
                   sc(i,j)=btemp
                   st(i,j)=btemp

                endif ! npix .gt. 1

             end if  ! enough data for num_lines .gt. 0

           endif   !r_llij's .ne. r_missing_data

cisid      if(sa(i,j).eq.r_missing_data)then
           if((sa(i,j).eq.r_missing_data).and.(max_wi.lt.10)
     &         .and.(max_wj.lt.10)) then
                print*,'found missing data in ir remapping'
           endif

           if(i .eq. i/200*200 .and. j .eq. j/200*200       
     1        .or.                   sa(i,j) .lt. 190.      
     1        .or.                   
     1       (sa(i,j) .gt. 500. .and. sa(i,j) .ne. r_missing_data)
     1                                 )then
                write(6,5555)i,j,wm,wc,npix,nwarm,sa(i,j)
     1                      ,r_llij_lut_ri(i,j),r_llij_lut_rj(i,j)
5555            format(1x,2i4,2f10.2,2i5,f10.2,2x,2f10.2)
           endif

   10     continue ! i,j
          write(6,*)'max num ir pix for avg: ',maxpix
          write(6,*)'number of laps gridpoints missing',
     &fcount

          call array_range(sa,imax,jmax,rmin,rmax,r_missing_data)
          write(6,*)' sa (non-missing) range is ',rmin,rmax

c       elseif(r_grid_ratio .lt. 1.0)then
c       ---------------------------------------

c input image resolution is somewhat large relative to output grid.
c this section uses bilinear interpolation to map
c the four surrounding input pixels to the output grid.
c

c         write(6,*)'image res .ge. output grid spacing'
c         write(6,*)'using bilinear interp for ir tb'
c         do 20 j=1,jmax
c         do 20 i=1,imax

c            if(i_found_it .ne. 1)then
c              write(6,*)'enter i and j point'
c              read(5,*)ipt,jpt
c              i_found_it = 1
c            end if
c           if(st(i,j).ne.0.) go to 20
c
c line/elem are floating point i/j positions for input lat/lon
c
c bilinear_interp_extrap checks on boundary conditions and
c uses r_missing_data if out of bounds.
c

c           if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
c    &         r_llij_lut_rj(i,j).ne.r_missing_data)then

c              call  bilinear_interp_extrap(
c    &              r_llij_lut_ri(i,j),
c    &              r_llij_lut_rj(i,j),
c    &              elem_dim,line_dim,image_ir,
c    &              result,istatus)

c              if(result .ne. r_missing_data .and.
c    &            result .gt. 0.0)then

c                 sa(i,j) = result

c              else

c                 sa(i,j) = r_missing_data

c              endif
c              sc(i,j)=sa(i,j)
c              st(i,j)=sa(i,j)

c            if(i.eq.ipt .and. j.eq.jpt)then
c               write(29,39)ipt,jpt
c               write(29,49)image_ir(ii,jj)
c               i_found_it = -1
c 39          format(1x,'ipt = ',i3,2x,'jpt = ',i3)
c 49      format(1x,'ir_count = ',1x,f6.1,/,'--------------------',//)
c            end if

c           endif  !lut's = r_missing_data

c  20     continue ! i,j
c          close(29)

        else      ! here the input data resolution is course relative to
c                   the analysis domain

          print*,'fine laps grid: use spline - gdtost'   

          do 30 j=1,jmax
          do 30 i=1,imax

           if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &        r_llij_lut_rj(i,j).ne.r_missing_data)then

               call gdtost(image_ir,elem_dim,line_dim,
     .             r_llij_lut_ri(i,j),r_llij_lut_rj(i,j)
     .            ,sa(i,j),.false.)
           
               sc(i,j)=sa(i,j)
               st(i,j)=sa(i,j)

           endif

   30     continue

        end if     ! r_image_ratio

c        write(6,1234) ib,i4vtime,ict
c1234       format(1x,'band ',i4,' count for i4time ',i10,' is ',i8)

        if(insufdata.gt.0)then
           print*,'found ',insufdata,' points that are too'
           print*,'close to data edge to compute average'
        endif
        istatus = 1

        goto 1000

999     print*,'sorry bad value for r_grid_ratio ', r_grid_ratio
c
1000    return
        end
