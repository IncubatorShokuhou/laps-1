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
       subroutine satdat2laps_sndr(imax,jmax,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  image_sndr,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  line_dim,elem_dim, ! input array dimensions
     &                  sa,sc,st,
     &                  istatus)

c.....  this routine can be called for any goes sounder channel data.
c.....  this is mcalbers version returing warmest and coldest pixels
c.....  the warmest go into st, the "extremum" to sc, the mean to sa
c       sc is the result of an edge enhancing filter that tries to filter
c       out inbetween pixels in favor of either the warmest or the coldest.
c
c       j. smart          jul 1995          orginal subroutine for goes 8
c.....          changes:  02-oct-1990       set up for prodgen.
c.....          changes:     sep-1993       add average (sa)
c.....  j. smart          dec 1996          modified the satdat2laps_ir for sounder data
c
      implicit none

      integer    max_elem,max_line
      integer    imax,jmax
      parameter (max_elem = 15)
      parameter (max_line = 15)
      integer    line_dim,elem_dim
      real     r_grid_ratio
      real     image_sndr(elem_dim,line_dim)
      real     sc(imax,jmax)
      real     sa(imax,jmax)
      real     st(imax,jmax)
      real r_llij_lut_ri(imax,jmax)
      real r_llij_lut_rj(imax,jmax)

      real line_mx,line_mn,elem_mx,elem_mn
      real t_array(max_elem*max_line)
      real wm, wc, btemp, tmean
      real frac
c     real fraci,fracj
      real pixsum 
      real r_missing_data
      real result

        integer i,j,ii,jj
        integer i1,j1,i2,j2
        integer istart,jstart
        integer iend,jend
        integer npix, nwarm
        integer maxpix
        integer ipix
        integer istatus
        integer qcstatus
        integer fcount
        integer icnt_out

        call initialize(sa,imax,jmax,r_missing_data)
        call initialize(sc,imax,jmax,r_missing_data)
        call initialize(st,imax,jmax,r_missing_data)

        istatus = -1
        qcstatus=0
        fcount=0
c
c       write(6,*)'   i   j   warmpix  coldpix  npix nwarm  cldtemp'
c
c the "10" loop represents input image resolution < output grid resolution such
c that there are enough pixels from the input image to get a representative
c mean value for the remapped output grid value
c

        icnt_out=0
        if(r_grid_ratio .lt. 0.5)then  !0.75)then

          write(6,*)'grid ratio .lt. 0.5'   !0.75'
          write(6,*)'use pixel avg to get sndr rad'

          do 10 j=1,jmax
          do 10 i=1,imax

             if(st(i,j).ne.0.) go to 10
c
c line/elem are floating point i/j positions in ispan grid for input lat/lon
c also, use the lat/lon to real i/j look up table (r_llij_lut) to map out points
c needed for satellite pixels.
c****************************************************************************

c            if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
c    &          r_llij_lut_rj(i,j).ne.r_missing_data)then

             elem_mx = r_llij_lut_ri(i,j) + ((1./r_grid_ratio) * 0.5)
             elem_mn = r_llij_lut_ri(i,j) - ((1./r_grid_ratio) * 0.5)
             line_mx = r_llij_lut_rj(i,j) + ((1./r_grid_ratio) * 0.5)
             line_mn = r_llij_lut_rj(i,j) - ((1./r_grid_ratio) * 0.5)
             jstart = nint(line_mn+0.5)
             jend   = int(line_mx)
             istart = nint(elem_mn+0.5)
             iend   = int(elem_mx)

             if(istart.le.0 .or. jstart.le.0 .or.
     &iend.gt.elem_dim .or. jend.gt.line_dim)then
             icnt_out=icnt_out+1

c            write(*,*)'insufficient data for lat/lon sector'
c               write(*,1020)i,j
c1020              format(1x,'laps grid (i,j) = ',i3,1x,i3)
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

                      if(image_sndr(ii,jj) .ne.
     &r_missing_data)then
                         npix=npix+1
                         t_array(npix) = image_sndr(ii,jj)

                      endif

    3           continue  
c
                if(npix.gt.1)then

c
c...  this section finds the warmest pixel, coldest pixel, and mean pixel temp.
c
                   do ii=1,npix

                      btemp  = t_array(ii)
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

c            endif

             endif  ! enough data for num_lines .gt. 0

c            if(i .eq. i/10*10 .and. j .eq. j/10*10)then
c               write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
c 5555         format(1x,2i4,2f10.2,2i5,f10.2)
c            endif

   10     continue ! i,j

          print*,'max num sndr pix for avg: ',maxpix
          print*,'number of laps gridpoints missing',
     &fcount

        else
c       ---------------------------------------
c input image resolution is large relative to output grid spacing
c this section uses bilinear interpolation to map
c the four surrounding input pixels to the output grid.
c
          write(6,*)'image res .ge. output grid spacing'
          write(6,*)'using bilinear interp for sndr rad'
          do j=1,jmax
          do i=1,imax
c            if(i_found_it .ne. 1)then
c              write(6,*)'enter i and j point'
c              read(5,*)ipt,jpt
c              i_found_it = 1
c            end if
c
            if(r_llij_lut_ri(i,j).ne.r_missing_data .and.
     &         r_llij_lut_rj(i,j).ne.r_missing_data)then

               call bilinear_laps(r_llij_lut_ri(i,j),r_llij_lut_rj(i,j),
     &       elem_dim,line_dim,image_sndr,result,istatus)

              sa(i,j) = result
              if(result.le.0)then
                 print*,'<= 0.0: sndr values'
                 i1=int(r_llij_lut_ri(i,j))
                 i2=int(r_llij_lut_ri(i+1,j))
                 j1=int(r_llij_lut_rj(i,j))
                 j2=int(r_llij_lut_rj(i,j+1))
           print*,'sndr(i,j)/sndr(i+1,j)/sndr(i,j+1)/sndr(i+1,j+1):'
                 print*,image_sndr(i1,j1),image_sndr(i1,j2)
     &,image_sndr(i2,j1),image_sndr(i2,j2) 

              endif
            else
                sa(i,j) = r_missing_data
                icnt_out=icnt_out+1
            endif

            sc(i,j)=sa(i,j)
            st(i,j)=sa(i,j)

c            if(i.eq.ipt .and. j.eq.jpt)then
c               write(29,39)ipt,jpt
c               write(29,49)image_sndr(ii,jj)
c               i_found_it = -1
c 39          format(1x,'ipt = ',i3,2x,'jpt = ',i3)
c 49      format(1x,'ir_count = ',1x,f6.1,/,'--------------------',//)
c            end if

           enddo
           enddo
c          close(29)

        end if     ! r_image_ratio

        if(icnt_out.gt.0)then
           print*,'found ', icnt_out,' gdpts out of domain'
        endif

        icnt_out=0
        do j=1,jmax
        do i=1,imax
           if(sa(i,j).le.0.0)then
              icnt_out=icnt_out+1
           endif
        enddo
        enddo
        if(icnt_out.gt.0)print*,'found ',icnt_out,' <= 0.0 pts'

c       write(6,1234) ib,i4vtime,ict
c1234       format(1x,'band ',i4,' count for i4time ',i10,' is ',i8)

        istatus = 1
c
        return
        end

c ------------------

        subroutine initialize(a,imax,jmax,r_missing_data)

        implicit none

        integer  imax,jmax
        real     a(imax,jmax)
        real     r_missing_data

        integer  i,j

        do j=1,jmax
        do i=1,imax

           a(i,j)=r_missing_data

        enddo
        enddo
        return
        end
