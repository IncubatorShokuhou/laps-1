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
cdis
cdis   
cdis
      subroutine ccpfil(field_in,mreg,nreg,scale_l_in,scale_h_in
     1                 ,colortable,n_image,scale,c5_sect
     1                 ,plot_parms,namelist_parms)       


      include 'lapsplot.inc'
c         

c define error file, fortran unit number, and workstation type,
c and workstation id.
c 
      parameter (ierrf=6, lunit=2, iwtype=1, iwkid=1)
      real xreg(mreg),yreg(nreg),zreg(mreg,nreg),field_in(mreg,nreg)
      character*(*)colortable
      character*5 c5_sect

      logical log_scaling, l_integral, l_discrete, l_divisible
      logical l_set_contours

      common /plot_field_cmn/ i_plotted_field

      i_plotted_field = 1

      write(6,*)' subroutine ccpfil for solid fill plot...'

      if(colortable(1:3) .eq. 'acc' .or. colortable(1:3) .eq. 'sno')then       
          l_set_contours = .true.             ! for testing
!         l_set_contours = .false.            ! operational
      else
          l_set_contours = .false.
      endif

      if( (colortable(1:3) .eq. 'acc' .or. colortable(1:3) .eq. 'sno')
     1                        .and. (.not. l_set_contours)         )then       
          log_scaling = .true.
      else
          log_scaling = .false.
      endif

      if(colortable .eq. 'haines' .or. colortable .eq. 'cwi')then
          l_integral = .true.
      else
          l_integral = .false.
      endif

!     if(colortable .eq. 'tpw')then
!         power = 0.7
!     else
          power = plot_parms%color_power
!     endif

      write(6,*)' colortable power = ',power

      n_image = n_image + 1

      call get_r_missing_data(r_missing_data,istatus)

      call array_range(field_in,mreg,nreg,rmin,rmax,r_missing_data)
      write(6,*)' data array range (before scaling): ',rmin,rmax
      write(6,*)' input scale: ',scale
      write(6,*)' scaled range: ',rmin/scale,rmax/scale
     1                                      ,(rmax/scale)**power

!     the 0.05 empirically relates to the inverse number of colors
!     as plots apparently must span a full color interval
      if(rmin .eq. 0. .and. (rmax/scale)**power .le. 0.05 
     1                .and. c5_sect .eq. 'xsect')then
          write(6,*)' note, xsect background may be non-black - return'
          return
      endif

      write(6,*)' scale_l_in,scale_h_in = ',scale_l_in,scale_h_in

!     if(n_image .gt. 1)then
!         write(6,*)' image was already plotted - returning from ccpfil'
!         return
!     endif

      if(scale_l_in .lt. scale_h_in)then
          ireverse = 0
          scale_l = scale_l_in * scale
          scale_h = scale_h_in * scale
      else
          ireverse = 1
          scale_l = scale_h_in * scale
          scale_h = scale_l_in * scale
      endif

      if(l_integral)then
          scale_l = scale_l - 0.5 * scale
          scale_h = scale_h + 0.5 * scale
      endif

!     apply scaling to the array
!     call addcon(field_in,-scale_l,zreg,mreg,nreg)

      if(log_scaling)then ! e.g. for precip
          scale_l = .01 * scale

          do i = 1,mreg
          do j = 1,nreg
            zreg(i,j) = alog10(max(field_in(i,j),scale_l))
          enddo ! j
          enddo ! i

          scale_l = alog10(scale_l)
          scale_h = alog10(scale_h)

          zreg = zreg - scale_l                          ! array subtraction

          scale_loc = scale_h - scale_l

      elseif(.not. l_set_contours)then
          scale_loc = (scale_h - scale_l)**power

          zreg = field_in - scale_l                      ! array subtraction

!         adjust field values (missing data or reverse cases or power law)
          do i = 1,mreg
          do j = 1,nreg
            if(field_in(i,j) .eq. r_missing_data)then

              if(c5_sect .eq. 'hsect')then

                if(colortable(1:3) .eq. 'lin')then
                  zreg(i,j) = scale_loc * 0.50 ! e.g. csc 

                elseif(colortable(1:3) .eq. 'cpe')then ! apply only to hsects
                  if(scale_h_in .eq. 7000.)then 
                      zreg(i,j) = scale_loc * 0.00 ! ! cape for hsect
                  elseif(scale_h_in .eq. 50.)then 
                      zreg(i,j) = scale_loc * 1.00 ! ! cin for hsect
                  endif

                elseif(l_integral)then ! e.g. haines
                  zreg(i,j) = scale_loc * 100. ! 1.2 

                else

                  zreg(i,j) = scale_loc * 0.96 
                endif

              endif

            else ! valid data 
              if(ireverse .eq. 1)then
                zreg(i,j) = scale_loc - (zreg(i,j)**power)
              else                       
                zreg(i,j) = zreg(i,j)**power
              endif

            endif

!           prevent overshoot beyond colortable (except for discrete)
            if(c5_sect .eq. 'hsect')then
              if(zreg(i,j) .gt. scale_loc     .and. 
!    1             colortable(1:3) .ne. 'cpe' .and.
     1             (.not. l_integral)               )then
                zreg(i,j) = scale_loc
              endif
            endif

            if(c5_sect .eq. 'xsect')then
              if(zreg(i,j) .gt. scale_loc      .and. 
     1           zreg(i,j) .ne. r_missing_data .and.
     1             (.not. l_integral)               )then
                zreg(i,j) = scale_loc
              endif
            endif

            if(zreg(i,j) .lt. 0.0      )then
              zreg(i,j) = 0.0
            endif

          enddo ! j
          enddo ! i

      else ! l_set_contours = .true.
          zreg = field_in 

      endif ! type of scaling / contours

      write(6,*)' scale_h/scale_l/scale_loc =',scale_h,scale_l,scale_loc      

      ireverse_colorbar = ireverse

      ireverse = 0  ! turn off later use of ireverse
c      
c get data array
c
!     call getdat(xreg,yreg,zreg,mreg,nreg)
c 
c open gks, open and activate a workstation.
c 
!     call gopks (ierrf, iszdm)
!     call gopwk (iwkid, lunit, iwtype)
!     call gacwk (iwkid)
c      
c call conpack color fill routine
c      
      icol_offset = 40 ! offset new colortable to preserve previous low end
      ncols = 20

      range = scale_loc/scale
      write(6,*)' scale_loc/scale/range = ',scale_loc,scale,range
      call get_colorbar_int(range,colorbar_int,l_divisible)

      if(log_scaling)then
          l_discrete = .false.
      else
          if(colortable .eq. 'vnt')then
              l_discrete = .true.
          elseif(colortable .eq. 'spectral' .or. 
     1           colortable .eq. 'spectralr')then
!             test for isotachs 
!             if(scale_h_in .eq. 50. .and. scale_l_in .eq. 0.)then 
!                 l_discrete = .true. 
!             else
                  l_discrete = plot_parms%l_discrete
!             endif
          elseif(colortable .eq. 'moist')then
              l_discrete = plot_parms%l_discrete
          elseif(colortable .eq. 'tpw')then
              l_discrete = plot_parms%l_discrete
          elseif(colortable .eq. 'upflux')then
              l_discrete = plot_parms%l_discrete
          elseif(colortable .eq. 'hues'  .or. 
     1           colortable .eq. 'omega' .or.
     1           colortable .eq. 'temp'
     1                                       )then
              if(l_divisible)then
                  l_discrete = plot_parms%l_discrete
              else
                  l_discrete = .false.
              endif
          elseif(colortable .eq. 'cpe')then       
              if(l_divisible)then
                  l_discrete = plot_parms%l_discrete
              else
                  l_discrete = .false.
              endif
          else
              l_discrete = .false.
          endif
      endif

      write(6,*)' colortable is ',colortable,scale_l,scale_h
     1                           ,ireverse_colorbar,l_discrete
     1                           ,log_scaling,l_integral,l_set_contours       

      if(l_discrete)then
!         set interval for writing numbers
          if(l_integral)then
              colorbar_int = 0.5
          else
              continue ! use pre-existing value of colorbar_int
          endif

          ncols = nint(range / colorbar_int)

          write(6,*)' l_discrete case: range/colorbar_int/ncols = '
     1             ,range,colorbar_int,ncols  
      endif

      lmap=mreg*nreg*256 ! 16000000
      lmap = min(lmap,32000000)
      call ccpfil_sub(zreg,mreg,nreg,-15,iwkid,scale_loc,scale,ireverse       
     1                               ,r_missing_data,plot_parms        ! i
     1                               ,namelist_parms                   ! i
     1                               ,lmap,log_scaling,l_set_contours  ! i
     1                               ,colortable                       ! i
     1                               ,ncols,power                      ! i/o
     1                               ,l_discrete,c5_sect               ! i
     1                               ,icol_offset)      
c      
c close frame
c      
!     call frame
c 
c deactivate and close workstation, close gks.
c 
!     call gdawk (iwkid)
!     call gclwk (iwkid)
!     call gclks

c     call local colorbar routine
      write(6,*)' drawing colorbar: ',mreg,nreg
      call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
      call colorbar(mreg, nreg, namelist_parms, plot_parms, 
     1              ncols, ireverse_colorbar, log_scaling, power,      ! i
     1              scale_l, scale_h, colortable, scale,icol_offset,
     1              c5_sect, l_discrete, l_integral, l_set_contours,
     1              colorbar_int)

      jdot = 1
      
      return
      end

      
      subroutine ccpfil_sub(zreg,mreg,nreg,ncl,iwkid,scale_loc,scale
     1                                ,ireverse       
     1                                ,r_missing_data,plot_parms       ! i
     1                                ,namelist_parms                  ! i
     1                                ,lmap,log_scaling,l_set_contours
     1                                ,colortable                      ! i
     1                                ,ncols,power                     ! i/o
     1                                ,l_discrete,c5_sect              ! i
     1                                ,icol_offset)      

      include 'lapsplot.inc'

      parameter (lrwk=600000,liwk=600000,nwrk=600000
     1          ,nogrps=5)       
      real zreg(mreg,nreg),rwrk(lrwk), xwrk(nwrk), ywrk(nwrk)
      integer mreg,nreg,iwrk(liwk)
      integer map(lmap),iarea(nogrps),igrp(nogrps)
      integer colia(mreg,nreg)
      character*(*) colortable
      character*5 c5_sect
      logical log_scaling, l_discrete, l_set_contours, l_raster
      logical l_uniform_col

      integer maxvals
      parameter (maxvals=100)

      real vals(maxvals),vals_scaled(maxvals)
      
      external fill
c      
c set up color table
      write(6,*)' ccpfil_sub - scale,scale_loc = ',scale,scale_loc
c      
      if(l_set_contours)then
          if(colortable(1:3) .eq. 'acc')then       
              if(colortable(4:7) .eq. '_inc')then
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals,1)
              else
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals
     1                             ,namelist_parms%i_pcp_sto_colorbar)
              endif
          endif

          if(colortable(1:3) .eq. 'sno')then       
              if(colortable(4:7) .eq. '_inc')then
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals,2)
              else
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals
     1                             ,namelist_parms%i_sno_sto_colorbar)
              endif
          endif

          ncols = nvals - 1
      endif

      call set_image_colortable(iwkid
     1                         ,plot_parms                           ! i
     1                         ,ncols                                ! i/o
     1                         ,l_discrete,ireverse,c5_sect
     1                         ,l_set_contours,colortable,power
     1                         ,mreg,nreg,log_scaling,icol_offset)
c      
c initialize areas
c      
      call arinam(map, lmap)

      if(l_set_contours)then
          icol_offset2 = icol_offset - 2

!         set up 'clv' array
          call cpseti('ncl',nvals+icol_offset2)

          do i = 1,icol_offset2 ! bogus in contour values in offset region
              call cpseti ('pai - parameter array index',i)
              call cpsetr ('clv - contour level'
     1                    ,vals_scaled(1)-float(i))      
          enddo ! i

          do i = 1,nvals
              is = i + icol_offset2
              vals_scaled(is) = vals(i) * scale
              write(6,*)'is/scaled value',is,vals_scaled(is)
              call cpseti ('pai - parameter array index',is)
              call cpsetr ('clv - contour level',vals_scaled(is))
          enddo ! i

          call cpseti('cls - contour level selection flag',0)

          l_uniform_col = .false.

      else
          col_offset = float(icol_offset) / float(ncols)
          write(6,*)' col_offset / scale_loc = ',col_offset,scale_loc       
          cis = abs(scale_loc) / float(ncols)

          icol_min = +1000000
          icol_max = -1000000

          do m = 1,mreg
          do n = 1,nreg
              zreg(m,n) = zreg(m,n) + (col_offset * scale_loc)

              colia_arg = (zreg(m,n)-cmn)/cis
              if(abs(colia_arg) .le. 1e9)then
                  colia(m,n) = int(colia_arg) + 3
              else
                  colia(m,n) = 0
              endif

!             colia(m,n) = int( (zreg(m,n)-cmn)/cis ) + 3

              icol_min = min(colia(m,n),icol_min)
              icol_max = max(colia(m,n),icol_max)
          enddo ! n
          enddo ! m

          write(6,*)' zreg range ',minval(zreg),maxval(zreg)
          write(6,*)' colia range (1) ',icol_min,icol_max

          if(icol_min .eq. icol_max)then
              l_uniform_col = .true.
          else
              l_uniform_col = .false.
          endif
c      
c         set number of contour levels and initialize conpack
c      
!         call cpseti('cls - contour level selection flag',ncl)

          call cpseti('cls - contour level selection flag',+1)
          call cpsetr('cis', cis)
          call cpsetr('cmn',(0.0           ) * abs(scale_loc) + 2.0*cis)
          call cpsetr('cmx',(1.0+col_offset) * abs(scale_loc) + 2.0*cis)       

      endif

      call cpsetr ('spv',r_missing_data)

!     nreg threshold =    0 means do raster image whenever possible
!     nreg threshold =  500 means do raster image for large domains (e.g. stmas)
!     nreg threshold = 9999 means never do raster image 
      if(nreg .ge. 300)then ! set raster default
          l_raster = .true.
          write(6,*)' large domain - turn on raster plot'
      else
          l_raster = .false.
      endif

      if(l_uniform_col)then
          write(6,*)' uniform colors detected...'
          l_raster = .true.
      endif

!     override default with parameter input
      if(plot_parms%iraster .eq. 1)then
          l_raster = .true.
          write(6,*)' turn on raster plot for this field'
      elseif(plot_parms%iraster .eq. -1)then
          l_raster = .false.
      endif

      if(plot_parms%zoom_wdw .gt. 1.0)then
          l_raster = .false.
          write(6,*)' turn off raster plot for zoomed field'
      endif       

      if(l_raster .and. (.not. l_set_contours))then ! do raster image

!         try and reset things borrowing from the regular image section
          call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)

          cmn = (0.0           ) * abs(scale_loc) + 2.0*cis
          cmx = (1.0+col_offset) * abs(scale_loc) + 2.0*cis       
          crange = cmx-cmn
          do m = 1,mreg
          do n = 1,nreg
              colia_arg = (zreg(m,n)-cmn)/cis
              if(abs(colia_arg) .le. 1e9)then
                  colia(m,n) = int(colia_arg) + 3
              else
                  colia(m,n) = 0
              endif
!             write(6,*)m,n,colia(m,n)
          enddo ! n
          enddo ! m

          write(6,*)' zreg range (2)',minval(zreg),maxval(zreg)
          write(6,*)' colia range (2) ',icol_min,icol_max

          call get_r_missing_data(r_missing_data,istatus)

          if(c5_sect .eq. 'xsect')then
              where(zreg(:,:) .eq. r_missing_data)
     1              colia(:,:) = 0
          endif

          write(6,*)' colia range (3) ' ! ,minval(colia),maxval(colia)

          write(6,*)' calling gca for raster fill plot'
          call get_border(mreg,nreg,x_1,x_2,y_1,y_2)
          call gca (x_1, y_1, x_2, y_2, mreg, nreg, 1,  1,
     1              mreg, nreg, colia)

      else                                           ! regular contour fill
          write(6,*)' regular image contour fill without raster ' 
     1             ,l_raster,l_set_contours

          call cprect(zreg, mreg, mreg, nreg, rwrk, lrwk, iwrk, liwk)
c      
c add contours to area map
c      
          call cpclam(zreg, rwrk, iwrk, map)
c      
c set fill style to solid, and fill contours
c      
          call gsfais(1)
          call arscam(map, xwrk, ywrk, nwrk, iarea, igrp, nogrps, fill)       

      endif
c      
c draw perimeter
c      
      call setusv_dum(2hin,34) ! gray
      call perim(1,1,1,1)
c      
c draw labels
c      
!     call cplbdr(zreg,rwrk,iwrk)
c      
c draw contours
c      
!     call cpcldr(zreg,rwrk,iwrk)
      
      return
      end
      
      subroutine set_image_colortable(iwkid,plot_parms,ncols
     1                               ,l_discrete,ireverse,c5_sect
     1                               ,l_set_contours,colortable,power
     1                               ,mreg,nreg,log_scaling,icol_offset)    

      include 'lapsplot.inc'

      character*(*) colortable
      character*5 c5_sect
      logical log_scaling, l_discrete, l_set_contours
c 
c background color
c black
c
!     call gscr(iwkid,0,0.,0.,0.)

      write(6,*)' subroutine set_image_colortable: ncols/mreg/nreg = '
     1         ,ncols,mreg,nreg

      if(colortable(1:3) .eq. 'lin')then
          if(colortable .eq. 'linear_reduced')then
              ncols = 10
          else
              ncols = 40
          endif

          rcols = ncols - 1
          do i = 1,ncols+10 ! 255-icol_offset
              if(colortable .eq. 'linear_reduced')then
!                 rintens = min(max( (float(i) / rcols) - 0.4 ,0.),1.)       
                  rintens = min(max(float(i-1) / rcols,0.),1.)
              else
                  rintens = min(max(float(i-1) / rcols,0.),1.)
              endif

!             rintens = 0.08 + (rintens * 0.92) ! brighten on the monitor
!             rintens = (rintens**0.76) * 0.95  ! brighten on the monitor
              rintens = (rintens**0.40) * 0.95  ! brighten on the monitor

              if(ireverse .eq. 1)rintens = 1.0 - rintens
              call gscr(iwkid, i+icol_offset, rintens, rintens, rintens)
          enddo ! i

      elseif(colortable .eq. 'hues' .or. colortable .eq. 'ref')then       

          if(.not. l_discrete)then
              ncols = 60
          endif

          call generate_colortable(ncols,'hues',iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

          if(colortable .eq. 'ref')then
              do i = 1,3
                  call gscr(iwkid, i+icol_offset, 0., 0., 0.)
              enddo 

              do i = ncols,ncols
                  call gscr(iwkid, i+icol_offset, 0.3, 0.3, 0.3)
              enddo
          endif

      elseif(colortable .eq. 'cpe')then       

          if(.not. l_discrete)then
              if(mreg*nreg .gt. 62500)then       
                  ncols = 24
              else
                  ncols = 48
              endif
          endif

          call generate_colortable(ncols,'cpe',iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

          do i = 1,1
              call gscr(iwkid, i+icol_offset, 0., 0., 0.)
          enddo 

!         if(.not. l_discrete)then
!             if(mreg*nreg .gt. 62500)then       
!                 ncols = 24
!             else
!                 ncols = 48
!             endif
!         endif

      elseif(colortable .eq. 'tpw' .or. colortable .eq. 'upflux')then       
          if(.not. l_discrete)then
              ncols = 179 ! previous good value was 79 (139)      
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(colortable .eq. 'moist')then
          if(.not. l_discrete)then
              ncols = 60
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(            colortable .eq. 'spectral' 
     1              .or. colortable .eq. 'spectralr'      )then       

          if(.not. l_set_contours)then
              if(.not. l_discrete)then
                  if(            colortable .eq. 'spectralr'     
!    1                      .or. mreg*nreg .gt.  300000)then       
     1                      .or. mreg*nreg .gt. 2000000)then       
                      ncols = 20
                  else
                      ncols = 60
                  endif
              endif
          endif

          call generate_colortable(ncols,'spectral',iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(colortable .eq. 'haines')then       
          ncols = 5 
          call color_ramp(1,1
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,0.6,0.7,0.4                 ! violet
     1                   ,0.6,0.7,0.4)                ! violet
          call color_ramp(2,2
     1                   ,iwkid,icol_offset,plot_parms       
     1                   ,1.0,0.85,0.55               ! blue
     1                   ,1.0,0.85,0.55)              ! blue
          call color_ramp(3,3
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,2.0,0.4,0.4                 ! green
     1                   ,2.0,0.4,0.4)                ! green
          call color_ramp(4,4
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,2.5,0.95,0.60               ! yellow
     1                   ,2.5,0.95,0.60)              ! yellow
          call color_ramp(5,5
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,3.0,0.9,0.7                 ! red
     1                   ,3.0,0.9,0.7)                ! red

!         extra color at the top for r_missing_data
          call color_ramp(6,6
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,1.0,0.0,1.0                 ! white
     1                   ,1.0,0.0,1.0)                ! white

      elseif(colortable .eq. 'cwi')then       
          ncols = 2 
          call color_ramp(1,1
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,0.6,0.7,0.4                 ! violet
     1                   ,0.6,0.7,0.4)                ! violet
          call color_ramp(2,2
     1                   ,iwkid,icol_offset,plot_parms
     1                   ,3.0,0.9,0.7                 ! red
     1                   ,3.0,0.9,0.7)                ! red

      elseif(colortable .eq. 'omega')then       
          if(l_discrete)then
              ncols = 19
          elseif(mreg*nreg .ge. 62500 .and. c5_sect .eq. 'hsect')then
              ncols = 19
          else
              ncols = 79
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(colortable .eq. 'vnt')then       
          if(.not. l_discrete)then
              ncols = 60
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(colortable .eq. 'temp')then       
          if(.not. l_discrete)then
              ncols = 76
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      elseif(colortable .eq. 'green' .or. colortable(1:3) .eq. 'acc'
     1  .or. colortable(1:3) .eq. 'sno'                            )then       
          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      else ! generic color table settings
          if(.not. l_discrete)then
              ncols = 179 ! previous good value was 79 (139)      
          endif

          call generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                            ,power,plot_parms,istatus)

      endif

      if(colortable(1:3) .eq. 'acc' .or. 
     1   colortable(1:3) .eq. 'sno')then ! set colortable ends
          do i = 1,1
!         do i = 1,3
              call gscr(iwkid, i+icol_offset, 0., 0., 0.)
          enddo 

!         do i = ncols,ncols
!             call gscr(iwkid, i+icol_offset, 0.3, 0.3, 0.3)
!         enddo
      endif
c 
      write(6,*)' leaving set_image_colortable: ncols = ',ncols

      return
      end

      subroutine color_ramp(ncol1,ncol2,iwkid,icol_offset       ! i
     1                     ,plot_parms
     1                     ,hue1,sat1,rintens1                  ! i
     1                     ,hue2,sat2,rintens2)                 ! i

      include 'lapsplot.inc'

      write(6,*)'  subroutine color_ramp.. ',ncol1,ncol2

      do icol = ncol1,ncol2
          if(ncol2 .ne. ncol1)then
              frac = float(icol-ncol1) / float(ncol2-ncol1)
          else
              frac = 0.0
          endif   

          hue     = (1.0 - frac) * hue1     + frac * hue2  
          sat     = (1.0 - frac) * sat1     + frac * sat2
          rintens = (1.0 - frac) * rintens1 + frac * rintens2

!         modify intensity according to user input
          rintens = rintens * plot_parms%rimage_intensity

          call hsi_to_rgb(hue,sat,rintens,red,grn,blu)

          write(6,1)icol,frac,hue,sat,rintens,red,grn,blu
 1        format(i5,f8.3,3x,'hsi:',3f8.3,3x,'rgb:',3f8.3)
          
          call gscr(iwkid,icol+icol_offset,red,grn,blu)
      enddo

      if(ncol2+2+icol_offset .le. 255)then
          call gscr(iwkid,ncol2+1+icol_offset,red,grn,blu)
          call gscr(iwkid,ncol2+2+icol_offset,red,grn,blu)
      endif

      return
      end

      subroutine hsi_to_rgb(hue,sat,rintens,red,grn,blu)

!     hue is 0:r, 1:b, 2:g, 3:r

      red1 = max(1.0 - abs(hue - 0.0),0.0)
      red2 = max(1.0 - abs(hue - 3.0),0.0)
      red = max(red1,red2)
      grn = max(1.0 - abs(hue  - 2.0),0.0)
      blu = max(1.0 - abs(hue  - 1.0),0.0)

!     normalize to the max intensity
      colmax = max(red,grn,blu)
      if(colmax .gt. 0.)then
          red = red/colmax
          grn = grn/colmax
          blu = blu/colmax
      endif

      red = (red*sat) + 1.0*(1.0-sat)
      grn = (grn*sat) + 1.0*(1.0-sat)
      blu = (blu*sat) + 1.0*(1.0-sat)

      red = red * rintens
      grn = grn * rintens
      blu = blu * rintens

      return
      end



      subroutine colorbar(ni,nj,namelist_parms,plot_parms 
     1                   ,ncols,ireverse,log_scaling,power
     1                   ,scale_l,scale_h
     1                   ,colortable,scale,icol_offset,c5_sect
     1                   ,l_discrete,l_integral,l_set_contours          ! i
     1                   ,colorbar_int)                                 ! i

      include 'lapsplot.inc'

      character*8 ch_low, ch_high, ch_mid, ch_frac
      character*(*)colortable
      character*5 c5_sect
      logical log_scaling,l_loop, l_discrete, l_integral, l_divisible       
      logical l_set_contours

      integer maxvals
      parameter (maxvals=100)

      real frac_a(maxvals)
      real vals(maxvals)

      write(6,*)' colorbar: scale_l,scale_h,scale',scale_l,scale_h,scale

      range = abs(scale_h - scale_l) / scale

      call get_colorbar_int(range,colorbar_int,l_divisible)

      if(scale_l .eq. -20. .or. colortable .eq. 'cpe' 
     1                     .or. nint(range) .eq. 8000    ! pbl
     1                     .or. l_integral               ! e.g. hah, ham, cwi
     1                     .or. l_discrete               
     1                     .or. colortable .eq. 'spectral' 
     1                     .or. colortable .eq. 'spectralr' 
     1                     .or. colortable .eq. 'tpw'    ! tpw
     1                     .or. colortable .eq. 'upflux' ! upslope moist flux
     1                     .or. colortable .eq. 'vnt'    ! vnt
     1                     .or. l_divisible              ! divisible colorbars
     1                     .or. colortable .eq. 'temp'
     1                     .or. range .eq. 100.    )then ! sfc t, td, rh, cape
          if(l_set_contours)then
              l_loop = .false.
          else
              l_loop = .true.
          endif
      else
          l_loop = .false.
      endif

      write(6,*)' l_integral,l_discrete,l_divisible,l_set_contours',
     1            l_integral,l_discrete,l_divisible,l_set_contours

!     if(colortable(1:3) .eq. 'lin')l_loop = .false.

      call get_border(ni,nj,x_1,x_2,y_1,y_2)

      xlow =  0.40 ! 0.35
      xhigh = xlow + 0.50

      if(c5_sect .eq. 'xsect')then
          y_2 = 0.81
      endif

      ylow =  y_2 + .01
      yhigh = y_2 + .03

!     set for zoom
      zoom_colorbar = plot_parms%zoom_wdw

!     note that "square" case works for aspect ratio up to 1.192
      write(6,*)'colorbar zoom/power = ',zoom_colorbar,power

      write(6,*)'colorbar location before zoom ',xlow,xhigh,ylow,yhigh

      call zoomit(plot_parms,zoom_colorbar,xlow,ylow,xlowz,ylowz)
      call zoomit(plot_parms,zoom_colorbar,xhigh,yhigh,xhighz,yhighz)

      write(6,*)'colorbar location after zoom ',
     1          xlowz,xhighz,ylowz,yhighz

!     set up for number of lines in colorbar
      ilow = 1
      ihigh = 2000

      xrange  = xhigh  - xlow 
      xrangez = xhighz - xlowz

      irange = ihigh - ilow

!     put colorbar
      do i = ilow,ihigh
          frac = float(i-ilow) / float(irange)
          x1   = xlowz + frac*xrangez 
          x2   = xlowz + frac*xrangez 

          if(ireverse .eq. 0)then
              rcol = 0.5 + float(ncols) * (frac**power)
          else
              rcol = 0.5 + float(ncols) * ((1.0 - frac)**power)
          endif

          icol = nint(rcol)

          call setusv_dum(2hin,icol+icol_offset)

          y1 = ylowz
          y2 = yhighz
          call line(x1,y1,x2,y2)
      enddo ! i

c     restore original color table
!     call color

!     write labels at middle and ends of colorbar
      call setusv_dum(2hin,34) ! gray

      call line(xlowz,ylowz,xhighz,ylowz)
      call line(xlowz,yhighz,xhighz,yhighz)
      call line(xlowz,ylowz,xlowz,yhighz)
      call line(xhighz,ylowz,xhighz,yhighz)

      call setusv_dum(2hin,7)  ! yellow

      rsize = .0072 ! .008
      iy = (y_2+.021) * 1024

      if(.not. l_integral .and. .not. l_set_contours)then

!         left edge
          if(log_scaling)then
              rlow = 0.
          else
              rlow = scale_l/scale
          endif

          if(abs(rlow) .gt. 0.0 .and. 
     1       abs(rlow) .le. 0.5                  )then
              write(ch_low,3)rlow
              call right_justify(ch_low)
          else
              write(ch_low, 1)nint(rlow)
              call right_justify(ch_low)
          endif

!         ixl = 353 + nint(.05 * 1024.)
          ixl = nint((xlow - .005) * 1024.)
!         write(6,*)' cpux = ',xlow,ixl,iy,cpux(ixl),cpux(iy)
          call zoomit(plot_parms,zoom_colorbar
     1               ,cpux(ixl),cpux(iy),rxl,ry)
          call pchiqu (rxl,ry,ch_low,rsize/zoom_colorbar ,0,+1.0)
          write(6,*)' colorbar left edge = ',ch_low

!         right edge
          if(log_scaling)then
              rhigh = (10.**scale_h) / scale
          else
              rhigh = scale_h / scale
          endif

          if(abs(rhigh) .ge. 1.0)then
              if(abs(rhigh-float(nint(rhigh))) .lt. .05 .or. 
     1                                     colortable .eq. 'haines')then       
                  write(ch_high,1)nint(rhigh)
 1                format(i8)
              else
                  write(ch_high,2)rhigh
              endif
          else
              write(ch_high,3)rhigh
 3            format(f8.2)
          endif
          call left_justify(ch_high)

          ixh = ixl + 525 ! 878
          call zoomit(plot_parms,zoom_colorbar
     1               ,cpux(ixh),cpux(iy),rxh,ry)
          call pchiqu (rxh,ry,ch_high,rsize/zoom_colorbar,0,-1.0)
          write(6,*)' colorbar right edge = ',ch_high

      endif ! l_integral

      if(.not. l_loop .and. .not. l_integral 
     1                .and. .not. l_set_contours
     1                .and. .not. log_scaling)then ! plot midpoint

          frac = 0.5
          x1   = xlow + frac*xrange 
          x2   = xlow + frac*xrange 
          call setusv_dum(2hin,0)

          y1 = ylow
          y2 = yhigh
          call line(x1,y1,x2,y2)

!         plot number
          call setusv_dum(2hin,7)  ! yellow
          if(log_scaling)then
              rmid = (10.** ((scale_l+scale_h) / 2.0) ) / scale
          else
              rmid = ((scale_l+scale_h) / scale)/2.0
          endif

          if( (abs(rmid) .gt. 1.0 .or. abs(rlow) .gt. 1.0
     1                            .or. abs(rhigh) .gt. 1.0 )
     1                            .and. 
     1                abs(rmid-float(nint(rmid))) .lt. .05   
     1                                                           )then       
              write(ch_mid,1)nint(rmid)
          elseif(abs(rhigh) .ge. 1.0)then
              write(ch_mid,2)rmid
2             format(f8.1)
          else
              write(ch_mid,3)rmid
          endif 
          call left_justify(ch_mid)
          call s_len(ch_mid,len_mid)

          ixm = (ixl+ixh)/2
          call zoomit(plot_parms,zoom_colorbar
     1               ,cpux(ixm),cpux(iy),rxm,ry)
          call pchiqu (rxm,ry,ch_mid(1:len_mid),rsize/zoom_colorbar
     1                ,0,0.0)       
          write(6,*)' colorbar midpoint = ',ch_mid(1:len_mid)

      endif

      ixl = 409
      ixh = 924

      if(l_loop)then ! plot additional numbers

!         set interval for writing numbers
          if(l_discrete)then
              continue ! use passed in value of colorbar_int
          elseif(l_integral)then
              colorbar_int = 0.5
          else
              call get_colorbar_int(range,colorbar_int,l_divisible)
          endif

!         interval for writing lines
          colorbar_int = colorbar_int * scale / 1.0

          write(6,*)' plotting colorbar',scale_l,colorbar_int,ixl,ixh          

          loop_count = 0

          do rarg = scale_l+colorbar_int,scale_h-.001*scale,colorbar_int
              frac = (rarg - scale_l) / (scale_h - scale_l)

              loop_count = loop_count + 1

              if(.not. l_integral .and. .not. l_discrete)then ! plot black line
                  x1   = xlow + frac*xrange 
                  x2   = xlow + frac*xrange 
                  call setusv_dum(2hin,25)                    ! dark gray

                  y1 = ylow
                  y2 = yhigh
                  call zoomit(plot_parms,zoom_colorbar,x1,y1,x1z,y1z)
                  call zoomit(plot_parms,zoom_colorbar,x2,y2,x2z,y2z)
                  call line(x1z,y1z,x2z,y2z)
              endif

              if(loop_count .eq. (loop_count/1) * 1 )then ! number every line
!                 plot number
                  call setusv_dum(2hin,7)  ! yellow
                  if(log_scaling)then
                      rlabel = (10.**(rarg)) / scale
                  else
                      rlabel = rarg / scale
                  endif

                  if(range .lt. 0.2)then
                      write(ch_frac,3)rlabel
                  elseif(range .lt. 2.0)then
                      write(ch_frac,2)rlabel
                  else
                      if( abs(rlabel- nint(rlabel)) .lt. .001)then
                          write(ch_frac,1)nint(rlabel)
                      else
                          write(ch_frac,2)rlabel
                      endif
                  endif

                  call left_justify(ch_frac)
                  call s_len(ch_frac,len_frac)

                  ixm = ixl + (ixh-ixl)*frac

                  call zoomit(plot_parms,zoom_colorbar
     1                       ,cpux(ixm),cpux(iy),rxm,ry)
                  if(l_integral)then
                      if(rlabel .eq. float(nint(rlabel)))then
                          call pchiqu (rxm,ry,ch_frac(1:len_frac)       
     1                                ,rsize/zoom_colorbar,0 , 0.0)    
                      endif
                  else
                      call pchiqu (rxm,ry,ch_frac(1:len_frac)          
     1                            ,rsize/zoom_colorbar,0 , 0.0)       
                  endif

              endif ! loop_count
          enddo ! rarg

      else      
!         other fractions

          if(l_set_contours)then
              if(colortable(1:3) .eq. 'acc')then
                if(colortable(4:7) .eq. '_inc')then
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals,1)
                else
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals
     1                             ,namelist_parms%i_pcp_sto_colorbar)
                endif
              elseif(colortable(1:3) .eq. 'sno')then
                if(colortable(4:7) .eq. '_inc')then
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals,2)
                else
                  call get_pcp_vals(maxvals,namelist_parms,nvals,vals
     1                             ,namelist_parms%i_sno_sto_colorbar)
                endif
              endif

              do i = 1,nvals
                  frac_a(i) = float(i-1) / float(nvals-1)
              enddo ! i
              nfrac = nvals

          elseif(log_scaling .and. 
     1          (colortable(1:3) .eq. 'acc' .or. 
     1           colortable(1:3) .eq. 'sno'     )  )then
              frac_a(1) = 0.125
              frac_a(2) = 0.230
              frac_a(3) = 0.330
              frac_a(4) = 0.436  ! .420 = .18
              frac_a(5) = 0.490  ! .500 = .32
              frac_a(6) = 0.567
              frac_a(7) = 0.670
              frac_a(8) = 0.77
              frac_a(9) = 0.824  ! .830 = 3.1
              frac_a(10) = 0.900
              nfrac = 10

          else
              frac_a(1) = 0.25
              frac_a(2) = 0.75
              nfrac = 2

          endif

          do ifrac = 1,nfrac
              frac = frac_a(ifrac)

              if(frac .ne. 0.0 .and. frac .ne. 1.0)then
!                 plot black line
                  x1   = xlow + frac*xrange 
                  x2   = xlow + frac*xrange 
                  call setusv_dum(2hin,0)

                  y1 = ylow
                  y2 = yhigh
                  call line(x1,y1,x2,y2)
              endif

!             plot number
              call setusv_dum(2hin,7)  ! yellow
              rarg = scale_l + (scale_h-scale_l) * frac

              if(l_set_contours)then
                  rlabel = vals(ifrac)
              elseif(log_scaling)then
                  rlabel = (10.**(rarg)) / scale
              else
                  rlabel = rarg / scale
              endif

              if(abs(rlabel) .lt. 0.999)then
                  rarg = rlabel*10.
                  if(rarg .ne. nint(rarg))then
                      write(ch_frac,3)rlabel
                  else
                      write(ch_frac,2)rlabel
                  endif
              elseif( (abs(rlabel) .lt. 2.0 
     1                .or. rlabel .ne. nint(rlabel) )       
     1                         .and.   colortable .ne. 'haines'    )then       
                  write(ch_frac,2)rlabel
              else
                  write(ch_frac,1)nint(rlabel)
              endif

              call left_justify(ch_frac)
              call s_len(ch_frac,len_frac)

              ixm = ixl + (ixh-ixl)*frac
              call zoomit(plot_parms,zoom_colorbar
     1                   ,cpux(ixm),cpux(iy),rxm,ry)
              call pchiqu (rxm,ry,ch_frac(1:len_frac)
     1                    ,rsize/zoom_colorbar,0 , 0.0)       
          enddo ! ifrac


      endif ! l_loop

      return
      end 


      subroutine generate_colortable(ncols,colortable,iwkid,icol_offset       
     1                              ,power,plot_parms,istatus)

!     generate colortable from ramp information stored in a file
!     number of colors is passed in

      include 'lapsplot.inc'

      character*200 path,filename
      character*(*)colortable
      real frac_a(300)
      real hue_a(300)
      real sat_a(300)
      real rint_a(300)

      write(6,*)' subroutine generate_colortable: ncols = ',ncols

      if(plot_parms%ncols .gt. 0)then
          ncols = plot_parms%ncols
          write(6,*)' resetting ncols to plot_parms value of ',ncols       
      endif

      lun = 41
      nramp = 0

      call get_directory('static',path,lenp)
      call s_len(colortable,lenc)

      lenc2 = lenc
      do i = 2,lenc
          if(colortable(i:i) .eq. '_')then ! truncate at the underscore
              lenc2 = i-1
          endif
      enddo

      filename = path(1:lenp)//'/www/'//colortable(1:lenc2)//'.lut'

      call s_len(filename,lenf)

      write(6,*)' reading colortable ',filename(1:lenf)
      open(lun,file=filename(1:lenf),status='old')

 10   read(lun,*,end=101)frac,hue,sat,rint

      nramp = nramp + 1

      frac_a(nramp) = frac**power
      write(6,*)' generate_colortable ',nramp,power,frac,frac_a(nramp)
      hue_a(nramp) = hue
      sat_a(nramp) = sat
      rint_a(nramp) = rint

      go to 10

 101  continue

      if(.true.)then ! for each color, interpolate from the ramp elements
          do i = 1,ncols
              frac_color = float(i-1) / float(ncols-1)
              do j = 1,nramp-1
                  if(frac_a(j)   .le. frac_color .and. 
     1               frac_a(j+1) .ge. frac_color      )then !interpolate
                      frac_ramp = (frac_color  - frac_a(j)) / 
     1                            (frac_a(j+1) - frac_a(j))
                      hue_col  = hue_a(j)    * (1. - frac_ramp) 
     1                         + hue_a(j+1)  *       frac_ramp
                      sat_col  = sat_a(j)    * (1. - frac_ramp) 
     1                         + sat_a(j+1)  *       frac_ramp
                      rint_col = rint_a(j)   * (1. - frac_ramp) 
     1                         + rint_a(j+1) *       frac_ramp
                      write(6,*)' i/frac_color = ',i,frac_color
                      call color_ramp(i,i,iwkid,icol_offset,plot_parms    
     1                               ,hue_col,sat_col,rint_col
     1                               ,hue_col,sat_col,rint_col)
                  endif
              enddo ! j
          enddo ! i
      else ! for each ramp element - assign a range of colors
          do i = 1,nramp-1
              i1 = 1 + nint(frac_a(i  ) * float(ncols-1))
              i2 = 1 + nint(frac_a(i+1) * float(ncols-1))
     
              call color_ramp(i1,i2,iwkid,icol_offset,plot_parms
     1                       ,hue_a(i  ),sat_a(i  ),rint_a(i  )
     1                       ,hue_a(i+1),sat_a(i+1),rint_a(i+1) )

          enddo ! i

      endif

      close(lun)

      return
      end


      subroutine get_colorbar_int(range,colorbar_int,l_divisible)

      logical l_divisible

      if(range .ge. 13000.)then
          colorbar_int = 2000.
      elseif(range .ge. 4000.)then
          colorbar_int = 1000.
      elseif(range .ge. 2500.)then
          colorbar_int = 500.
      elseif(range .ge. 1200.)then
          colorbar_int = 400.
      elseif(range .ge. 700.)then
          colorbar_int = 100.
      elseif(range .ge. 200.)then
          colorbar_int = 50.
      elseif(range .ge. 190.)then
          colorbar_int = 20.
      elseif(range .ge. 45.)then
          colorbar_int = 10.
      elseif(range .ge. 25.)then
          colorbar_int = 5.
      elseif(range .ge. 10.)then
          colorbar_int = 2.
      elseif(range .ge. 6.)then
          colorbar_int = 1.
      elseif(range .ge. 1.5)then
          colorbar_int = 0.5
      else !
          colorbar_int = 0.1
      endif

      rints = range / colorbar_int
      resid = abs(rints - nint(rints))

      if(resid .lt. .001)then
          l_divisible = .true.
      else
          l_divisible = .false.
      endif

      write(6,*)' returning from get_colorbar_int ',range,colorbar_int       

      return
      end

     
      subroutine get_pcp_vals(maxvals,namelist_parms,nvals,vals,itype)

      include 'lapsplot.inc'

      real vals(maxvals)

      write(6,*)' get_pcp_vals: itype = ',itype

      if(itype .eq. 1)then ! rain
        if(namelist_parms%c_units_type .eq. 'english')then ! (0-10 inches)
          vals(1) = 0.
          vals(2) = .01
          vals(3) = .02
          vals(4) = .05
          vals(5) = .10
          vals(6) = .20
          vals(7) = .50
          vals(8) = 1.0
          vals(9) = 1.5
          vals(10) = 2.0
          vals(11) = 3.0
          vals(12) = 4.0
          vals(13) = 5.0
          vals(14) = 7.0
          vals(15) = 10.0

        else ! metric (mm)
          vals(1) = 0.
          vals(2) = .2
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 5.0
          vals(7) = 10.
          vals(8) = 20.
          vals(9) = 30.
          vals(10) = 50.
          vals(11) = 75.
          vals(12) = 100.
          vals(13) = 150.
          vals(14) = 200.
          vals(15) = 400.

        endif

      elseif(itype .eq. 2)then ! snow 
        if(namelist_parms%c_units_type .eq. 'english')then ! (0-40 inches)
          vals(1) = 0.
          vals(2) = .1
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 3.0
          vals(7) = 4.0
          vals(8) = 6.0
          vals(9) = 8.0
          vals(10) = 10.0
          vals(11) = 12.0
          vals(12) = 15.0
          vals(13) = 20.0
          vals(14) = 30.0
          vals(15) = 40.0

        else ! metric (mm)
          vals(1) = 0.
          vals(2) = .2
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 5.0
          vals(7) = 10.
          vals(8) = 20.
          vals(9) = 30.
          vals(10) = 50.
          vals(11) = 75.
          vals(12) = 100.
          vals(13) = 150.
          vals(14) = 200.
          vals(15) = 400.

        endif

      elseif(itype .eq. 3)then ! heavy rain
        if(namelist_parms%c_units_type .eq. 'english')then ! (0-25 inches)
          vals(1) = 0.
          vals(2) = .1
          vals(3) = .2
          vals(4) = .5
          vals(5) = 1.0
          vals(6) = 2.0
          vals(7) = 3.0
          vals(8) = 4.0
          vals(9) = 6.0
          vals(10) = 8.0
          vals(11) = 10.0
          vals(12) = 12.0
          vals(13) = 15.0
          vals(14) = 20.0
          vals(15) = 25.0

        else ! metric (mm)
          vals(1) = 0.
          vals(2) = .2
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 5.0
          vals(7) = 10.
          vals(8) = 20.
          vals(9) = 30.
          vals(10) = 50.
          vals(11) = 75.
          vals(12) = 100.
          vals(13) = 150.
          vals(14) = 200.
          vals(15) = 400.

        endif

      elseif(itype .eq. 4)then ! heavy snow
        if(namelist_parms%c_units_type .eq. 'english')then ! (0-100 inches)
          vals(1) = 0.
          vals(2) = .1
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 4.0
          vals(7) = 6.0
          vals(8) = 10.
          vals(9) = 15.
          vals(10) = 20.
          vals(11) = 30.
          vals(12) = 40.
          vals(13) = 50.
          vals(14) = 70.
          vals(15) = 100.

        else ! metric (mm)
          vals(1) = 0.
          vals(2) = .2
          vals(3) = .5
          vals(4) = 1.0
          vals(5) = 2.0
          vals(6) = 5.0
          vals(7) = 10.
          vals(8) = 20.
          vals(9) = 30.
          vals(10) = 50.
          vals(11) = 75.
          vals(12) = 100.
          vals(13) = 150.
          vals(14) = 200.
          vals(15) = 400.

        endif
      endif

      nvals = 15

      return
      end

      subroutine zoomit(plot_parms,zoom,x,y,x_out,y_out)

      include 'lapsplot.inc'

      frame_factx = 1.0  ! / 0.75
      zxcen = (0.5 + ((plot_parms%xcen - 0.5) * frame_factx))
      zycen = (0.5 + ((plot_parms%ycen - 0.5) * frame_factx))
      x_out = zxcen + ((x-0.5) / zoom)
      y_out = zycen + ((y-0.5) / zoom)

      write(6,*)' zoomit before ',zoom,x,y 
      write(6,*)' zoomit after ',zoom,x_out,y_out 

      return
      end
