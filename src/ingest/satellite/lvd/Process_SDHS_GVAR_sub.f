      subroutine process_sdhs_gvar_sub(filename,chtype,isat,jtype,
     &nelem,nlines,image_out,nlfi,nefi,depth,width,istatus)

c**********************************************************************
c  purpose:  convert the sdhs image data from a cellular format to an 
c  image format
c
c  method:  read the image header information to determine how many
c  64 scanline by 256 pixel cells make up the image.  read in the image
c  into a singlely dimnesioned array and convert the data into a two
c  dimensional array where the first index indicates the row (scanline
c  number) and the second dimension the column(pixel number).  some 
c  fortran 90 constructs are used to allow dynamic array sizing
c  but the majority of the code is in fortran 77. 
c  no attempt is made to remove the blank parts of cells
c  which may be in the cells send to satisfy an image request.  also
c  since aix fortran will not allow an unsigned byte variable the data
c  is read into an integer*1 variable and the bit string picked off and
c  stored into an integer*2 variable.  without using this method byte 
c  values such as 255 (all bits on) were being interpreted as -128.
c  
c  references:
c  1.  final interface specification for the satellite data handling
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro
c
c**********************************************************************

      implicit none

      character*255 filename
      character*3   chtype

      integer cell_width !number of pixels in width of a cell
      integer cell_depth !number of scanline in depth of cell
      integer header_size  !size of gvar header
      integer pixels_per_cell !number of pixels per cell
      integer nelem,nlines
      integer nefi,nlfi
      integer imagelen
      integer nf
      integer lend
      integer width      !number of cells in width of image
      integer depth      !number of cells in depth of image
      integer itempintgr

      parameter (header_size = 256) 
      parameter (cell_width = 256)
      parameter (cell_depth = 64)
      parameter (pixels_per_cell = cell_width*cell_depth)

c     parameter (width=5,depth=13)  !this for the test images so far.
c     parameter (imagelen=width*cell_width*depth*cell_depth)

c     parameter (imagelen_ir=nlines_ir_max*nelem_ir_max)
c     parameter (imagelen_vis=nlines_vis_max*nelem_vis_max)
c     parameter (imagelen_wv=nlines_wv_max*nelem_wv_max)

      integer size       !size of the sdhs image including the header
c     integer cell_num   !the number of the cell
c     integer cell_row   !the row number of a cell
c     integer cell_column !the column number of a cell
c     integer xs         !pixels in width of image
c     integer ys         !scanlines depth of image

      integer header(header_size)!array used to read over the header
      integer iostatus   !i/o status variable
      integer istatus    !return status
      integer i,j        !do loop indices
      integer rnum       !the row number (scanline) of the converted
                           !image
      integer cnum       !the column number(pixel) of the converted 
                           !image
      integer ii,jj,m
      integer ispec
      integer istat
      integer istart,jstart
      integer iend,jend
      integer isat,jtype

c**********************************************************************
c  fortran 90 used to declare allocatable arrays 
c**********************************************************************

c     integer*1,dimension(:),allocatable:: imagei1 
c     integer*2,dimension(:,:),allocatable:: image 
c static case: jsmart 5-12-97
c     integer*1 imagei1(width*cell_width*depth*cell_depth)

      character      imagec(nlfi*nefi)*1
      integer        image_decell_2d(nefi,nlfi)
      integer        image_temp(nefi,nlfi)
      integer        image_decell_1d(nefi*nlfi+header_size*4)
      real           image_out(nelem,nlines)
      real           r8to10

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus=0

      call lvd_file_specifier(chtype,ispec,istat)
      if(ispec.eq.1)then
       istart=i_start_vis(jtype,isat)
       jstart=j_start_vis(jtype,isat)
       iend=i_end_vis(jtype,isat)
       jend=j_end_vis(jtype,isat)
      elseif(ispec.eq.3)then
       istart=i_start_wv(jtype,isat)
       jstart=j_start_wv(jtype,isat)
       iend=i_end_wv(jtype,isat)
       jend=j_end_wv(jtype,isat)
      elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then
       istart=i_start_ir(jtype,isat)
       jstart=j_start_ir(jtype,isat)
       iend=i_end_ir(jtype,isat)
       jend=j_end_ir(jtype,isat)
      endif
c
c**********************************************************************
c  read the header information and pass back the number of 64 scanline
c  by 256 pixel cells in the width and depth of the image.
c**********************************************************************

c     call readhdr(filename,width, depth)

c**********************************************************************
c  compute the number of pixels in the width of the image and the 
c  number of scanline in the depth of the image.  also compute the size
c  of the image file including the header. 
c**********************************************************************
c original constructs:
c     xs =width*cell_width
c     ys = depth*cell_depth
c     size = ys*xs+header_size

      size = nlfi*nefi+header_size*4   !note header_size = 256*4=1024
c
c**********************************************************************
c  open the file that contains the gvar header and pixel data.  the
c  record length in the direct access read is set equal to the size
c**********************************************************************
c
      nf=index(filename,' ')-1
      open(unit=8, file=filename, err=100, iostat=iostatus,
     & access='direct', form='unformatted',recl= size)

      if(l_cell_afwa)then
         imagelen=(nlfi*nefi)/4
         read (8,rec=1,err=99)header,imagec
      else
         imagelen=nlfi*nefi
         call read_binary_field(image_decell_1d,1,4,imagelen+
     +header_size*4,filename,nf)
      endif

      close(8)

      if(l_cell_afwa)then

         call decellularize_image(imagec,imagelen,width,depth,
     &pixels_per_cell,cell_width,cell_depth,nefi,nlfi,image_decell_2d)

c back as floating 10-bit info for the sector in domain.

         jj=0
         do j=jstart,jend
            jj=jj+1
            ii=0
         do i=istart,iend
            ii=ii+1
            image_out(ii,jj)=float(image_decell_2d(i,j))*4.0
         enddo
         enddo

      else ! afwa data not cellularized but bits need moving

c goes data = 10 bit; meteosat data = 8 bit.
         r8to10=4.0
         if(isat.eq.2 .and. jtype.eq.4)r8to10=1.0

         m=header_size*4
         do j=1,nlfi
         do i=1,nefi
            m=m+1
            image_decell_2d(i,j)=image_decell_1d(m)
         enddo
         enddo

         if(isat.eq.2 .and. jtype.eq.4)then

c meteosat origin is se corner (wrt ri/rj lut). make it sw corner
            do j=1,nlfi
            do i=1,nefi
               image_temp(nefi-i+1,j)=image_decell_2d(i,j)
            enddo
            enddo
            do j=1,nlfi
            do i=1,nefi
               image_decell_2d(i,j)=image_temp(i,j)
            enddo
            enddo

         endif
c
c now load that part within the laps domain; convert to 10-bit.
c
         jj=0
         do j=jstart,jend
            jj=jj+1
            ii=0
         do i=istart,iend
            ii=ii+1
            itempintgr=ibits(image_decell_2d(i,j),0,8)
            image_out(ii,jj)=float(itempintgr)*r8to10
         enddo
         enddo

      endif

      istatus = 1
      goto 1000

99    write(6,*)'error reading gwc file'
      write(6,*)'filename = ',filename(1:nf)
      close(8)
      return 

100   if (iostatus .ne. 0)then
        print *, 'error reading ',filename, ' io status is', 
     &  iostatus
        istatus = -1
      end if

1000  return
      end
