      subroutine decellularize_image(imagein,imagelen,
     &width,depth,pixels_per_cell,cell_width,cell_depth,
     &nelem,nlines,imageout)

      implicit none

      integer nelem,nlines
      integer width,depth
      integer imagelen
      integer cell_width,cell_depth
      integer pixels_per_cell
      integer i,j
      integer cell_num,cell_row,cell_column
      integer rnum,cnum

      integer      imagein(imagelen)

      integer        imageout(nelem,nlines)

      do j=1,nlines
      do i=1,nelem
         imageout(i,j)=0
      enddo
      enddo

      do i=1,nlines*nelem

        cell_num = (i-1)/pixels_per_cell + 1
        cell_row = (i-1)/(pixels_per_cell*width) +1
        cell_column = mod( (i-1)/pixels_per_cell, width ) + 1

        rnum = (i-((cell_num-1)*pixels_per_cell)-1)/cell_width
     &         + (cell_row-1)*cell_depth + 1

        cnum = mod( (i-((cell_num-1)*pixels_per_cell)-1), cell_width)
     &         + 1 + (cell_column-1)*256

cc        imageout(cnum,rnum)=ibits(imagein(i),0,8)
        imageout(cnum,rnum)=ibits(imagein(1+(i-1)/4),24-8*mod(i-1,4),8)

      end do

      return
      end
