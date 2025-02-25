      subroutine get_ascii_dimensions(cdatapath,
     &c_channel_type,nelems,nlines,istatus)
c
      implicit none

      integer     nlines,nelems

      integer     i,j,k
      integer     iprev
      integer     n,nc

      real        itb
      real        img_lin
      real        img_lin_prev
      real        img_ele

      integer     istatus

      logical       eof

      real        xlat,xlon

      character     c_date*5
      character     c_time*4
      character     c_filetime*9
      character     c_filename*255
      character     c_channel_type*3
      character     cdatapath*(*)

      istatus=1
c
c build filename and open ascii file
c
      n = index(cdatapath,' ')-1
      nc= index(c_channel_type,' ')-1
      if(nc.le.0)nc=3
      c_filename = cdatapath(1:n)//'laps_'//c_channel_type(1:nc)
      open(19,file=c_filename,form='formatted',status='old',
     &err=900)
c
c first line is date/time
c
      read(19,102,err=901)c_date,c_time
102   format(1x,a5,1x,a4)
      c_filetime=c_date//c_time
      write(6,*)c_filetime
c
c read second line to set up the line and element counting
c
      read(19,*,err=901)xlat,xlon,img_lin_prev,img_ele,itb
      eof=.false.
      iprev=1
      j=1

      do while (.not.eof)

         read(19,*,end=29)xlat,xlon,img_lin,img_ele,itb
         if(img_lin.ne.img_lin_prev)then
            j=j+1
            if(iprev.gt.i)i=iprev
            iprev=0
            img_lin_prev=img_lin
         endif
 
         iprev=iprev+1
         goto 200

29       eof=.true.

200   enddo

      close(19)
c
      goto 995

900   write(6,*)'error opening data file',c_filename
      istatus=-1
      goto 1000

901   write(6,*)'error - initial read data file ',c_filename
      istatus=-1
      goto 1000

995   write(6,*)'finished. nelem/nlines: ',i,j
      write(6,*)

      nelems=i
      nlines=j

1000  close(19)
      return
      end
