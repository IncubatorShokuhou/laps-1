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


      subroutine conrec_line(field,ni,nii,nj
     1                      ,clow_in,chigh_in,cint_in,plot_parms
     1                      ,i1,i2,i3,i4)

!     97-aug-17     ken dritz     commented out assignment to r_missing_data
!                                 and inserted call to get_r_missing_data

      include 'lapsplot.inc'

      data ihl/0/
      save ihl

      real field(ni,nj)

      character*48 c_blank

      c_blank = '   '

      nf = 51

      call gflas1(nf)
      call gflas2()

!     call color

!     r_missing_data = 1.e37
      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'error getting r_missing_data in conrec_line'
         stop
      endif

!     ihl = 2
      ihl = ihl + 1

      if(chigh_in .eq. clow_in)then ! special contouring (non-equal intervals)
          clow = clow_in
          chigh = chigh_in
          cint = cint_in
          lis = 1 ! labelling interval
!         ihl = ihl - 1
      else
!         ihl = 1
          clow = clow_in
          chigh = chigh_in
          cint = cint_in
          lis = 2
      endif

      write(6,*)' conrec_line: calling plot_contour'
     1                ,clow,chigh,cint,lis,ihl

      call plot_contour
     1  (nf,ni,nii,nj,field,c_blank,cint,cl,ch,c_blank,4,ihl
     1      ,clow,chigh,plot_parms,lis,r_missing_data)

      return
      end



      subroutine plot_contour(nf, mx, nx, ny, f, lab1, fint, fmin, fmax,
     + timelabel,nc,ihl,clow,chigh,plot_parms,lis,r_missing_data)
c
c purpose       plot contours of f(nx,ny).
c
c arguments     nf --- i  flash buffer to be used
c               mx --- i  maximum size of first dimension of f
c               nx --- i  actual size of first dimension of f
c               ny --- i  actual size of second dimension of f
c               f ---- i  two dimensional array to be contoured
c               lab1 - i  user-supplied plot label
c               fint - i  contour interval (if 0., .1*(fmax-fmin) is used)
c               fmin -  o minimum grid value
c               fmax -  o maximum grid value
c               timelabel - i string
c               nc   - i  contour color index
c               ihl  - i  flag as whether to plot high/low labels
c
c


      include 'lapsplot.inc'

      parameter (n2x=400)
      integer nf, mx, nx, ny, ivar
      real f1(n2x,n2x),hold(n2x,2)
      real f(mx,ny), fint, fmin, fmax, machine_epsilon_p

      character*48 lab3
      character*48 timelabel
      character*48 lab1

      parameter (lrwk=64000, liwk=64000)

      real rwrk(lrwk)
      integer iwrk(liwk)

      integer hic,hlc

      real scale
      data scale /1.1/

      parameter (machine_epsilon_p = 1.19e-07)  ! from iftran.im - bj
      common /savmap/
     .  mx1, mx2, my1, my2,
     .  u1,  u2,  v1,  v2
      common /fxfy1/ xa, ya, ua, va, u1a, v1a, dudi, dvdj

        common /nxny/ nx1,ny1
        common /icol_index/ icol_current
        common /zoom/       zoom

c ... compare,spv def'n from iftran.im/b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      spv(a) = compare(a,spval_p,0.01)

        spval_p=r_missing_data

        fmin0 = fmin
        fmax0 = fmax

        nx1 = nx
        ny1 = ny
        nx2x = nx*2-1
        ny2x = ny*2-1

!       maxchr = 'x'
!       minchr = 'n'

c
c     calc max and min of f
c
      fmax = -1.e8
      fmin =  1.e8
      do 55 j=1,ny
        do 60 i=1,nx
!         if (spv(f(i,j)) .eq. no_p)
          if (spv(f(i,j)) .eq. 0)
     .    then
!          if(f(i,j) .ne. r_missing_data)then
            if (f(i,j).gt.fmax) fmax=f(i,j)
            if (f(i,j).lt.fmin) fmin=f(i,j)
!          endif
          end if
60      continue
55    continue
c
c     test for zero fint
c
      if (fint .le. machine_epsilon_p) then
        df = .1*(fmax-fmin)
      else
        df = fint
      end if

      write(6,5) fmin,fmax,df
5     format(/,' min = ',g12.5,5x,'max = ',
     .  g12.5,5x,' plot interval = ',g12.5)
c
c    plot if there is adequate field variation
c

        write(lab3,101)  fmax,fmin,df
  101   format(2x,'max=',f10.3,2x,'min=',f10.3,2x,'int=',f8.3,2x)

        call gflas3(nf)

!       call set (0.,1.,0.,1.,0.,1.,0.,1.,1)

        if (nc.eq.1) call plchhq (0.0,0.88,timelabel,0.018,0.,-1.)
!       call pcseti ('cc - character color', 7+nc)
!       call plchhq (0.45,0.14-(nc-1)*0.02, lab3   , .01,
!    1       0., -1.)
        call plchhq (0.05,0.14-(nc-1)*0.02, lab1   , .01,
     1       0., -1.)

!       fmin = ceiling2(fmin,df,machine_epsilon_p)   ! minimum contour
        spvalu=spval_p
!       call fxfy(nf)                               ! initialize fx and fy fcns



        nulbll = 1         ! number of unlabelled lines between labelled ones
        nhi = -1

c --- do contouring
!       if (ihl.ge.1) then

        if (.true.) then
!         call cpsetr('hls  - high/low label size',.020)
          call cpsetr('hls  - high/low label size',.013)
          call cpsetc('ilt',' ')

          if(fmax .eq. 100. .and. fmin .le. 5.0)then ! rh case
              nsd = 3
              nls = 1
          elseif(fmax .ge. 1000. .and. fmax .le. 1100.
     1     .and. fmin .gt. 500.  .and. fmin .lt. 1000.)then ! stnp case
              nsd = 3
              nls = 1
          else ! general case
              nsd = 2
              nls = 1
          endif

          call cpseti('nsd',nsd)
          call cpseti('nls',nls)

          write(6,*)' fmax/nsd/nls = ',fmax,nsd,nls

          call cpseti ('hic', icol_current)
          call cpseti ('loc', icol_current)
          call cpseti ('hlc', icol_current)

!         if (ihl.ge.2) call cpsetc('hlt',' ')

        end if

!       call cpsetr ('lls',.040)
        call cpseti ('cls - contour level selection flag',20)
        call cpsetr ('cis',df)
        call cpseti ('lis - label interval specifier',lis)
!       call cpseti ('llp', 2)
!       call cpsetr ('cmn',fmin)
!       call cpsetr ('cmx',fmax)
        call cpsetr ('cmn',clow)
        call cpsetr ('cmx',chigh)
        call cpsetr ('spv',spvalu)
!       call cpsetr ('lls - line label size',.025/sqrt(zoom))
        call cpsetr ('lls - line label size',.025)
        call cpgetr ('lls - line label size',clls)           
        call cpgetr ('hls',hls)           
        call cpsetr ('cwm',1.00/sqrt(zoom))


        call cprect (f,nx,nx,ny,rwrk,lrwk,iwrk,liwk)
        call cppkcl (f,rwrk,iwrk)
       if (ihl.ge.1) then
        call cplbdr (f,rwrk,lwrk)
       end if


        call cpgeti ('ncl - number of contour levels', ncon)
        do 111 i=1,ncon
          call cpseti ('pai - parameter array index', i)
          call cpgetr ('clv - contour level values',cval)
          if (cval.lt.0.) then
            call cpseti ('cld',61166)
          else
            call cpseti ('cld',65535)
          end if
          call cpseti ('cll',plot_parms%contour_line_width)       ! line width
!         call cpseti ('clc - contour line color index', 7+nc)
111     continue

        call cpcldr (f,rwrk,iwrk)

        isize=7
        call pwrit (mx1,my2,' ',1,0,0,0)
        call sflush                                  ! flush frame buffer

        call cpgeti ('liu', liu)
        call cpgeti ('llp', llp)
        call cpgeti ('hic', hic)
        call cpgeti ('hlc', hlc)
        call cpgeti ('loc', loc)

        write(6,*)'ncon/ihl / hls / hic/hlc/loc/llp/lls/lis/liu/icol = '
        write(6,*) ncon,ihl,hls,hic,hlc,loc,llp,clls,lis,liu
     1             ,icol_current


        write(6,155)
155     format(' the field has been plotted')


      return
      end

cccccccccccccccccccccccccccccc  ceiling2  cccccccccccccccccccccccccccc
c
c  ceiling2 - ceiling of x, over interval y, to accuracy z
c
        real function ceiling2(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0) then                ! ceiling(x,y,z) = trunc(x+y-z,y)
          ceiling2=trunc(x+y-z,y)
        else if (x.lt.0) then           ! ceiling(x,y,z) = -trunc(-x,y)
          ceiling2= -trunc(-x,y)
        else                            ! ceiling(x,y,z) = 0.
          ceiling2 = 0.
        endif
c
        return
        end


        subroutine fill (xwrk,ywrk,nwrk,iarea,igrp,ngrps)
c
        dimension xwrk(*),ywrk(*),iarea(*),igrp(*)

        do 10, i=1,ngrps
          if (igrp(i).eq.3) iarea3=iarea(i)
 10     continue

        if (iarea3 .gt. 0) then
c if the area is defined by 3 or more points, fill it
           call gsfaci(iarea3+1)
           call gfa(nwrk,xwrk,ywrk)
        endif
        return
        end

      subroutine color(iwhite)
!     write(6,*)' white or black background. [enter 1/2]    ?   :'
!     accept 5, iwhite
!5     format (i)
!     iwhite = 2

c    background color
c  the background is white here for better visibility on paper

      if (iwhite.eq.1) then
        write(6,*)' color: white background'
        call gscr (1,0,1.,1.,1.)
        call gscr (1,1,0.8,0.8,1.)
      else
        write(6,*)' color: black background'
        call gscr (1,0,0.,0.,0.)
        call gscr (1,1,0.,0.,0.)
      end if

c
c     background color
c     black
c     call gscr(1,0,0.,0.,0.)
c
c     forground colors
c white
      call gscr(1,  1, 1.0, 1.0, 1.0)
        call gscr (1,1,0.8,0.8,1.)
c white
      call gscr(1,  2, 1.0, 1.0, 1.0)
c red
      call gscr(1,  3, 1.0, 0.0, 0.0)
c orangered
      call gscr(1,  4, 1.0, 0.30, 0.0)
c orange
      call gscr(1,  5, 1.0, 0.65, 0.0)
c gold
      call gscr(1,  6, 1.0, 0.85, 0.0)
c yellow
      call gscr(1,  7, 1.0, 1.0, 0.0)
c greenyellow
      call gscr(1,  8, 0.7, 1.0, 0.2)
c chartreuse
      call gscr(1,  9, 0.5, 1.0, 0.0)
c celeste
      call gscr(1, 10, 0.2, 1.0, 0.5)
c green
!     call gscr(1, 11, 0.2, 0.8, 0.2)
      call gscr(1, 11, 0.1, 0.9, 0.1)
c aqua
      call gscr(1, 12, 0.0, 0.9, 1.0)
c deepskyblue
      call gscr(1, 13, 0.0, 0.75, 1.0)
!     call gscr(1, 13, 0.0, 0.65, 0.9)
c royalblue
      call gscr(1, 14, 0.25, 0.45, 0.95)
c slateblue
      call gscr(1, 15, 0.4, 0.35, 0.8)
c darkviolet
      call gscr(1, 16, 0.6, 0.0, 0.8)
c lavender
      call gscr(1, 17, 1.00, 0.65, 1.0)
c black
      call gscr (1,21,0.,0.,0.)
c quasi black
      call gscr (1,22,0.05,0.05,0.05)
c medium gray
      call gscr (1,23,.3,.3,.3)
! unused...
      call gscr (1,24,.0,.7,.0)
! dark gray
      call gscr (1,25,.1,.1,.1)
! pale orange
      call gscr (1,26,.6,.3,.1)
c orchid
      call gscr (1,27, 0.85, 0.45, 0.8)
! unused...
      call gscr (1,28,.0,.0,.95)
! paler yellow
      call gscr (1,29,.50,.38,.0)
! dark pink
      call gscr (1,30,.5,.30,.5)
      call gscr (1,31,.2,.7,.7)

! pale yellow
      call gscr (1,32,.6,.6,.0)
! pale red
      call gscr (1,33,.5,.0,.0)
! grey
      call gscr (1,34,.5,.5,.5)
! darkviolet
      call gscr(1, 35, 0.6, 0.0, 0.8)

! radar vel color table (about 20 values)
      i_velcol_offset = 100
      do icol = i_velcol_offset,i_velcol_offset+10
          ivel = i_velcol_offset+10 - icol
          amp = float(ivel) / 10.
          amp_gry = 0.3
          amp_grn = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,amp_gry,amp_grn,amp_gry)
!         call gscr(1,icol,0.,1.0,0.)
!         write(6,*)' setting color ' ,icol
      enddo
      do icol = i_velcol_offset+11,i_velcol_offset+20
          ivel = icol - (i_velcol_offset+10)
          amp = float(ivel) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,amp_red,amp_gry,amp_gry)
!         write(6,*)' setting color ' ,icol
      enddo

! radar ref color table
      i_refcol_offset = 180
      do idbz = 0,20
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz) / 20.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_blu = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,amp_red,amp_gry,amp_blu)
      enddo
      do idbz = 21,30
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-20) / 10.
          amp_gry = 0.3
          amp_grn = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,0.,amp_grn,0.)
      enddo
      do idbz = 31,40
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-30) / 10.
          amp_gry = 0.3
          amp_blu = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,0.,0.,amp_blu)
      enddo
      do idbz = 41,50
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-40) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_grn = amp_red * 0.5
          call gscr(1,icol,amp_red,amp_grn,0.)
      enddo
      do idbz = 51,60
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-50) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_grn = amp_red
          call gscr(1,icol,amp_red,amp_grn,0.)
      enddo
      do idbz = 61,70
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-60) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call gscr(1,icol,amp_red,0.,0.)
      enddo

!     white
      call gscr(1, 251, 1.0, 1.0, 1.0)

c done.
c
        return
c
      end
