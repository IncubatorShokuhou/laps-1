      module g2grids
!$$$  subprogram documentation block
!                .      .    .                                       .
! module:    g2grids 
!   prgmmr: gilbert         org: w/np11    date: 2004-04-27
!
! abstract: this fortran module allows access to predefined grib2 grid 
!   definition templates stored in a file.  the gdts are represented by 
!   a predefined number or a character abbreviation.
!
!   at the first request, all the grid gdt entries in the file associated
!   with input fortran file unit number, lunit, are read into a linked list
!   named gridlist.  this list is searched for the requested entry.
!
!   users of this fortran module should only call routines getgridbynum
!   and getgridbyname.
!
!   the format of the file scanned by routines in this module is as follows.
!   each line contains one grid entry containing five fields, each separated
!   by a colon, ":".  the fields are:
!      1) - predefined grid number
!      2) - up to an 8 character abbreviation
!      3) - grid definition template number
!      4) - number of entries in the grid definition template
!      5) - a list of values for each entry in the grid definition template.
!
!   as an example, this is the entry for the 1x1 gfs global grid 
!   3:gbl_1deg:  0:19: 0 0 0 0 0 0 0 360 181 0 0 90000000 0 48 -90000000 359000000 1000000 1000000 0
!
!   comments can be included in the file by specifying the symbol "#" as the
!   first character on the line.  these lines are ignored.
!
!
! program history log:
! 2004-04-27  gilbert
!
! usage:    use g2grids
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      integer,parameter :: maxtemp=200

      type,private :: g2grid
          integer :: grid_num
          integer :: gdt_num
          integer :: gdt_len
          integer,dimension(maxtemp) :: gridtmpl
          character(len=8) :: cdesc
          type(g2grid),pointer :: next
      end type g2grid

      type(g2grid),pointer,private :: gridlist
      integer :: num_grids=0

      contains


         integer function readgrids(lunit)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    readgrids
!   prgmmr: gilbert         org: w/np11    date: 2001-06-28
!
! abstract: this function reads the list of gdt entries in the file 
!   associated with fortran unit, lunit.  all the entries are stored in a
!   linked list called gridlist.
!
! program history log:
! 2001-06-28  gilbert
!
! usage:    number=readgrids(lunit)
!   input argument list:
!     lunit   - fortran unit number associated the the gdt file.
!
! returns:  the number of grid definition templates read in.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: lunit

           integer,parameter :: linelen=1280
           character(len=8) :: desc
           character(len=linelen) :: cline
           integer  ient,igdtn,igdtmpl(200),igdtlen
           integer :: pos1,pos2,pos3,pos4

           type(g2grid),pointer :: gtemp
           type(g2grid),pointer :: prev
           integer count

           count=0

           !   for each line in the file....
           do
             !  read line into buffer
             !
             cline(1:linelen)=' '
             read(lunit,end=999,fmt='(a)') cline

             !
             !  skip line if commented out
             !
             if (cline(1:1).eq.'#') cycle

             !
             !  find positions of delimiters, ":"
             !
             pos1=index(cline,':')
             cline(pos1:pos1)=';'
             pos2=index(cline,':')
             cline(pos2:pos2)=';'
             pos3=index(cline,':')
             cline(pos3:pos3)=';'
             pos4=index(cline,':')
             if ( pos1.eq.0 .or. pos2.eq.0 .or. pos3.eq.0 .or. 
     &            pos4.eq.0) cycle

             !
             !  read each of the five fields.
             !
             read(cline(1:pos1-1),*) ient
             read(cline(pos1+1:pos2-1),*) desc
             read(cline(pos2+1:pos3-1),*) igdtn
             read(cline(pos3+1:pos4-1),*) igdtlen
             read(cline(pos4+1:linelen),*) (igdtmpl(j),j=1,igdtlen)

             !
             !  allocate new type(g2grid) variable to store the gdt
             !
             allocate(gtemp,stat=iom)
             count=count+1
             gtemp%grid_num=ient
             gtemp%gdt_num=igdtn
             gtemp%gdt_len=igdtlen
             gtemp%gridtmpl=igdtmpl
             gtemp%cdesc=desc
             nullify(gtemp%next)              ! defines end of linked list.
             if ( count .eq. 1 ) then
                gridlist => gtemp
             else                       ! make sure previous entry in list
                prev%next => gtemp      ! points to the new entry,
             endif
             prev => gtemp

           enddo

 999       readgrids=count
           return

         end function


         subroutine getgridbynum(lunit,number,igdtn,igdtmpl,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getgridbynum 
!   prgmmr: gilbert         org: w/np11    date: 2004-04-26
!
! abstract: this subroutine searches a file referenced by fortran unit lunit
!   for a grid definition template assigned to the requested number. 
!   the input file format is described at the top of this module.
!
! program history log:
! 2004-04-26  gilbert
!
! usage:    call getgridbynum(lunit,number,igdtn,igdtmpl,iret)
!   input argument list:
!     lunit    - unit number of file containing grid definitions 
!     number   - grid number of the requested grid definition 
!
!   output argument list:      
!     igdtn    - nn, indicating the number of the grid definition 
!                template 3.nn
!     igdtmpl()- an array containing the values of each entry in
!                the grid definition template.
!     iret     - error return code.
!                0 = no error
!               -1 = undefined grid number.
!                3 = could not read any grids from file.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: lunit,number
           integer,intent(out) :: igdtn,igdtmpl(*),iret

           type(g2grid),pointer :: tempgrid

           iret=0
           igdtn=-1
           !igdtmpl=0

           !
           !  if no grids in list, try reading them from the file.
           !
           if ( num_grids .eq. 0 ) then
              num_grids=readgrids(lunit)
           endif

           if ( num_grids .eq. 0 ) then
              iret=3                         ! problem reading file
              return
           endif

           tempgrid => gridlist

           !
           !  search through list
           !
           do while ( associated(tempgrid) )
               if ( number .eq. tempgrid%grid_num ) then
                  igdtn=tempgrid%gdt_num
                  igdtmpl(1:tempgrid%gdt_len)=
     &                        tempgrid%gridtmpl(1:tempgrid%gdt_len)
                  return
               else
                  tempgrid => tempgrid%next
               endif
           enddo 
 
           iret=-1
           return
 
         end subroutine


         subroutine getgridbyname(lunit,name,igdtn,igdtmpl,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getgridbyname 
!   prgmmr: gilbert         org: w/np11    date: 2004-04-26
!
! abstract: this subroutine searches a file referenced by fortran unit lunit
!   for a grid definition template assigned to the requested name. 
!   the input file format is described at the top of this module.
!
! program history log:
! 2004-04-26  gilbert
!
! usage:    call getgridbyname(lunit,name,igdtn,igdtmpl,iret)
!   input argument list:
!     lunit    - unit number of file containing grid definitions 
!     name     - grid name of the requested grid definition 
!
!   output argument list:      
!     igdtn    - nn, indicating the number of the grid definition 
!                template 3.nn
!     igdtmpl()- an array containing the values of each entry in
!                the grid definition template.
!     iret     - error return code.
!                0 = no error
!               -1 = undefined grid number.
!                3 = could not read any grids from file.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: lunit
           character(len=8),intent(in) :: name
           integer,intent(out) :: igdtn,igdtmpl(*),iret

           type(g2grid),pointer :: tempgrid

           iret=0
           igdtn=-1
           !igdtmpl=0

           !
           !  if no grids in list, try reading them from the file.
           !
           if ( num_grids .eq. 0 ) then
              num_grids=readgrids(lunit)
           endif

           if ( num_grids .eq. 0 ) then
              iret=3                         ! problem reading file
              return
           endif

           tempgrid => gridlist

           !
           !  search through list
           !
           do while ( associated(tempgrid) )
               if ( name .eq. tempgrid%cdesc ) then
                  igdtn=tempgrid%gdt_num
                  igdtmpl(1:tempgrid%gdt_len)=
     &                     tempgrid%gridtmpl(1:tempgrid%gdt_len)
                  return
               else
                  tempgrid => tempgrid%next
               endif
           enddo 
 
           iret=-1
           return
 
         end subroutine


      end

