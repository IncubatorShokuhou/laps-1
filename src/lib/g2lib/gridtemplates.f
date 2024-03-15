      module gridtemplates
!$$$  subprogram documentation block
!                .      .    .                                       .
! module:    gridtemplates 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-09
!
! abstract: this fortran module contains info on all the available 
!   grib2 grid definition templates used in section 3 (gds).
!   each template has three parts: the number of entries in the template
!   (mapgridlen);  a map of the template (mapgrid), which contains the
!   number of octets in which to pack each of the template values; and
!   a logical value (needext) that indicates whether the template needs 
!   to be extended.  in some cases the number of entries in a template 
!   can vary depending upon values specified in the "static" part of 
!   the template.  ( see template 3.120 as an example )
!
!   this module also contains two subroutines.  subroutine getgridtemplate
!   returns the octet map for a specified template number, and
!   subroutine extgridtemplate will calculate the extended octet map
!   of an appropriate template given values for the "static" part of the 
!   template.  see docblocks below for the arguments and usage of these 
!   routines.
!
!   note:  array mapgrid contains the number of octets in which the 
!   corresponding template values will be stored.  a negative value in
!   mapgrid is used to indicate that the corresponding template entry can
!   contain negative values.  this information is used later when packing
!   (or unpacking) the template data values.  negative data values in grib
!   are stored with the left most bit set to one, and a negative number
!   of octets value in mapgrid() indicates that this possibility should
!   be considered.  the number of octets used to store the data value
!   in this case would be the absolute value of the negative value in 
!   mapgrid().
!  
!
! program history log:
! 2000-05-09  gilbert
! 2003-09-02  gilbert   -  added gdt 3.31 - albers equal area
!
! usage:    use gridtemplates
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      integer,parameter :: maxlen=200,maxtemp=23

      type gridtemplate
          integer :: template_num
          integer :: mapgridlen
          integer,dimension(maxlen) :: mapgrid
          logical :: needext
      end type gridtemplate

      type(gridtemplate),dimension(maxtemp) :: templates

      data templates(1)%template_num /0/     !  lat/lon 
      data templates(1)%mapgridlen /19/
      data templates(1)%needext /.false./
      data (templates(1)%mapgrid(j),j=1,19) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1/

      data templates(2)%template_num /1/     !  rotated lat/lon 
      data templates(2)%mapgridlen /22/
      data templates(2)%needext /.false./
      data (templates(2)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,4/

      data templates(3)%template_num /2/     !  stretched lat/lon 
      data templates(3)%mapgridlen /22/
      data templates(3)%needext /.false./
      data (templates(3)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,-4/

      data templates(4)%template_num /3/     !  stretched & rotated lat/lon 
      data templates(4)%mapgridlen /25/
      data templates(4)%needext /.false./
      data (templates(4)%mapgrid(j),j=1,25) 
     &       /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,4,-4,4,-4/

      data templates(5)%template_num /10/     !  mercator
      data templates(5)%mapgridlen /19/
      data templates(5)%needext /.false./
      data (templates(5)%mapgrid(j),j=1,19)
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,-4,-4,4,1,4,4,4/

      data templates(6)%template_num /20/     !  polar stereographic
      data templates(6)%mapgridlen /18/
      data templates(6)%needext /.false./
      data (templates(6)%mapgrid(j),j=1,18) 
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,-4,4,4,4,1,1/

      data templates(7)%template_num /30/     !  lambert conformal
      data templates(7)%mapgridlen /22/
      data templates(7)%needext /.false./
      data (templates(7)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,-4,4,4,4,1,1,-4,-4,-4,4/

      data templates(8)%template_num /40/     !  gaussian lat/lon
      data templates(8)%mapgridlen /19/
      data templates(8)%needext /.false./
      data (templates(8)%mapgrid(j),j=1,19) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1/

      data templates(9)%template_num /41/     !  rotated gaussian lat/lon
      data templates(9)%mapgridlen /22/
      data templates(9)%needext /.false./
      data (templates(9)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,4/

      data templates(10)%template_num /42/     !  stretched gaussian lat/lon
      data templates(10)%mapgridlen /22/
      data templates(10)%needext /.false./
      data (templates(10)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,-4/

      data templates(11)%template_num /43/     !  strtchd and rot'd gaus lat/lon
      data templates(11)%mapgridlen /25/
      data templates(11)%needext /.false./
      data (templates(11)%mapgrid(j),j=1,25) 
     &          /1,1,4,1,4,1,4,4,4,4,4,-4,4,1,-4,4,4,4,1,-4,4,4,-4,4,-4/

      data templates(12)%template_num /50/    !  spherical harmonic coefficients
      data templates(12)%mapgridlen /5/
      data templates(12)%needext /.false./
      data (templates(12)%mapgrid(j),j=1,5) /4,4,4,1,1/

      data templates(13)%template_num /51/   !  rotated spherical harmonic coeff
      data templates(13)%mapgridlen /8/
      data templates(13)%needext /.false./
      data (templates(13)%mapgrid(j),j=1,8) /4,4,4,1,1,-4,4,4/

      data templates(14)%template_num /52/   !  stretch spherical harmonic coeff
      data templates(14)%mapgridlen /8/
      data templates(14)%needext /.false./
      data (templates(14)%mapgrid(j),j=1,8) /4,4,4,1,1,-4,4,-4/

      data templates(15)%template_num /53/   !  strch and rot spher harm coeffs
      data templates(15)%mapgridlen /11/
      data templates(15)%needext /.false./
      data (templates(15)%mapgrid(j),j=1,11) /4,4,4,1,1,-4,4,4,-4,4,-4/

      data templates(16)%template_num /90/     !  space view perspective
      data templates(16)%mapgridlen /21/
      data templates(16)%needext /.false./
      data (templates(16)%mapgrid(j),j=1,21) 
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,4,4,4,4,1,4,4,4,4/

      data templates(17)%template_num /100/    !  triangular grid (icosahedron)
      data templates(17)%mapgridlen /11/
      data templates(17)%needext /.false./
      data (templates(17)%mapgrid(j),j=1,11) /1,1,2,1,-4,4,4,1,1,1,4/

      data templates(18)%template_num /110/ !  equatorial azimuthal equidistant
      data templates(18)%mapgridlen /16/
      data templates(18)%needext /.false./
      data (templates(18)%mapgrid(j),j=1,16) 
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,4,4,1,1/

       data templates(19)%template_num /120/     !  azimuth-range 
       data templates(19)%mapgridlen /7/
       data templates(19)%needext /.true./
       data (templates(19)%mapgrid(j),j=1,7) /4,4,-4,4,4,4,1/

       data templates(20)%template_num /1000/     !  cross section grid 
       data templates(20)%mapgridlen /20/
       data templates(20)%needext /.true./
       data (templates(20)%mapgrid(j),j=1,20) 
     &              /1,1,4,1,4,1,4,4,4,4,-4,4,1,4,4,1,2,1,1,2/

       data templates(21)%template_num /1100/     !  hovmoller diagram grid 
       data templates(21)%mapgridlen /28/
       data templates(21)%needext /.false./
       data (templates(21)%mapgrid(j),j=1,28) 
     &    /1,1,4,1,4,1,4,4,4,4,-4,4,1,-4,4,1,4,1,-4,1,1,-4,2,1,1,1,1,1/

       data templates(22)%template_num /1200/     !  time section grid 
       data templates(22)%mapgridlen /16/
       data templates(22)%needext /.true./
       data (templates(22)%mapgrid(j),j=1,16) 
     &              /4,1,-4,1,1,-4,2,1,1,1,1,1,2,1,1,2/

      data templates(23)%template_num /31/     !  albers equal area
      data templates(23)%mapgridlen /22/
      data templates(23)%needext /.false./
      data (templates(23)%mapgrid(j),j=1,22) 
     &              /1,1,4,1,4,1,4,4,4,-4,4,1,-4,4,4,4,1,1,-4,-4,-4,4/

      contains


         integer function getgridindex(number)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getgridindex
!   prgmmr: gilbert         org: w/np11    date: 2001-06-28
!
! abstract: this function returns the index of specified grid
!   definition template 3.nn (nn=number) in array templates.
!
! program history log:
! 2001-06-28  gilbert
!
! usage:    index=getgridindex(number)
!   input argument list:
!     number   - nn, indicating the number of the grid definition
!                template 3.nn that is being requested.
!
! returns:  index of gdt 3.nn in array templates, if template exists.
!           = -1, otherwise.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: number

           getgridindex=-1

           do j=1,maxtemp
              if (number.eq.templates(j)%template_num) then
                 getgridindex=j
                 return
              endif
           enddo

         end function


         subroutine getgridtemplate(number,nummap,map,needext,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getgridtemplate 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-09
!
! abstract: this subroutine returns grid template information for a 
!   specified grid definition template 3.nn.
!   the number of entries in the template is returned along with a map
!   of the number of octets occupied by each entry.  also, a flag is
!   returned to indicate whether the template would need to be extended.
!
! program history log:
! 2000-05-09  gilbert
!
! usage:    call getgridtemplate(number,nummap,map,needext,iret)
!   input argument list:
!     number   - nn, indicating the number of the grid definition 
!                template 3.nn that is being requested.
!
!   output argument list:      
!     nummap   - number of entries in the template
!     map()    - an array containing the number of octets that each 
!                template entry occupies when packed up into the gds.
!     needext  - logical variable indicating whether the grid defintion
!                template has to be extended.  
!     ierr     - error return code.
!                0 = no error
!                1 = undefine grid template number.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: number
           integer,intent(out) :: nummap,map(*),iret
           logical,intent(out) :: needext

           iret=0

           index=getgridindex(number)

           if (index.ne.-1) then
              nummap=templates(index)%mapgridlen
              needext=templates(index)%needext
              map(1:nummap)=templates(index)%mapgrid(1:nummap)
           else
             nummap=0
             needext=.false.
             print *,'getgridtemplate: grid template ',number,
     &               ' not defined.'
             iret=1
           endif

         end subroutine


         subroutine extgridtemplate(number,list,nummap,map)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    extgridtemplate 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-09
!
! abstract: this subroutine generates the remaining octet map for a 
!   given grid definition template, if required.  some templates can 
!   vary depending on data values given in an earlier part of the 
!   template, and it is necessary to know some of the earlier entry
!   values to generate the full octet map of the template.
!
! program history log:
! 2000-05-09  gilbert
!
! usage:    call extgridtemplate(number,list,nummap,map)
!   input argument list:
!     number   - nn, indicating the number of the grid definition 
!                template 3.nn that is being requested.
!     list()   - the list of values for each entry in 
!                the grid definition template.
!
!   output argument list:      
!     nummap   - number of entries in the template
!     map()    - an array containing the number of octets that each 
!                template entry occupies when packed up into the gds.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: number,list(*)
           integer,intent(out) :: nummap,map(*)

           index=getgridindex(number)
           if (index.eq.-1) return

           if ( .not. templates(index)%needext ) return
           nummap=templates(index)%mapgridlen
           map(1:nummap)=templates(index)%mapgrid(1:nummap)

           if ( number.eq.120 ) then
              n=list(2)
              do i=1,n
                map(nummap+1)=2
                map(nummap+2)=-2
                nummap=nummap+2
              enddo
           elseif ( number.eq.1000 ) then
              n=list(20)
              do i=1,n
                map(nummap+1)=4
                nummap=nummap+1
              enddo
           elseif ( number.eq.1200 ) then
              n=list(16)
              do i=1,n
                map(nummap+1)=4
                nummap=nummap+1
              enddo
           endif

         end subroutine

         integer function getgdtlen(number)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getgdtlen
!   prgmmr: gilbert         org: w/np11    date: 2004-05-11
!
! abstract: this function returns the initial length (number of entries) in
!   the "static" part of specified grid definition template 3.number.
!
! program history log:
! 2004-05-11  gilbert
!
! usage:    call getgdtlen(number)
!   input argument list:
!     number   - nn, indicating the number of the grid definition
!                template 3.nn that is being requested.
!
! returns:     number of entries in the "static" part of gdt 3.number
!              or returns 0, if requested template is not found.
!
! remarks: if user needs the full length of a specific template that
!    contains additional entries based on values set in the "static" part
!    of the gdt, subroutine extgridtemplate can be used.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: number

           getgdtlen=0

           index=getgridindex(number)

           if (index.ne.-1) then
              getgdtlen=templates(index)%mapgridlen
           endif

         end function


      end

