      module drstemplates
!$$$  subprogram documentation block
!                .      .    .                                       .
! module:    drstemplates 
!   prgmmr: gilbert         org: w/np11    date: 2001-04-03
!
! abstract: this fortran module contains info on all the available 
!   grib2 data representation templates used in section 5 (drs).
!   each template has three parts: the number of entries in the template
!   (mapgridlen);  a map of the template (mapgrid), which contains the
!   number of octets in which to pack each of the template values; and
!   a logical value (needext) that indicates whether the template needs 
!   to be extended.  in some cases the number of entries in a template 
!   can vary depending upon values specified in the "static" part of 
!   the template.  ( see template 5.1 as an example )
!
!   this module also contains two subroutines.  subroutine getdrstemplate
!   returns the octet map for a specified template number, and
!   subroutine extdrstemplate will calculate the extended octet map
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
! 2000-05-11  gilbert
! 2002-12-11  gilbert - added templates for jpeg2000 and png encoding
!
! usage:    use drstemplates
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      integer,parameter :: maxlen=200,maxtemp=9

      type drstemplate
          integer :: template_num
          integer :: mapdrslen
          integer,dimension(maxlen) :: mapdrs
          logical :: needext
      end type drstemplate

      type(drstemplate),dimension(maxtemp) :: templates

      data templates(1)%template_num /0/     !  simple packing
      data templates(1)%mapdrslen /5/ 
      data templates(1)%needext /.false./
      data (templates(1)%mapdrs(j),j=1,5) 
     &                             /4,-2,-2,1,1/

      data templates(2)%template_num /2/     !  complex packing
      data templates(2)%mapdrslen /16/
      data templates(2)%needext /.false./
      data (templates(2)%mapdrs(j),j=1,16)
     &                        /4,-2,-2,1,1,1,1,4,4,4,1,1,4,1,4,1/

      data templates(3)%template_num /3/     !  complex packing - spatial diff
      data templates(3)%mapdrslen /18/
      data templates(3)%needext /.false./
      data (templates(3)%mapdrs(j),j=1,18)
     &                        /4,-2,-2,1,1,1,1,4,4,4,1,1,4,1,4,1,1,1/

      data templates(4)%template_num /50/     !  simple packing - spectral data
      data templates(4)%mapdrslen /5/
      data templates(4)%needext /.false./
      data (templates(4)%mapdrs(j),j=1,5)
     &                         /4,-2,-2,1,4/

      data templates(5)%template_num /51/    !  complex packing - spectral data
      data templates(5)%mapdrslen /10/
      data templates(5)%needext /.false./
      data (templates(5)%mapdrs(j),j=1,10)
     &                         /4,-2,-2,1,-4,2,2,2,4,1/

      data templates(6)%template_num /40000/     !  jpeg2000 encoding
      data templates(6)%mapdrslen /7/ 
      data templates(6)%needext /.false./
      data (templates(6)%mapdrs(j),j=1,7) 
     &                             /4,-2,-2,1,1,1,1/

      data templates(7)%template_num /40010/     !  png encoding
      data templates(7)%mapdrslen /5/ 
      data templates(7)%needext /.false./
      data (templates(7)%mapdrs(j),j=1,5) 
     &                             /4,-2,-2,1,1/

      data templates(8)%template_num /40/     !  jpeg2000 encoding
      data templates(8)%mapdrslen /7/ 
      data templates(8)%needext /.false./
      data (templates(8)%mapdrs(j),j=1,7) 
     &                             /4,-2,-2,1,1,1,1/

      data templates(9)%template_num /41/     !  png encoding
      data templates(9)%mapdrslen /5/ 
      data templates(9)%needext /.false./
      data (templates(9)%mapdrs(j),j=1,5) 
     &                             /4,-2,-2,1,1/

!      data templates(5)%template_num /1/      !  simple packing - matrix
!      data templates(5)%mapdrslen /15/ 
!      data templates(5)%needext /.true./
!      data (templates(5)%mapdrs(j),j=1,15)
!     &                        /4,-2,-2,1,1,1,4,2,2,1,1,1,1,1,1/


      contains

         integer function getdrsindex(number)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getdrsindex 
!   prgmmr: gilbert         org: w/np11    date: 2001-06-28
!
! abstract: this function returns the index of specified data 
!   representation template 5.nn (nn=number) in array templates.
!
! program history log:
! 2001-06-28  gilbert
!
! usage:    index=getdrsindex(number)
!   input argument list:
!     number   - nn, indicating the number of the data representation 
!                template 5.nn that is being requested.
!
! returns:  index of drt 5.nn in array templates, if template exists.
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

           getdrsindex=-1

           do j=1,maxtemp
              if (number.eq.templates(j)%template_num) then
                 getdrsindex=j
                 return
              endif
           enddo

         end function


         subroutine getdrstemplate(number,nummap,map,needext,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getdrstemplate 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-11
!
! abstract: this subroutine returns drs template information for a 
!   specified data representation template 5.nn.
!   the number of entries in the template is returned along with a map
!   of the number of octets occupied by each entry.  also, a flag is
!   returned to indicate whether the template would need to be extended.
!
! program history log:
! 2000-05-11  gilbert
!
! usage:    call getdrstemplate(number,nummap,map,needext,iret)
!   input argument list:
!     number   - nn, indicating the number of the data representation 
!                template 5.nn that is being requested.
!
!   output argument list:      
!     nummap   - number of entries in the template
!     map()    - an array containing the number of octets that each 
!                template entry occupies when packed up into the drs.
!     needext  - logical variable indicating whether the data representation
!                template has to be extended.  
!     ierr     - error return code.
!                0 = no error
!                1 = undefined data representation template number.
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

           index=getdrsindex(number)

           if (index.ne.-1) then
              nummap=templates(index)%mapdrslen
              needext=templates(index)%needext
              map(1:nummap)=templates(index)%mapdrs(1:nummap)
           else
             nummap=0
             needext=.false.
             print *,'getdrstemplate: drs template ',number,
     &               ' not defined.'
             iret=1
           endif

         end subroutine

         subroutine extdrstemplate(number,list,nummap,map)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    extdrstemplate 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-11
!
! abstract: this subroutine generates the remaining octet map for a
!   given data representation template, if required.  some templates can
!   vary depending on data values given in an earlier part of the 
!   template, and it is necessary to know some of the earlier entry
!   values to generate the full octet map of the template.
!
! program history log:
! 2000-05-11  gilbert
!
! usage:    call extdrstemplate(number,list,nummap,map)
!   input argument list:
!     number   - nn, indicating the number of the data representation 
!                template 5.nn that is being requested.
!     list()   - the list of values for each entry in the 
!                the data representation template 5.nn.
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

           index=getdrsindex(number)
           if (index.eq.-1) return

           if ( .not. templates(index)%needext ) return
           nummap=templates(index)%mapdrslen
           map(1:nummap)=templates(index)%mapdrs(1:nummap)

           if ( number.eq.1 ) then
              n=list(11)+list(13)
              do i=1,n
                map(nummap+i)=4
              enddo
              nummap=nummap+n
           endif

         end subroutine

      end module



