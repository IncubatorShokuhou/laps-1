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

        subroutine      upcase(input,output)

!       this routine only handles strings up to 500 characters long.

        character*(*)   input,
     1          output
cc        character*500   string

        integer       nchar,
!       1               lnblnk,
!     1          l,
     1          len,
     1          chr,
     1          i

cc        string=input

!       nchar=lnblnk(string)

cc        if(string(500:500) .ne. ' ')then
cc            write(6,*)'string truncated to 500 characters.'
cc        endif

cc        do i = 500,1,-1
cc            if(string(i:i) .ne. ' ')then
cc                nchar = i
cc                go to 10
cc            endif
cc        enddo

10      continue

        len_in =len(input)
        len_out=len(output)
        nchar=min(len_in,len_out)

        do i=1,nchar
                chr=ichar(input(i:i))
                if (chr .ge. 97 .and. chr .le. 122) chr=chr-32
                output(i:i)=char(chr)
        enddo

        do i=nchar+1,len_out
          output(i:i)=' '
        enddo

        return

        end

