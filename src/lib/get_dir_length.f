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

c
        subroutine get_directory_length(c_fname,lenf)
c
c********************************************************************
c
cdoc    this routine takes as input a full path filename and returns an
cdoc    index that points at the end of the directory portion of the pathname.
c       simple-minded algorithm just searches backwards for the first
c       occurance of a `/' (for unix) or a `]' (for vms).
c
c       input/output:
c
c       name            type    i/o     description
c       ----            ---     --      -----------
c       c_fnames        char    i       file name.
c       lenf            i       o       index to end of directory
c
c********************************************************************
c
        character c_fname*(*)
        integer lenf
c
        integer i, strlen
c
c****************************
c
        strlen = len(c_fname)
c
        i = strlen
        do while (i .gt. 0)
        if( (c_fname(i:i) .ne. ']')
     1.and. (c_fname(i:i) .ne. '/') )then
           i = i-1
        else
           goto 100
        endif
        enddo
c
100     lenf = i
c
        return
        end

c
        subroutine get_time_length(c_fname,lenf)
c
c********************************************************************
c
cdoc    this routine takes as input a full path filename and returns an
cdoc    index that points at the end of the filetime portion of the pathname.
c       simple-minded algorithm just searches backwards for the first
c       occurance of a `.'.
c
c       input/output:
c
c       name            type    i/o     description
c       ----            ---     --      -----------
c       c_fnames        char    i       file name.
c       lenf            i       o       index to end of filetime
c
c********************************************************************
c
        character c_fname*(*)
        integer lenf
c
        integer i, strlen
c
c****************************
c
        call s_len(c_fname,strlen)
c
        i = strlen
        do while (i .gt. 0)
        if (c_fname(i:i) .ne. '.')then
           i = i-1
        else
           goto 100
        endif
        enddo
c
100     lenf = i - 1
c
        return
        end

c
c#####################################################
c
        subroutine s_len(string,s_length)

c*********************************************************************
c
cdoc    this routine receives a fortran string, and
cdoc    returns the number of characters in the string
cdoc    before the first "space" is encountered.  
c       it considers ascii characters 33 to 126 to be valid
c       characters, and ascii 0 to 32, and 127 to be "space"
c       characters.
c
c       name            type      i/o     description
c       ----            ---       --      -----------
c       string          char       i       string
c       s_length        integer  o       valid number characters
c                                            in string

        implicit none

        character*(*)   string
        integer       s_length, i, len_str, aval
        logical         space

        space = .false.
        i = 1
        len_str = len(string)
        s_length = len_str      !default, used if string
                                !  is full, with no "spaces"
        do while ((i .le. len_str) .and. (.not. space))
          aval = ichar(string(i:i))
          if ((aval .lt. 33) .or. (aval .gt. 126)) then
            s_length = i - 1
            space = .true.
          endif
          i = i + 1
        enddo

        return
        end

c####################################################################
c
        subroutine get_fname_length(c_fname,lenf)
c
c********************************************************************
c
cdoc    this routine takes as input a full path filename and returns
cdoc    the length of the filename portion of the path.
c       simple-minded algorithm just searches backwards for the first
c       occurance of a `/'. the number of characters searched up to
c       that point indicates the filename length.
c
c       input/output:
c
c       name            type    i/o     description
c       ----            ---     --      -----------
c       c_fnames        char    i       file name.
c       lenf            i       o       length of filename part of c_fnames
c
c********************************************************************
c
        character c_fname*(*)
        integer lenf
c
        integer i, strlen, i_char_len
c
c****************************
c
        i_char_len = len(c_fname)
        call s_len(c_fname,strlen)
c
        i = i_char_len
        do while (i .gt. 0)
        if(c_fname(i:i) .ne. '/')then
           i = i-1
        else
           goto 100
        endif
        enddo
c
100     lenf = strlen - i
c
        return
        end
c
c#####################################################################
c
        subroutine get_filetime_length(c_fname,lent)
c
c*********************************************************************
c
cdoc    this routine takes as input a full path filename and returns
cdoc    the length of the file time length portion of the path.
c       simple-minded algorithm uses other simple minded algorithms
c       found in this file to determine the length of the time portion
c       of the character string.
c

        character*(*) c_fname
        integer lend,lent,lenf

        call get_directory_length(c_fname,lend)
        call s_len(c_fname,len_fname)

        lenf = len_fname - lend

        if(lenf .ge. 13 .and. c_fname(lend+9:lend+9).eq.'_')then
           lent=13

        elseif(lenf .ge. 9)then !assume 9 chars for time portion of filename
           lent=9

        else
           lent = lenf

        endif

c       write(6,*)'time portion of string = ',c_fname(lend+1:lend+lent)       

        return
        end

c
c#####################################################################
c
        subroutine get_filetime_type(c_fname,c20_type,leni,lent)
c
c*********************************************************************
c
cdoc    this routine takes as input a full path filename and returns
cdoc    the length and type of the file time length portion of the path.
c       simple-minded algorithm uses other simple minded algorithms
c       found in this file to determine the length and type of the time portion
c       of the character string.
c
c       1998          steve albers

        character*(*) c_fname
        character*20 c20_type

        integer lend ! directory length including the last 'slash'. 
        integer lenf ! length of the filename excluding the path.
        integer lent ! length of the filetime portion.
        integer leni ! length of the initial non-filetime portion.
                       ! this is often but not always the directory length.

        call filter_non_numeric_fnames(c_fname,1,num_out,1
     1                                    ,istatus)

c initialize
        c20_type = 'unknown'

        if(num_out.eq.1)then

           call s_len(c_fname,len_fname)
           call get_directory_length(c_fname,lend)
           lenf = len_fname - lend

           if(lenf .eq. 20)then
            if(c_fname(lend+14:lend+14) .eq. '_')then
               c20_type = 'yyyyjjjhhmmss'                         ! rsa radar type
               leni = lend
               lent=13
               return
            endif
           endif

           if(lenf .eq. 13)then
            if(c_fname(lend+1:lend+5) .eq. 'raob.')then
               c20_type = 'yymmddhh'                              ! afwa raob
               leni = lend+5
               lent=8
               return
c
c - repositioned to "lenf .ge. 13" switch below. j.smart 3-20-00 
c           elseif(c_fname(lend+1:lend+2) .eq. 'nf'  )            !.or.
c     +             c_fname(lend+1:lend+2) .eq. 're' )then         ! taiwan fa model
c              c20_type = 'ymmddhh'
c              leni = lend+2
c              lent = 7
c              return

            endif
           endif

           if(lenf .ge. 13)then
            if(c_fname(lend+9:lend+9) .eq. '_')then
               c20_type = 'yyyymmdd_hhmm'                         ! wfo type
               leni = lend
               lent=13
               return
            elseif(c_fname(lend+1:lend+2) .eq. 'nf'.or.
     +             c_fname(lend+1:lend+2) .eq. 're'.or.
     +             c_fname(lend+1:lend+2) .eq. 'gb'.or.
     +             c_fname(lend+1:lend+2) .eq. 'sb'.or.
     +             c_fname(lend+1:lend+2) .eq. 'gs')then
               c20_type = 'yyyymmddhh'                            !taiwan/cwb nfs, gfs, or sb (tropical cyclone) model
               leni = lend+2
               lent = 10
               return
            endif
           endif

           if(lenf .eq. 16)then
            if(c_fname(lend+1:lend+4) .eq. 'temp')then
               c20_type = 'yymmddhh'                              ! cwb raob
               leni = lend+4
               lent=8
               return
            endif
           endif

           if(lenf .ge. 9)then ! assume 9 chars for time portion of filename
              c20_type = 'yyjjjhhmm'                             ! nimbus/laps
              leni = lend
              lent=9
              return
           endif

           c20_type = 'unknown'
           leni = lend
           lent = lenf

        endif

!       write(6,*)'time portion of string = ',c20_type
!    1                                       ,c_fname(leni+1:leni+lent)
        
        return
        end
c
c#####################################################
c
        subroutine filter_string(string)

c       it considers ascii characters 33 to 126 to be valid
c       characters, and ascii 0 to 32, and 127 to be "space"
c       characters.

c*********************************************************************
c
cdoc    this routine filters a fortran string of unprintable characters
c
c       name            type      i/o     description
c       ----            ---       --      -----------
c       string          char       i       string

        implicit none

        character*(*)   string
        integer       i, len_str, aval
        logical         space

        space = .false.
        len_str = len(string)

        do i = 1,len_str
          aval = ichar(string(i:i))
          if ((aval .lt. 33) .or. (aval .gt. 126)) then
            space = .true.
            string(i:i) = ' '
          endif
        enddo

        return
        end


        function l_string_contains(string,substring,istatus)

cdoc    returns boolean on whether 'string' contains 'substring'
cdoc    the 'string' and 'substring' should be free of blanks.

        logical l_string_contains

        character*(*) string
        character*(*) substring
 
        l_string_contains = .false.

        call s_len(string,len1)

        call s_len(substring,len2)

        if(len1 .ge. len2)then
            i_search_end = len1 - len2 + 1

            do i = 1,i_search_end
                if(string(i:i+len2-1) .eq. substring(1:len2))then
                    l_string_contains = .true.
                endif
            enddo ! i

            istatus = 1
            return

        else
            istatus = 0
            return

        endif

        end


        function l_parse(string1,string2)

cdoc    returns boolean on whether 'string1' contains 'string2'
cdoc    similar to 'l_string_contains' except that blanks are allowed

        logical l_parse

        character*(*) string1,string2

!       integer slen1,slen2

        len1 = len(string1)
        len2 = len(string2)

!       call s_len(string1,slen1)
!       call s_len(string2,slen2)

        l_parse = .false.

        if(len2 .gt. len1)return

        i_offset_max = len1-len2

        do i = 0,i_offset_max
            if(string1(i+1:i+len2) .eq. string2(1:len2))then
                l_parse = .true.
            endif ! match is found
        enddo ! i             

        return
        end        


        subroutine s_len2(string,len_string)

cdoc    this routine finds the length of the string counting intermediate
cdoc    blanks

        character*(*) string

        len1 = len(string)

        len_string = 0

        do i = 1,len1
            call s_len(string(i:i),len2)
            if(len2 .eq. 1)len_string = i
        enddo ! i

        return
        end

