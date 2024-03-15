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
      subroutine rcp_to_remote(i4time,ext
     1                   ,destination_node,destination_dir)

!     you guessed it, copy forbidden data outside the lab

      character*150 dir
      character*(*) destination_dir
      character*(*) destination_node
      character*200 command
!     character*200 remote_command

      character*9 gtime
      character*31 ext,ext_in

      character*91 full_fname
      integer len_fname

      write(6,*)' subroutine rcp_to_remote'

      call s_len (destination_node,len_node)

      call s_len (destination_dir,len_dir)

!     get ascii 9 time
      call make_fnam_lp(i4time,gtime,istatus)

!     get directory to copy file from
      call downcase(ext,ext_in)
      call s_len(ext_in,len_ext)
      dir='../lapsprd/'//ext_in(1:len_ext)//'/'

!     get full file name of file being copied
      call cvt_fname_data(dir,gtime,ext_in,full_fname,istatus)
      call s_len(full_fname,len_fname)

      write(6,*)full_fname

!     generate full command

!     option 1 - use simple rcp command
      command = 'gzip '//full_fname(1:len_fname)//';'
     1             //' rcp '//full_fname(1:len_fname)//'* '
     1             //destination_node(1:len_node)
     1             //destination_dir(1:len_dir)

!     option 2 - use dan's tar command
!     remote_command = '(cd '//destination_dir(1:len_dir)
!    1                       //'; tar -xf -)'
!     command = 'gzip '//full_fname(1:len_fname)
!    1                 //' | 'tar -cf | rsh '
!    1                 //destination_node(1:len_node)
!    1                 //remote_command 

      write(6,*)' calling system with this command:'
      write(6,*)command

      call system(command)

      return
      end
