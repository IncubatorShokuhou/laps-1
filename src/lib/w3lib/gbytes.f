      subroutine gbytes(ipackd,iunpkd,noff,nbits,iskip,iter)
c
c this program written by.....
c             dr. robert c. gammill, consultant
c             national center for atmospheric research
c             may 1972
c
c             changes for silicongraphics iris-4d/25
c             silicongraphics 3.3 fortran 77
c             march 1991, russell e. jones
c             national weather service
c
c this is the fortran version of gbytes.
c
c***********************************************************************
c
c subroutine gbytes (ipackd,iunpkd,noff,nbits,iskip,iter)
c
c purpose                to unpack a series of bytes into a target
c                        array.  each unpacked byte is right-justified
c                        in its target word, and the remainder of the
c                        word is zero-filled.
c
c usage                  call gbytes (ipackd,iunpkd,noff,nbits,nskip,
c                                     iter)
c
c arguments
c on input                ipackd
c                           the word or array containing the packed
c                           bytes.
c
c                         iunpkd
c                           the array which will contain the unpacked
c                           bytes.
c
c                         noff
c                           the initial number of bits to skip, left
c                           to right, in 'ipackd' in order to locate
c                           the first byte to unpack.
c
c                        nbits
c                          number of bits in the byte to be unpacked.
c                          maximum of 64 bits on 64 bit machine, 32
c                          bits on 32 bit machine.
c
c                         iskip
c                           the number of bits to skip between each byte
c                           in 'ipackd' in order to locate the next byte
c                           to be unpacked.
c
c                         iter
c                           the number of bytes to be unpacked.
c
c arguments
c on output               iunpkd
c                           contains the requested unpacked bytes.
c***********************************************************************

      integer    ipackd(*)

      integer    iunpkd(*)
      integer    masks(64)
c
      save
c
      data ifirst/1/
      if(ifirst.eq.1) then
         call w3fi01(lw)
         nbitsw = 8 * lw
         jshift = -1 * nint(alog(float(nbitsw)) / alog(2.0))
         masks(1) = 1
         do i=2,nbitsw-1
            masks(i) = 2 * masks(i-1) + 1
         enddo
         masks(nbitsw) = -1
         ifirst = 0
      endif
c
c nbits must be less than or equal to nbitsw                                    
c
      icon   = nbitsw - nbits
      if (icon.lt.0) return
      mask   = masks(nbits)
c
c index tells how many words into the array 'ipackd' the next byte
c appears.         
c
      index  = ishft(noff,jshift)
c
c ii tells how many bits the byte is from the left side of the word.
c
      ii     = mod(noff,nbitsw)
c
c istep is the distance in bits from the start of one byte to the next.
c
      istep  = nbits + iskip      
c
c iwords tells how many words to skip from one byte to the next.                
c
      iwords = istep / nbitsw    
c
c ibits tells how many bits to skip after skipping iwords.                      
c
      ibits  = mod(istep,nbitsw) 
c
      do 10 i = 1,iter
c
c mover specifies how far to the right a byte must be moved in order            
c
c    to be right adjusted.                                                      
c
      mover = icon - ii
c                                                                               
c the byte is split across a word break.                 
c                       
      if (mover.lt.0) then                                                  
        movel   = - mover                                                       
        mover   = nbitsw - movel                                                
        iunpkd(i) = iand(ior(ishft(ipackd(index+1),movel),
     &            ishft(ipackd(index+2),-mover)),mask)
c
c right adjust the byte.
c
      else if (mover.gt.0) then
        iunpkd(i) = iand(ishft(ipackd(index+1),-mover),mask)
c                                             
c the byte is already right adjusted.
c
      else
        iunpkd(i) = iand(ipackd(index+1),mask)
      endif
c                                                                               
c increment ii and index.
c
        ii    = ii + ibits
        index = index + iwords
        if (ii.ge.nbitsw) then
          ii    = ii - nbitsw
          index = index + 1
        endif
c
   10 continue
        return
      end
