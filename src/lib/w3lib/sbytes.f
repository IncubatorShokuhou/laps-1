      subroutine sbytes(ipackd,iunpkd,noff,nbits,iskip,iter)                            
c this program written by.....                                                  
c             dr. robert c. gammill, consultant                                 
c             national center for atmospheric research                          
c             july 1972                                                         
c this is the fortran versions of sbytes. 
c
c             changes for silicongraphics iris-4d/25  
c             silicongraphics 3.3 fortran 77
c             march 1991  russell e. jones
c             national weather service
c
c***********************************************************************
c
c subroutine sbyte (ipackd,iunpkd,noff,nbits,iskip,iter)
c
c purpose                given a byte, right-justified in a word, to
c                        pack the byte into a target word or array.
c                        bits surrounding the byte in the target
c                        area are unchanged.
c
c usage                  call sbyte (ipackd,iunpkd,noff,nbits)
c
c arguments
c on input               ipackd
c                          the word or array which will contain the
c                          packed byte.  byte may cross word boundaries.
c
c                        iunpkd
c                          the word containing the right-justified byte
c                          to be packed.
c
c                        noff
c                          the number of bits to skip, left to right,
c                          in 'ipackd' in order to locate where the
c                          byte is to be packed.
c
c                        nbits
c                          number of bits in the byte to be packed.
c                          maximum of 64 bits on 64 bit machine, 32
c                          bits on 32 bit machine.
c
c                         iskip
c                           the number of bits to skip between each byte
c                           in 'iunpkd' in order to locate the next byte
c                           to be packed.
c
c                         iter
c                           the number of bytes to be packed.
c
c on output              ipackd
c                          word or consecutive words containing the
c                          requested byte.
c
c***********************************************************************

      integer     iunpkd(*)
      integer     ipackd(*)
      integer     masks(64)
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
      icon = nbitsw - nbits
      if (icon.lt.0) return   
      mask   = masks(nbits)
c
c index tells how many words into iout the next byte is to be stored.           
c
      index  = ishft(noff,jshift)  
c
c ii tells how many bits in from the left side of the word to store it.         
c
      ii     = mod(noff,nbitsw)
c
c istep is the distance iunpkd bits from one byte position to the next.             
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
        j     = iand(mask,iunpkd(i))    
        movel = icon - ii         
c                                                                               
c byte is to be stored in middle of word.  shift left.                          
c
        if (movel.gt.0) then
          msk           = ishft(mask,movel)  
          ipackd(index+1) = ior(iand(not(msk),ipackd(index+1)),
     &    ishft(j,movel))
c                                                                               
c the byte is to be split across a word break.                                  
c
        else if (movel.lt.0) then
          msk           = masks(nbits+movel)      
          ipackd(index+1) = ior(iand(not(msk),ipackd(index+1)),
     &    ishft(j,movel))  
          itemp         = iand(masks(nbitsw+movel),ipackd(index+2))
          ipackd(index+2) = ior(itemp,ishft(j,nbitsw+movel))
c             
c byte is to be stored right-adjusted.                                          
c
        else
          ipackd(index+1) = ior(iand(not(mask),ipackd(index+1)),j)
        endif
c     
        ii    = ii + ibits 
        index = index + iwords    
        if (ii.ge.nbitsw) then
          ii    = ii - nbitsw 
          index = index + 1
        endif
c
10    continue
c
      return
      end
