      subroutine getbit(ibm,ibs,ids,len,mg,g,ground,gmin,gmax,nbit)
c$$$  subprogram documentation block
c
c subprogram:    getbit      compute number of bits and round field.
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: the number of bits required to pack a given field
c   for particular binary and decimal scalings is computed.
c   the field is rounded off to the decimal scaling for packing.
c   the minimum and maximum rounded field values are also returned.
c   grib bitmap masking for valid data is optionally used.
c
c program history log:
c   96-09-16  iredell
c
c usage:    call gtbits(ibm,ibs,ids,len,mg,g,gmin,gmax,nbit)
c   input argument list:
c     ibm      - integer bitmap flag (=0 for no bitmap)
c     ibs      - integer binary scaling
c                (e.g. ibs=3 to round field to nearest eighth value)
c     ids      - integer decimal scaling
c                (e.g. ids=3 to round field to nearest milli-value)
c                (note that ids and ibs can both be nonzero,
c                 e.g. ids=1 and ibs=1 rounds to the nearest twentieth)
c     len      - integer length of the field and bitmap
c     mg       - integer (len) bitmap if ibm=1 (0 to skip, 1 to keep)
c     g        - real (len) field
c
c   output argument list:
c     ground   - real (len) field rounded to decimal and binary scaling
c                (set to zero where bitmap is 0 if ibm=1)
c     gmin     - real minimum valid rounded field value
c     gmax     - real maximum valid rounded field value
c     nbit     - integer number of bits to pack
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension mg(len),g(len),ground(len)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  round field and determine extremes where bitmap is on
      s=2.**ibs*10.**ids
      if(ibm.eq.0) then
        ground(1)=nint(g(1)*s)/s
        gmax=ground(1)
        gmin=ground(1)
        do i=2,len
          ground(i)=nint(g(i)*s)/s
          gmax=max(gmax,ground(i))
          gmin=min(gmin,ground(i))
        enddo
      else
        i1=1
        dowhile(i1.le.len.and.mg(i1).eq.0)
          i1=i1+1
        enddo
        if(i1.le.len) then
          do i=1,i1-1
            ground(i)=0.
          enddo
          ground(i1)=nint(g(i1)*s)/s
          gmax=ground(i1)
          gmin=ground(i1)
          do i=i1+1,len
            if(mg(i).ne.0) then
              ground(i)=nint(g(i)*s)/s
              gmax=max(gmax,ground(i))
              gmin=min(gmin,ground(i))
            else
              ground(i)=0.
            endif
          enddo
        else
          do i=1,len
            ground(i)=0.
          enddo
          gmax=0.
          gmin=0.
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute number of bits
      nbit=log((gmax-gmin)*s+0.9)/log(2.)+1.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
