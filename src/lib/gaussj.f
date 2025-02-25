
      subroutine gaussj(a,n,np,b,m,mp,ierr)

c  purpose: solution of the system of linear equations ax = b by
c     gauss-jordan elimination, where a is a matrix of order n and b is
c     an n x m matrix.  on output a is replaced by its matrix inverse
c     and b is preplaced by the corresponding set of solution vectors.

c  source: w.h. press et al, "numerical recipes," 1989, p. 28.

c  modifications: 
c     1. double  precision.
c     2. error parameter ierr included.  0 = no error. 1 = singular 
c        matrix encountered; no inverse is returned.

c  prepared by j. applequist, 8/17/91.

c        set largest anticipated value of n.

      dimension a(np,np),b(np,mp),ipiv(n),indxr(n),indxc(n)

      ierr=0
      do 11 j=1,n
      ipiv(j)=0
 11   continue
      do 22 i=1,n
      big=0.d0
      do 13 j=1,n
      if (ipiv(j).ne.1) then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
      if (abs(a(j,k)).ge.big) then
      big=abs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      ierr=1
      return
      endif
 12   continue
      endif
 13   continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
 14   continue
      do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
 15   continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.d0) then
      ierr=1
      return
      endif
      pivinv=1.d0/a(icol,icol)
      a(icol,icol)=1.d0
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
 16   continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
 17   continue
      do 21 ll=1,n
      if (ll.ne.icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.d0
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
 18   continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
 19   continue
      endif
 21   continue
 22   continue
      do 24 l=n,1,-1
      if (indxr(l).ne.indxc(l)) then
      do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
 23   continue
      endif
 24   continue
      return
      end
