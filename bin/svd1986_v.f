      module svd1986
      use ANALYSIS
      implicit none

      contains

C    I add a subroutine that does (most) of the
C    pre- and post- processing for the SVD when
C    in use for:
C
C      Creating graphs from pairwise distances
C
      SUBROUTINE PD2G(PD,Ninputs,G)
      implicit none
      integer,intent(in) :: Ninputs
      real(dp),dimension(Ninputs,Ninputs),intent(in) :: PD
      real(dp),dimension(Ninputs,2),intent(out) :: G
      real(dp),dimension(Ninputs,Ninputs) :: A, V, counts
      real(dp),dimension(Ninputs) :: w
      real(dp) :: maxW1, maxW2
      integer :: indexMaxW1, indexMaxW2
      integer :: i,j

C     Assume that PD is already centered

C     Normalize A
      counts = - 1.0d0 / Ninputs
      do i = 1, Ninputs
          counts(i,i) = counts(i,i) + 1.0d0
      end do
      A = -0.5d0 * matmul(counts,matmul(PD**2,counts))

C     Get the SVD decomposition. If A is
C     symmetric then the left and right
C     matrices should be the same. A should
C     by symmetric because PD is as well.
      call svdcmp(A,Ninputs,Ninputs,
     *              Ninputs,Ninputs,w,V)

C     Find the two largest eigenvalues
      maxW1 = 0.0d0
      maxW2 = 0.0d0
      indexMaxW1 = 1
      indexMaxW2 = 1

      do i = 1, Ninputs
          if (w(i) > maxW1) then
              maxW1 = w(i)
              indexMaxW1 = i
              maxW2 = maxW1
              indexMaxW2 = indexMaxW1
          elseif (w(i) > maxW2) then
              maxW2 = w(i)
              indexMaxW2 = i
          end if
      end do

C     Output coordinates for the two largest eigenvalues
      G(1:Ninputs,1) = V(1:Ninputs,indexMaxW1) * sqrt(maxW1)
      G(2:Ninputs,2) = V(1:Ninputs,indexMaxW2) * sqrt(maxW2)

      return
      END SUBROUTINE PD2G

c  SVDcmp subroutine
c    Given a matrix (1:m,1:n) with physical dimensions mp by np,
c    this routine computes its singular value decomposition,
c    A = U W VT.  The matrix U replaces a on output.  The diagonal
c    matrix of singular values W is output as a vector w(1:n)
c    The matrix V (not the transpose VT) is the output as v(1:n,1:n) 
c
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      integer,intent(in) :: m,mp,n,np
      real(dp),intent(inout) :: a(mp,np)
      real(dp),intent(out) :: v(np,np),w(np)
      integer,parameter :: NMAX = 10000
CU    USES pythag
      integer :: i,its,j,jj,k,l,nm
!     REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      real(dp) :: anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX)
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=pythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.


C
C

      FUNCTION pythag(a,b)
      real(dp) :: a,b,pythag
      real(dp) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          pythag=0.0d0
        else
          pythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      end module svd1986
