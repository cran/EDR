      subroutine edrstp3a(d,n,dp1,wij,kksi,kksii,mat,wi,s,work,
     1                     iwork)
C
C        length(work) : 10*(3*dp1+imax(n,6*dp1))
C        length(iwork): 8*dp1
C
      integer n,d,dp1,iwork(*)
      double precision wij(n,d),kksi(*),kksii(n),wi(n),mat(dp1,n),
     1       s(dp1),work(*)
      integer i,j,k,ll,lwork,info
      double precision u,vt,omkksii
      external dgesdd
      DO i=1,n
         call extrdist(n,kksi,i,kksii)
C        this provides distances from point X_i
         ll=0
         DO j=1,n
            if(kksii(j).gt.1.d0) CYCLE
            omkksii=1.d0-kksii(j)
            ll=ll+1
            mat(1,ll)=1.d0
            DO k=1,d
               mat(k+1,ll)=(wij(j,k)-wij(i,k))*omkksii
            END DO
         END DO
         if(ll.le.dp1) THEN
            wi(i) = 0.d0
         ELSE
            lwork=3*dp1+max(ll,6*dp1)
         call dgesdd('N',dp1,ll,mat,dp1,s,u,1,vt,1,work,10*lwork,
     1           iwork,info)
         if(info.ne.0) call intpr("info",4,info,1)
            wi(i)=s(dp1)/s(1)
         ENDIF
      END DO
      RETURN
      END
      subroutine edrstp3b(d,n,dp1,wij,kksi,y,kksii,mat,s,u,vt,work,
     1                     iwork,fx,fw,lll,yw)
C
C        length(work) : 5*(3*dp1*dp1+max(n,4*dp1*dp1+4*dp1))
C        length(iwork): 8*dp1
C
      implicit logical (a-z)
      integer n,d,dp1,iwork(*)
      double precision wij(n,d),kksi(*),kksii(n),mat(dp1,n),s(dp1),
     1       work(*),fx(n,d),lll,y(n),yw(n),fw(n)
      integer i,j,k,m,ll,lwork,info
      double precision u(dp1,dp1),vt(dp1,n),omkksii,z,z1
      external dgesdd
      lll=0.d0
      DO i=1,n
         call extrdist(n,kksi,i,kksii)
C        this provides distances from point X_i
         ll=0
         DO j=1,n
            if(kksii(j).gt.1.d0) CYCLE
            omkksii=1.d0-kksii(j)
            ll=ll+1
            mat(1,ll)=omkksii
            DO k=1,d
               mat(k+1,ll)=(wij(j,k)-wij(i,k))*omkksii
            END DO
            yw(ll)=y(j)*omkksii
            lll=lll+omkksii*omkksii
         END DO
         if(ll.gt.dp1) THEN
            lwork=3*dp1*dp1+max(n,4*dp1*dp1+4*dp1)
         call dgesdd('S',dp1,ll,mat,dp1,s,u,dp1,vt,dp1,work,5*lwork,
     1           iwork,info)
            if(info.ne.0) call intpr("info",4,info,1)
            if(s(dp1).le.1d-5) CYCLE 
            DO m=1,d
               z=0.d0
               DO j=1,ll
                  z1=0.d0
                  DO k=1,dp1
                     z1=z1+u(m+1,k)*vt(k,j)/s(k)
                  END DO
                  z=z+z1*yw(j)
               END DO
               fx(i,m)=z
            END DO
            z=0.d0
            DO j=1,ll
               z1=0.d0
               DO k=1,dp1
                  z1=z1+u(1,k)*vt(k,j)/s(k)
               END DO
               z=z+z1*yw(j)
            END DO
            fw(i)=z
         ELSE
            fw(i)=y(i)
         ENDIF
      END DO
      RETURN
      END
      subroutine extrdist(n,dist,j,dj)
      integer n,j
      double precision dist(*),dj(n)
      integer i
      IF(j.ne.1) THEN
         DO i=1,j-1
            dj(i)=dist(j-i*(i+1)/2+(i-1)*n)
         END DO
      END IF
      dj(j)=0.d0
      IF(j.ne.n) THEN
         DO i=j+1,n
            dj(i)=dist((j-1)*(2*n-j)/2+i-j)
         END DO
      END IF
      RETURN
      END
      subroutine edrstp3c(d,n,nest,dp1,x,xest,kksi,y,mat,s,u,vt,work,
     1                     iwork,fw,yw)
C
C        length(work) : 5*(3*dp1*dp1+max(n,4*dp1*dp1+4*dp1))
C        length(iwork): 8*dp1
C
      implicit logical (a-z)
      integer n,nest,d,dp1,iwork(*)
      double precision x(n,d),xest(nest,d),kksi(n,nest),mat(dp1,n),
     1       s(dp1),work(*),y(n),yw(n),fw(nest)
      integer i,j,k,ll,lwork,info
      double precision u(dp1,dp1),vt(dp1,n),omkksii,z,z1,kksij,swj,
     1      swjy,kksimin
      external dgesdd
      DO i=1,nest
         ll=0
         swj=0.d0
         swjy=0.d0
         DO j=1,n
            kksij=kksi(j,i)
            if(kksij.gt.1.d0) CYCLE
            omkksii=1.d0-kksij
            ll=ll+1
            mat(1,ll)=omkksii
            swj=swj+omkksii
            DO k=1,d
               mat(k+1,ll)=(x(j,k)-xest(i,k))*omkksii
            END DO
            yw(ll)=y(j)*omkksii
            swjy=swjy+yw(ll)
         END DO
         if(ll.gt.dp1) THEN
            lwork=3*dp1*dp1+max(n,4*dp1*dp1+4*dp1)
         call dgesdd('S',dp1,ll,mat,dp1,s,u,dp1,vt,dp1,work,5*lwork,
     1           iwork,info)
            if(info.ne.0.or.s(dp1).le.1d-5)  THEN
               IF(swj.gt.0.d0) THEN
C  use kernel estimate
                  fw(i)=swjy/swj
               ELSE
C  search for nearest element
                  ll=1
                  kksimin=kksi(1,i)
                  DO j=2,n
                     if(kksi(j,i).lt.kksimin) THEN
                        kksimin=kksi(j,i)
                        ll=j
                     END IF
                  END DO
                  fw(i)=y(ll)
               ENDIF
            ELSE
               z=0.d0
               DO j=1,ll
                  z1=0.d0
                  DO k=1,dp1
                     z1=z1+u(1,k)*vt(k,j)/s(k)
                  END DO
                  z=z+z1*yw(j)
               END DO
               fw(i)=z
            END IF
         ELSE IF(swj.gt.0.d0) THEN
C  use kernel estimate
            fw(i)=swjy/swj
         ELSE
C  search for nearest element
            ll=1
            kksimin=kksi(1,i)
            DO j=2,n
               if(kksi(j,i).lt.kksimin) THEN
                  kksimin=kksi(j,i)
                  ll=j
               END IF
            END DO
            fw(i)=y(ll)
         ENDIF
      END DO
      RETURN
      END
