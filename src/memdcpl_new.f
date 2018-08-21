      subroutine  memdcpl_new(isw1,isw2,fus,minll,maxll)

      ! dimension information stored in a header file
      include 'dimension.h'

      common/eig1/nord1,jcom1,lord1,wcom1,qbar1,cgp1,avert1,ahor1,phis1
     1          ,eif1(222,6)
      common/eig2/nord2,jcom2,lord2,wcom2,qbar2,cgp2,avert2,ahor2,phis2
     1          ,eif2(222,6)

      common/nond/rn,wn,vn,gn,rhobar
!      common/int/rot,asp
      common/int/rot
      common/int0/asp0,asp2
!      common/matr/w0,z
      common/matr0/z0,z1,z2

      complex    iunt,w0
      complex*16 zr(-ML:ML,-ML:ML),z(-ML:ML,-ML:ML),tt(0:MLL,-MLL:MLL),t
      complex*16 zr0(-ML:ML,-ML:ML),z0(-ML:ML,-ML:ML),tt0(0:MLL,-MLL:MLL)
      complex*16 zr1(-ML:ML,-ML:ML),z1(-ML:ML,-ML:ML)
      complex*16 z2(-ML:ML,-ML:ML),tt2(0:MLL,-MLL:MLL),t2,t0


      double precision rot,asp(0:MLL,2*MLL+1),sj(-ML:ML,-ML:ML)
      double precision asp0(0:MLL,2*MLL+1),asp2(0:MLL,2*MLL+1)


      data iunt/(0.,1.)/
      data capom/7.292115e-05/       ! siderial angular velocity
      data pai/3.14159 62535 9/

      s(l,m)=sqrt(float((l+m)*(l-m))/float((2*l+1)*(2*l-1)))


      w0r=cabs(w0)

      w0n2=(w0r/wn)**2

c     do the radial integrations

      
      call intmdcpl_new(isw2,fus,minll,maxll)


      if(jcom1.eq.1)jcom1=3
      if(jcom2.eq.1)jcom2=3

      do 50 m1=-lord1,lord1
      do 50 m2=-lord2,lord2
      zr0(m1,m2) = 0.d0
      zr1(m1,m2) = 0.d0
!      zr(m1,m2)=0.d0
      z0(m1,m2) = 0.d0
      z1(m1,m2) = 0.d0
      z2(m1,m2) = 0.d0
!   50 z(m1,m2)=0.d0
 50   continue

      fli=float(lord1)
      flii=float(lord2)
      fl3i=fli*(fli+1.)


cc  coupling due to rotation
cc     z will be 0 unless:mord1=mord2 and
cc                         jcom1.eq.jcom2 and lord1=lord2 , or
cc                         jcom1.ne.jcom2 and lord1=lord2(+,-)1

      lmin=min(lord1,lord2)

      if(lord1.ne.lord2.or.jcom1.ne.jcom2)goto 110

c     see equations (A17) to (A20) in Woodhouse 1980 for the
c     required formulae - note that rot is the returned by 
c     the routine intmdcpl, and is equal to the integral
c     of rho_{0} C^{(+)}r^{2}. 

      do 100 m1=-lmin,lmin
      fm=float(m1)
      if(jcom1.eq.3) then
!         zr(m1,m1)=dcmplx((-(2./3.)*capom*fl3i+2.*w0*fm)*capom*rot)      
         zr0(m1,m1) = dcmplx((-(2./3.)*capom*fl3i)*capom*rot)
         zr1(m1,m1) = dcmplx(2.*fm*capom*rot)  
      end if
      if(nord1.eq.nord2)then
       if(jcom1.eq.2)  then
!          zr(m1,m1)=dcmplx(2.*fm*capom*w0/fl3i)
          zr1(m1,m1) = dcmplx(2.*fm*capom/fl3i)
       end if
       if(jcom1.eq.3) then
!          zr(m1,m1)=zr(m1,m1)+dcmplx((2./3.)*capom**2)
          zr0(m1,m1)=zr0(m1,m1)+dcmplx((2./3.)*capom**2)
       end if
      endif
  100 continue

  110 if(jcom1.eq.jcom2)goto 199
      if(lord1.ne.(lord2-1).and.lord1.ne.(lord2+1))goto 199
      lmax=max(lord1,lord2)
      do m1=-lmin,lmin
!         zr(m1,m1)=dcmplx(iunt*s(lmax,m1)*2.*capom*rot*w0)
         zr1(m1,m1)=dcmplx(iunt*s(lmax,m1)*2.*capom*rot)
      end do
  199 continue


cc  coupling due to asphericity

      fac=sqrt((2.*fli+1.)*(2.*flii+1.)/pai)/2.
      ml2=2*ML+1

      do 200 ll=minll,maxll

      call wig2(lord1,ll,lord2,sj(-lord1,-lord2),ml2)

!      tt(ll,0)=asp(ll,1)
      tt0(ll,0)=asp0(ll,1)
      tt2(ll,0)=asp2(ll,1)
ca      if(ll.eq.0)goto 200
      do mm=-ll,-1
         ma=-2*mm
!         tt(ll,mm)=(-1)**mm*0.5*(asp(ll,ma)+asp(ll,ma+1)*iunt)
         tt0(ll,mm)=(-1)**mm*0.5*(asp0(ll,ma)+asp0(ll,ma+1)*iunt)
         tt2(ll,mm)=(-1)**mm*0.5*(asp2(ll,ma)+asp2(ll,ma+1)*iunt)
      end do
      do mm=1,ll
         ma=2*mm
!         tt(ll,mm)=0.5*(asp(ll,ma)-asp(ll,ma+1)*iunt)
         tt0(ll,mm)=0.5*(asp0(ll,ma)-asp0(ll,ma+1)*iunt)
         tt2(ll,mm)=0.5*(asp2(ll,ma)-asp2(ll,ma+1)*iunt)
      end do

      do 210 m2=-lord2,lord2
      do 210 m1=-lord1,lord1
cc       output of wig2 will be 0,unless: mm=m1-m2
      if(iabs(m1-m2).gt.ll)goto 210
!      t=tt(ll,m1-m2)*dcmplx(sqrt(2.*float(ll)+1.)*fac*sj(m1,m2))
      t0=tt0(ll,m1-m2)*dcmplx(sqrt(2.*float(ll)+1.)*fac*sj(m1,m2))
      t2=tt2(ll,m1-m2)*dcmplx(sqrt(2.*float(ll)+1.)*fac*sj(m1,m2))
      if(jcom1.eq.jcom2) then
!         z(m1,m2)=z(m1,m2)+dcmplx(t)
         z0(m1,m2)=z0(m1,m2)+dcmplx(t0)
         z2(m1,m2)=z2(m1,m2)+dcmplx(t2)
      end if
      if(jcom1.ne.jcom2) then
!         z(m1,m2)=z(m1,m2)+dcmplx(t*iunt)
         z0(m1,m2)=z0(m1,m2)+dcmplx(t0*iunt)
         z2(m1,m2)=z2(m1,m2)+dcmplx(t2*iunt)
      end if
  210 continue
      

  200 continue

      do 300 m2=-lord2,lord2
      do 300 m1=-lord1,lord1
         if(isw1.eq.1) then
!            z(m1,m2)=z(m1,m2)+zr(m1,m2)
            z0(m1,m2)=z0(m1,m2)+zr0(m1,m2)
            z1(m1,m2)=zr1(m1,m2)
         end if
c  300 write(7,"(2e12.4)")z(m1,m2)
  300 continue
  400  return
      end
