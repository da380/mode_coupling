      subroutine intmdcpl_new(isw,fus,mins,maxs)
      parameter(MLL=20)
      parameter(MM=2*MLL+1)

!      common/int/rot,htg
      common/int/rot

      common/int0/htg0,htg2

      common/grav/grv,qgrv

      common/hetswitch/ifhet

      real lcon,ncon,llp,nnp,ll,nn
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/nond/rn,wn,vn,gn,rhobar
      common/eig1/nord1,jcom1,lord1,wcom1,qbar1,cgp1,avert1,ahor1,phis1
     1          ,eif1(222,6)
      common/eig2/nord2,jcom2,lord2,wcom2,qbar2,cgp2,avert2,ahor2,phis2
     1          ,eif2(222,6)

c  common block for harmonic coefficients of the model
c  perturbations -- perturbations in A, C, L, N, F, Rho, Phi0, Phi0'
c  In this program Rho, Phi0, Phi0' are not perturbed -- upgrade in the future.

      common/vioem/da(0:MLL,MM,222,5),dc(0:MLL,MM,222,5),df(0:MLL,MM,222,
     1          5),dl(0:MLL,MM,222,5),dn(0:MLL,MM,222,5),dr(0:MLL,MM,222,
     2          5),dp(0:MLL,MM,222,5),dpp(0:MLL,MM,222,5),
     3             rhoin(222,5),grin(222,5),elin(222,5),etin(222,5)
      common/hetmdl/iw(5),crt(81)

      dimension u1(222),up1(222),v1(222),vp1(222),p1(222),pp1(222),
     1          u2(222),up2(222),v2(222),vp2(222),p2(222),pp2(222)
      dimension dta(0:MLL,MM),dtc(0:MLL,MM),dtf(0:MLL,MM),dtl(0:MLL,MM),
     1          dtn(0:MLL,MM),dtr(0:MLL,MM),dtp(0:MLL,MM),dtpp(0:MLL,MM),
     2          h(0:MLL,MM)
      dimension q1(3),qp1(3),q2(3),qp2(3)
      dimension xi(5),wt(5),a1(3),b1(3),a2(3),b2(3)

      equivalence (eif1(1,1),u1(1)),(eif1(1,2),up1(1)),(eif1(1,3),v1(1)),
     1            (eif1(1,4),vp1(1)),(eif1(1,5),p1(1)),(eif1(1,6),pp1(1)),
     2            (eif2(1,1),u2(1)),(eif2(1,2),up2(1)),(eif2(1,3),v2(1)),
     3            (eif2(1,4),vp2(1)),(eif2(1,5),p2(1)),(eif2(1,6),pp2(1))
      equivalence (q1(1),uri),(q1(2),vri),(q1(3),pri),
     1            (qp1(1),upri),(qp1(2),vpri),(qp1(3),ppri),
     2            (q2(1),urii),(q2(2),vrii),(q2(3),prii),
     3            (qp2(1),uprii),(qp2(2),vprii),(qp2(3),pprii)

      double precision temp(2,7),acum(2,7),add,rot,b(5),
     1                 asp(0:MLL,MM,7),ast(0:MLL,MM,7),htg(0:MLL,MM),
     1                 asp0(0:MLL,MM,7),ast0(0:MLL,MM,7),htg0(0:MLL,MM), 
     1                 asp2(0:MLL,MM,7),ast2(0:MLL,MM,7),htg2(0:MLL,MM),
     1                 temp0(2,7),temp2(2,7),add0,add2,acum0(2,7),
     1                 acum2(2,7)
      real grv(222,0:MLL),qgrv(3,222,0:MLL),gint(0:MLL)

c Gauss-Legendre absicae and weights from Abromowitz & Stegun

      data xi/-0.90617 98459 38664,-0.53846 93101 05683
     1       , 0.00000 00000 00000, 0.53846 93101 05683
     2       , 0.90617 98459 38664/
      data wt/ .23692 68850 56189, 0.47862 86704 99366
     1       , .56888 88888 88889, 0.47862 86704 99366
     2       , .23692 68850 56189/


      data n670/180/,n220/202/,nmidc/216/
      data capom/7.292115e-05/
      data coe/1.05688 76266 5/
c             coe=(sqpai/sqrt(5.))*4/3.
      data nint/2/
!      data drn,drnsl,drmidc,drmoho
!     1   ,droc,druc,drlc
!     2   ,dvpoc,dvpuc,dvplc
!     3   ,dvsoc,dvsuc,dvslc/
!     4     3.5028e-02, 2.1017e-01,-3.7886e-01,-7.8982e-01
!     5   ,-3.5028e-04, 0.        , 0.
!     6   , 0.        ,-2.1017e-02, 3.5028e-03
!     7   , 0.        ,-5.2542e-03, 1.7514e-03/

!    set crustal model to zero
      data drn,drnsl,drmidc,drmoho
     1   ,droc,druc,drlc
     2   ,dvpoc,dvpuc,dvplc
     3   ,dvsoc,dvsuc,dvslc/
     4     0., 0., 0., 0.
     5   , 0., 0., 0.
     6   , 0.        , 0. , 0.
     7   , 0.        , 0.,  0./







      call gravk(mins,maxs)


      if(jcom1.eq.0.or.jcom2.eq.0) return

      if(jcom1.eq.1)jcom1=3
      if(jcom2.eq.1)jcom2=3
      jcom=3
      if(jcom1.eq.2.and.jcom2.eq.2)jcom=2

      fli=float(lord1)
      flii=float(lord2)
      fl3i=fli*(fli+1.)
      fl3ii=flii*(flii+1.)

      capon2=(capom/wn)**2
      delv=1000./vn
      deldis=1000./rn
      wcom=sqrt(wcom1*wcom2)
      omn2=wcom1*wcom2/wn**2


      j11=1
      j12=3
      if(jcom1.eq.2) then
        j11=2
        j12=2
      endif

      j21=1
      j22=3
      if(jcom2.eq.2) then
        j21=2
        j22=2
      endif

      do 101 i=1,nint
      do 101 j=1,7
      acum0(i,j)=0.    
      acum2(i,j)=0.    
 101  continue

      do 102 i=0,MLL
      do 102 j=1,2*i+1
      do 102 k=1,7
      asp0(i,j,k)=0.
      asp2(i,j,k)=0.
 102  continue


      nstart=1
      if(jcom1.eq.2.or.jcom2.eq.2) nstart=noc



      do 1000 iq=nstart,n

      do 1001 i=1,nint
      do 1001 j=1,7
      temp0(i,j)=0.
      temp2(i,j)=0.
 1001 continue

      do 1002 i=mins,maxs
      do 1002 j=1,2*i+1
      do 1002 k=1,7
      ast0(i,j,k)=0.
      ast2(i,j,k)=0.
 1002  continue

      if(iq.eq.n) goto 1200

      call corfac(iq,wcom,jcom,xac,xf,xln)

      iq1=iq+1
      r1=r(iq)
      r2=r(iq1)
      hn=r2-r1
      hnh=hn*.5
      if(hn.lt.1.e-4) goto 1200
      hr=1./hn
      hsq=hr*hr
      hcu=hr*hsq

      do 1003 i=j11,j12
      i1=2*i-1
      i2=i1+1
      a1(i)=(eif1(iq,i2)+eif1(iq1,i2))*hsq
     1     +2.*(eif1(iq,i1)-eif1(iq1,i1))*hcu
 1003 b1(i)=-(2.*eif1(iq,i2)+eif1(iq1,i2))*hr
     1     -3.*(eif1(iq,i1)-eif1(iq1,i1))*hsq

      do 1004 i=j21,j22
      i1=2*i-1
      i2=i1+1
      a2(i)=(eif2(iq,i2)+eif2(iq1,i2))*hsq
     1     +2.*(eif2(iq,i1)-eif2(iq1,i1))*hcu
 1004 b2(i)=-(2.*eif2(iq,i2)+eif2(iq1,i2))*hr
     1     -3.*(eif2(iq,i1)-eif2(iq1,i1))*hsq


      do 1100 il=1,5

      t=.5*hn*(xi(il)+1.)
      rr=r1+t
      rr2=rr*rr

      do 1101 i=j11,j12
      i1=2*i-1
      i2=i1+1
      q1(i)=eif1(iq,i1)+t*(eif1(iq,i2)+t*(b1(i)+t*a1(i)))
 1101 qp1(i)=(eif1(iq,i2)+t*(2.*b1(i)+t*3.*a1(i)))*rr

      do 1102 i=j21,j22
      i1=2*i-1
      i2=i1+1
      q2(i)=eif2(iq,i1)+t*(eif2(iq,i2)+t*(b2(i)+t*a2(i)))
 1102 qp2(i)=(eif2(iq,i2)+t*(2.*b2(i)+t*3.*a2(i)))*rr


      el=elin(iq,il)
      et=etin(iq,il)
      gr=grin(iq,il)

      if(jcom1.eq.2.or.jcom2.eq.2) goto 1103
      aa=xac*(acon(iq)+t*(qacon(1,iq)+t*(qacon(2,iq)+t*qacon(3,iq))))
      cc=xac*(ccon(iq)+t*(qccon(1,iq)+t*(qccon(2,iq)+t*qccon(3,iq))))
      ff=xf*(fcon(iq)+t*(qfcon(1,iq)+t*(qfcon(2,iq)+t*qfcon(3,iq))))
      aap=xac*(qacon(1,iq)+t*(2.*qacon(2,iq)+t*3.*qacon(3,iq)))
      ccp=xac*(qccon(1,iq)+t*(2.*qccon(2,iq)+t*3.*qccon(3,iq)))
      ffp=xf*(qfcon(1,iq)+t*(2.*qfcon(2,iq)+t*3.*qfcon(3,iq)))
 1103 ll=xln*(lcon(iq)+t*(qlcon(1,iq)+t*(qlcon(2,iq)+t*qlcon(3,iq))))
      nn=xln*(ncon(iq)+t*(qncon(1,iq)+t*(qncon(2,iq)+t*qncon(3,iq))))
      rrho=rho(iq)+t*(qrho(1,iq)+t*(qrho(2,iq)+t*qrho(3,iq)))
      llp=xln*(qlcon(1,iq)+t*(2.*qlcon(2,iq)+t*3.*qlcon(3,iq)))
      nnp=xln*(qncon(1,iq)+t*(2.*qncon(2,iq)+t*3.*qncon(3,iq)))
      rhop=qrho(1,iq)+t*(2.*qrho(2,iq)+t*3.*qrho(3,iq))

      do 1104 i=mins,maxs
      do 1104 j=1,2*i+1
      dta(i,j)=da(i,j,iq,il)
      dtc(i,j)=dc(i,j,iq,il)
      dtf(i,j)=df(i,j,iq,il)
      dtl(i,j)=dl(i,j,iq,il)
      dtn(i,j)=dn(i,j,iq,il)
      dtr(i,j)=dr(i,j,iq,il)
      dtp(i,j)=dp(i,j,iq,il)
 1104 dtpp(i,j)=dpp(i,j,iq,il)


 6000 if(isw.eq.0)goto 5001
      dta(2,1)=dta(2,1)+coe*rr*el*aap
      dtc(2,1)=dtc(2,1)+coe*rr*el*ccp
      dtf(2,1)=dtf(2,1)+coe*rr*el*ffp
      dtl(2,1)=dtl(2,1)+coe*rr*el*llp
      dtn(2,1)=dtn(2,1)+coe*rr*el*nnp
      dtr(2,1)=dtr(2,1)+coe*rr*el*rhop      
      dtp(2,1)=dtp(2,1)+coe*0.5*capon2*rr2 
      dtpp(2,1)=dtpp(2,1)+coe*capon2*rr


 5001 continue

      do i=mins,maxs
         gint(i)=grv(iq,i)+t*(qgrv(1,iq,i)+t*(qgrv(2,iq,i)+t*qgrv(3,iq,i)))
      enddo

      xxi=vpri-vri
      xxii=vprii-vrii
      vr2=vri*vrii

      if(jcom1.eq.2) goto 1105
      fi=2.*uri-fl3i*vri
      xxi=xxi+uri
 1105 if(jcom2.eq.2) goto 1110
      fii=2.*urii-fl3ii*vrii
      xxii=xxii+urii

 1110 xx2=xxi*xxii
      if(jcom.eq.2)goto 1106
      uivii=uri*vrii
      uiivi=urii*vri
      fivii=fi*vrii
      fiivi=fii*vri
      upivii=upri*vrii
      upiivi=uprii*vri
      vpiuii=vpri*urii
      vpiiui=vprii*uri



 1106 if(jcom1.ne.jcom2)goto 1120

      ccc=vr2
      rkl1=xx2
      rkn2=vr2
      rkr10 = 0.
      rkr12 = -vr2*rr2
c ---- I think the following two values have to be set to zero here -----
      rkp1=0.
      rkpp1=0.

      if(jcom1.eq.2) goto 1150

      ur2=uri*urii
      fiuii=fi*urii
      fiiui=fii*uri

      ccc=vr2+uivii+uiivi
      rka0=fi*fii
      rkc0=upri*uprii
      rkf0=uprii*fi+upri*fii
      rkn0=-rka0
      rkr00 = ((8.*rrho)*ur2*rr+ppri*urii+pprii*uri
     1     -.5*gr*(4.*ur2+fiuii+fiiui))*rr 
      rkr02 = (-ur2*rr)*rr 
      rkr10=rkr10+(pri*vrii+prii*vri+.5*gr*(uivii+uiivi))*rr
      rkp0a=rrho*ur2
      rkp1a=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
      rkp1b=.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
      rkpp0=-rrho*(fiuii+fiiui)*rr
      rkpp1a=.5*rrho*rr*uivii
      rkpp1b=.5*rrho*rr*uiivi

      goto 1150

 1120 rkl1=xx2
      rkn2=vr2
      rkr10 = 0.
      rkr12 = -vr2*rr2
      if(jcom1.eq.2)then
        ccc=.5*((fl3i+fl3ii-2.)*vr2-(fl3i-fl3ii+2.)*uiivi)
        rkr10=rkr10+(prii+.5*gr*urii)*vri*rr
        rkp1=-.5*rrho*(vpiuii+uiivi-upiivi-2.*fiivi)
        rkpp1=-.5*rrho*uiivi*rr
      else
        ccc=-.5*((fl3ii+fl3i-2.)*vr2-(fl3ii-fl3i+2.)*uivii)
        rkl1=-rkl1
        rkn2=-rkn2
        rkr10=-rkr10-(pri+.5*gr*uri)*vrii*rr
        rkr12=-rkr12
        rkp1=.5*rrho*(vpiiui+uivii-upivii-2.*fivii)
        rkpp1=.5*rrho*uivii*rr
      endif

 1150 continue



c
c   normalization integral
c
      if(jcom1.ne.jcom2) goto  1160
      if(lord1.ne.lord2) goto  1160
      t1=vr2*fl3i
      if(jcom1.ne.2) t1=t1+ur2
      add0=wt(il)*rr2*t1*rrho*hnh
      temp0(1,1)=temp0(1,1)+add0
 1160 continue

c
c   rotational splitting
c
      add0=wt(il)*rr2*ccc*rrho*hnh
      temp0(1,2)=temp0(1,2)+add0



c
c   aspherical perturbation (including ellipticity)
c

      do 1170 i=mins,maxs
      do 1170 j=1,2*i+1
      add0 = wt(il)*(rkl1*dtl(i,j)+rkr10*dtr(i,j)+rkp1*dtp(i,j)+rkpp1*dtpp(
     1    i,j))*hnh
      add2 = wt(il)*(rkr12*dtr(i,j))*hnh
      ast0(i,j,2)=ast0(i,j,2)+add0
      ast2(i,j,2)=ast2(i,j,2)+add2
      add0=wt(il)*rkn2*dtn(i,j)*hnh
      ast0(i,j,3)=ast0(i,j,3)+add0
      add0=wt(il)*gint(i)*dtr(i,j)*rr2*hnh
      ast0(i,j,7)=ast0(i,j,7)+add0

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1170
      add0=wt(il)*(rka0*dta(i,j)+rkc0*dtc(i,j)+rkf0*dtf(i,j)+rkn0*dtn(i,j
     1         )+rkr00*dtr(i,j)+rkpp0*dtpp(i,j))*hnh
      add2=wt(il)*(rkr02*dtr(i,j))*hnh
      ast0(i,j,1)=ast0(i,j,1)+add0
      ast2(i,j,1)=ast2(i,j,1)+add2
      add0=wt(il)*rkp0a*dtp(i,j)*hnh
      ast0(i,j,4)=ast0(i,j,4)+add0
      add0=wt(il)*(rkp1a*dtp(i,j)+rkpp1a*dtpp(i,j))*hnh
      ast0(i,j,5)=ast0(i,j,5)+add0
      add0=wt(il)*(rkp1b*dtp(i,j)+rkpp1b*dtpp(i,j))*hnh
      ast0(i,j,6)=ast0(i,j,6)+add0
 1170 continue

c
c    crustal correction
c
 5003 if(iq.lt.moho) goto 1184
      vs=sqrt(ll/rrho)
      vvp=sqrt(aa/rrho)
      if(iq.gt.nmidc) goto 1181
      ddrho=drlc
      ddvs2=2.*vs*dvslc*delv
      ddvp2=2.*vvp*dvplc*delv
      goto 1183
 1181 if(iq.gt.nsl) goto 1182
      ddrho=druc
      ddvs2=2.*vs*dvsuc*delv
      ddvp2=2.*vvp*dvpuc*delv
      goto 1183
 1182 ddrho=droc
      ddvs2=2.*vs*dvsoc*delv
      ddvp2=2.*vvp*dvpoc*delv
 1183  continue
      add0=wt(il)*hnh*(rrho*ddvs2*rkl1+(ll*rkl1+rrho*rkr10)*ddrho/rrho)
      add2=wt(il)*hnh*(rrho*rkr12*ddrho/rrho)
      temp0(2,2)=temp0(2,2)+add0
      temp2(2,2)=temp2(2,2)+add2
      add0=wt(il)*hnh*(rrho*ddvs2+nn*ddrho/rrho)*rkn2
      temp0(2,3)=temp0(2,3)+add0
      if(jcom1.eq.2.or.jcom2.eq.2)goto 1184
      etan=ff/(aa-2.*ll)
      add0=wt(il)*hnh*(rrho*(ddvp2*(rka0+rkc0+etan*rkf0)+ddvs2*(rkn0-2.*
     1 etan*rkf0))+(aa*rka0+cc*rkc0+ff*rkf0+nn*rkn0+rrho*rkr00)*ddrho/rrho)
      add2=wt(il)*hnh*((rrho*rkr02)*ddrho/rrho)
      temp0(2,1)=temp0(2,1)+add0
      temp2(2,1)=temp2(2,1)+add2
 1184 continue

 1100 continue

      goto  1300


c
c  discontinuity contribution
c

 1200 rr=r(iq)
      rr2=rr*rr

      do 1299 idis=1,3,2

      iqt=iq+(idis-1)/2
      if(iqt.gt.n) goto 1299

      call corfac(iqt,wcom,jcom,xac,xf,xln)

      do 1201 i=j11,j12
      i1=2*i-1
      i2=i1+1
      q1(i)=eif1(iqt,i1)
 1201 qp1(i)=eif1(iqt,i2)*rr

      do 1202 i=j21,j22
      i1=2*i-1
      i2=i1+1
      q2(i)=eif2(iqt,i1)
 1202 qp2(i)=eif2(iqt,i2)*rr

      el=ell(iqt)

      if(jcom.eq.2) goto 1203
      if(jcom1.ne.2.and.jcom2.ne.2)  aa=xac*acon(iqt)
      gr=g(iqt)
      cc=xac*ccon(iqt)
      ff=xf*fcon(iqt)
 1203 ll=xln*lcon(iqt)
      nn=xln*ncon(iqt)
      rrho=rho(iqt)

      do i=mins,maxs
         gint(i)=grv(iqt,i)
      end do

      vr2=vri*vrii
      xxi=vpri-vri
      xxii=vprii-vrii

      if(jcom1.eq.2) goto 4405
      xxi=xxi+uri
      fi=2.*uri-fl3i*vri
 4405 if(jcom2.eq.2) goto 4440
      xxii=xxii+urii
      fii=2.*urii-fl3ii*vrii

4440  xx2=xxi*xxii

      if(jcom1.ne.jcom2)goto 1220

      rkl1=xx2-vpri*xxii-vprii*xxi
      rkn2=vr2
      rkr10 = 0.
      rkr12 = -vr2*rr2

      if(jcom1.eq.2)goto 1230

      fi=2.*uri-fl3i*vri
      fii=2.*urii-fl3ii*vrii
      ur2=uri*urii

      rka0=fi*fii
      rkc0=-upri*uprii
      rkc11=upri*vrii
      rkc12=uprii*vri
      rkf11=fi*vrii
      rkf12=fii*vri
      rkn0=-rka0
      rkr00 = ((8.*rrho)*ur2*rr+(ppri*urii+pprii*uri)
     1     -.5*gr*(4.*ur2+fi*urii+fii*uri))*rr 
      rkr02 = (-ur2*rr)*rr 
      rkr10 = rkr10 + ((pri*vrii+prii*vri)+.5*gr*(uri*vrii+urii*vri))*rr
      goto 1230

 1220 rkr10 = 0.
      rkr12 = -vr2*rr2
      if(jcom1.eq.2)then
        rkc1=-uprii*vri
        rkf1=-fii*vri
        rkl1=xx2-vprii*xxi-vpri*xxii
        rkn2=vr2
        rkr10=rkr10+(prii+.5*gr*urii)*vri*rr
      else
        rkc1=upri*vrii
        rkf1=fi*vrii
        rkl1=-xx2+vpri*xxii+xxi*vprii
        rkn2=-vr2
        rkr10=-rkr10-(pri+.5*gr*uri)*vrii*rr
        rkr12=-rkr12
      endif


c
c    aspherical perturbation(only ellipticity)
c


 1230 if(isw.eq.0)goto 1240

      do i=mins,maxs
         do j=1,2*i+1
            h(i,j)=0.
      enddo
      enddo

      h(2,1)=-coe*rr*el

      do 1239 i=mins,maxs
      do 1239 j=1,2*i+1

      ttt=float(idis-2)*h(i,j)

      add0=-ttt*(ll*rkl1+rrho*rkr10)
      add2=-ttt*(rrho*rkr12)
      ast0(i,j,2)=ast0(i,j,2)+add0
      ast2(i,j,2)=ast2(i,j,2)+add2
      add0=-ttt*nn*rkn2
      ast0(i,j,3)=ast0(i,j,3)+add0
      add0=-ttt*gint(i)*rrho*rr2
      ast0(i,j,7)=ast0(i,j,7)+add0

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1231
      add0=-ttt*(aa*rka0+cc*rkc0+nn*rkn0+rrho*rkr00)
      add2=-ttt*(rrho*rkr02)
      ast0(i,j,1)=ast0(i,j,1)+add0
      ast2(i,j,1)=ast2(i,j,1)+add2
      add0 = -ttt*(cc*rkc11+ff*rkf11)
      ast0(i,j,5)=ast0(i,j,5)+add0
      add0=-ttt*(cc*rkc12+ff*rkf12)
      ast0(i,j,6)=ast0(i,j,6)+add0

 1231 if(jcom1.eq.jcom2)goto 1239
      add0=-ttt*(cc*rkc1+ff*rkf1)
      ast0(i,j,2)=ast0(i,j,2)+add0

 1239 continue

c
c  crustal thickness
c
 1240 if(iq.ne.nmidc.and.iq.ne.moho.and.iq.ne.nsl.and.iq.ne.n) goto 1249
      if(iq.eq.n) delh=drn*deldis
      if(iq.eq.nsl) delh=drnsl*deldis
      if(iq.eq.nmidc) delh=drmidc*deldis
      if(iq.eq.moho) delh=drmoho*deldis
      ttt=float(idis-2)*delh

      add0=-ttt*(ll*rkl1+rrho*rkr10)
      add2 = -ttt*(rrho*rkr12)
      temp0(2,2)=temp0(2,2)+add0
      temp2(2,2)=temp2(2,2)+add2

      add0=-ttt*nn*rkn2
      temp0(2,3)=temp0(2,3)+add0

      if(jcom1.eq.2.or.jcom2.eq.2)goto 1241
      add0=-ttt*(aa*rka0+cc*rkc0+nn*rkn0+rrho*rkr00)
      add2=-ttt*(rrho*rkr02)
      temp0(2,1)=temp0(2,1)+add0
      temp2(2,1)=temp2(2,1)+add2
      add0=-ttt*(cc*rkc11+ff*rkf11)
      temp0(2,5)=temp0(2,5)+add0
      add0=-ttt*(cc*rkc12+ff*rkf12)
      temp0(2,6)=temp0(2,6)+add0

 1241 if(jcom1.eq.jcom2)goto 1249
      add0=-ttt*(cc*rkc1+ff*rkf1)
      temp0(2,2)=temp0(2,2)+add0
 1249 continue

 1299 continue

 1300 do 1301 i=1,nint
      do 1301 j=1,7
      acum0(i,j)=acum0(i,j)+temp0(i,j)
      acum2(i,j)=acum2(i,j)+temp2(i,j)
 1301 continue

      do 1302 i=mins,maxs
      do 1302 j=1,2*i+1
      do 1302 k=1,7
      asp0(i,j,k)=asp0(i,j,k)+ast0(i,j,k)
      asp2(i,j,k)=asp2(i,j,k)+ast2(i,j,k)
 1302 continue


 1000 continue

      anorm=acum0(1,1)*omn2
      rot=acum0(1,2)*omn2

      ind=mins**2
      do 2000 i=mins,maxs

      call bcoff1(lord1,i,lord2,b)

      if(i.le.8) then
      do 2002 j=1,2*i+1
      ind=ind+1
      htg0(i,j)=crt(ind)*fus
      htg2(i,j)=crt(ind)*fus
      do 2002 k=1,6
      asp0(i,j,k)=asp0(i,j,k)+acum0(2,k)*htg0(i,j)
      asp2(i,j,k)=asp2(i,j,k)+acum2(2,k)*htg2(i,j)
 2002 continue
      endif

      if(jcom1.ne.jcom2)goto 2020
      fl3=float(i*(i+1))
      do 2010 j=1,2*i+1
      htg0(i,j)=((asp0(i,j,1)+asp0(i,j,4)*fl3+.5*(asp0(i,j,2)*(fl3i+fl3ii-
     1   fl3)+asp0(i,j,5)*(fl3ii+fl3-fl3i)+asp0(i,j,6)*(fl3i+fl3-fl3ii))
     2            )*b(1)+asp0(i,j,3)*b(3)+asp0(i,j,7)*b(1))*omn2*wn**2
      htg2(i,j)=((asp2(i,j,1)+asp2(i,j,4)*fl3+.5*(asp2(i,j,2)*(fl3i+fl3ii-
     1   fl3)+asp2(i,j,5)*(fl3ii+fl3-fl3i)+asp2(i,j,6)*(fl3i+fl3-fl3ii))
     2            )*b(1)+asp2(i,j,3)*b(3)+asp2(i,j,7)*b(1))*omn2*wn**2
 2010 continue
      goto 2000

 2020 do 2021 j=1,2*i+1
      htg0(i,j)=(asp0(i,j,2)*b(4)+asp0(i,j,3)*b(5)+asp0(i,j,7)*b(4))*omn2*wn**2
      htg2(i,j)=(asp2(i,j,2)*b(4)+asp2(i,j,3)*b(5)+asp2(i,j,7)*b(4))*omn2*wn**2
 2021 continue

 2000 continue


      return
      end
