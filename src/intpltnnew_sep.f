c-----------------------------------------------
c This routine evaluates the harmonic coefficients
c of delta(A), delta(C) etc. using the input heterogeneous
c models. It will need to be changed to accomodate different
c model parameterisations etc.
c

      subroutine intpltnnew_sep(hetmodel_vs,hetmodel_vp,hetmodel_rho)

      parameter(MLL=20)
      parameter(MM=2*MLL+1)
      parameter (MXPARM=30)
      parameter (MXLENY=(MLL+1)**2)
      parameter (MXMDLL=MXLENY*MXPARM)


      real lcon,ncon,ll,nn
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)
      common/nond/rn,wn,vn,gn,rhobar
      common/vioem/da(0:MLL,MM,222,5),dc(0:MLL,MM,222,5),df(0:MLL,MM,222,
     1          5),dl(0:MLL,MM,222,5),dn(0:MLL,MM,222,5),dr(0:MLL,MM,222,
     2          5),dp(0:MLL,MM,222,5),dpp(0:MLL,MM,222,5),
     3             rhoin(222,5),grin(222,5),elin(222,5),etin(222,5)
      common/lmmdl/iz(4),p(588),b1(588)   ! L02.56
      common/hetmdl/iw(5),s(588),b2(588)  ! M84C




      character*80 file

      dimension dvs2(0:MLL,MM,4),dvp(0:MLL,MM,5),xi(5)

      data xi/-0.90617 98459 38664,-0.53846 93101 05683
     1       , 0.00000 00000 00000, 0.53846 93101 05683
     2       , 0.90617 98459 38664/
      data thrd,fot/.33333 33333 33333,1.33333 33333 33333/
      data n670/180/,n220/202/,nmidc/216/


      character*(*) hetmodel_vs
      character*(*) hetmodel_vp
      character*(*) hetmodel_rho
      character*2 ftype
      common/nmdl/lmaxn,nstr,nlp,anmdl(MXMDLL),bnmdl(MXMDLL)
      common/nmdl_vs/anmdl_vs(MXMDLL),bnmdl_vs(MXMDLL)
      common/nmdl_vp/anmdl_vp(MXMDLL),bnmdl_vp(MXMDLL)
      common/nmdl_rho/anmdl_rho(MXMDLL),bnmdl_rho(MXMDLL)
      dimension mask_vs(MXPARM),mask1_vs(MXPARM),
     1     lask_vs(0:MLL),lask1_vs(0:MLL)
      dimension mask_vp(MXPARM),mask1_vp(MXPARM),
     1     lask_vp(0:MLL),lask1_vp(0:MLL)
      dimension mask_rho(MXPARM),mask1_rho(MXPARM),
     1     lask_rho(0:MLL),lask1_rho(0:MLL)
      dimension dvsr(0:MLL,MM,MXPARM)
      dimension dvpr(0:MLL,MM,MXPARM)
      dimension drhor(0:MLL,MM,MXPARM)

c Legendre polynomials

      p1(x)=sqrt(1.5)*x
      p2(x)=sqrt(5./8.)*(3.*x*x-1.)
      p3(x)=sqrt(7./8.)*(5.*x**3-3.*x)
      p4(x)=sqrt(9./128.)*(35.*x**4-3.*x**2+3.)
      p0=sqrt(.5)



      call rdmdl(3)

      open(27,file=hetmodel_vs)
      rewind(27)
      read(27,'(a)') ftype
      close(27)

      if(ftype.eq.'LG') then

         stop 'can''t read in this format!!'


      else

ca-------------------------------------------------
ca Read in mantle model in splines (John's format)
ca-------------------------------------------------


!      call heread(1,hetmodel,bnmdl,nstr,lmaxn,mask,lask,0)

!     read in vs model                                         
      call heread(1,hetmodel_vs,bnmdl_vs,nstr,lmaxn,mask_vs,lask_vs,0)

!     read in vp model                                         
      call heread(1,hetmodel_vp,bnmdl_vp,nstr,lmaxn,mask_vp,lask_vp,0)

!     read in rho model                                         
      call heread(1,hetmodel_rho,bnmdl_rho,nstr,lmaxn,mask_rho, 
     1     lask_rho,0)

      do i=1,nstr
         mask1_vs(i)=1
         mask1_vp(i)=1
         mask1_rho(i)=1
      enddo
      do i=0,lmaxn
         lask1_vs(i)=1
         lask1_vp(i)=1
         lask1_rho(i)=1
      enddo

      call mskmdl(bnmdl_vs,nstr,lmaxn,mask_vs,lask_vs
     1       ,anmdl_vs,nstr,lmaxn,mask1_vs,lask1_vs) 
      call mskmdl(bnmdl_vp,nstr,lmaxn,mask_vp,lask_vp
     1       ,anmdl_vp,nstr,lmaxn,mask1_vp,lask1_vp) 
      call mskmdl(bnmdl_rho,nstr,lmaxn,mask_rho,lask_rho
     1       ,anmdl_rho,nstr,lmaxn,mask1_rho,lask1_rho) 
      leny=(lmaxn+1)**2

c     The 'sph' model format is defined in terms of splines
c     in depth. The argument of these spline functions is defined
c     to be -1 at CMB and +1 at moho (in PREM).
c     To set up the parameters needed to evaluate the splines:

      call splhsetup()

c     splines can now be evaluated using
c     value =splh(ispl, -1+2.0*(r-rcmb)/(rmoho-rcmb)  )
c     with ispl = 0,1,2,3, .... 20.

c     The value represented in 'sph' format is delta(vs)/vs
c

      fac=1000./vn

      do k=4,nstr ! in sph format dvs/vs is defined by parameters 4 to 24
        ind=0
        do i=0,MLL
          do j=1,2*i+1
            ind=ind+1
            dvsr(i,j,k) = anmdl_vs(leny*(k-1)+ind)  ! actually dvs/vs
            dvpr(i,j,k) = anmdl_vp(leny*(k-1)+ind)  ! actually dvp/vp
            drhor(i,j,k) = anmdl_rho(leny*(k-1)+ind)  ! actually drho/rho
         enddo
        enddo
      enddo

      endif
      print *, r(1),r(222),rho(1)
      do 1000 iq=1,n  ! do loop over PREM levels iq=1 (centre) to iq=n (surface)
      iq1=iq+1
      r1=r(iq)
      r2=r(iq1)
      hn=r2-r1
      hnh=hn*.5
      if(hn.lt.1.e-4)goto 1000
      hr=1./hn
      hsq=hr*hr
      hcu=hr*hsq

      elld=0
      if(iq.ne.1) elld=ell(iq)*eta(iq)/r(iq)
      elld1=ell(iq1)*eta(iq1)/r(iq1)
      ae=(elld+elld1)*hsq+2.*(ell(iq)-ell(iq1))*hcu
      be=-(2.*elld+elld1)*hr-3.*(ell(iq)-ell(iq1))*hsq

      gd=fot*rho(iq)
      if(iq.ne.1) gd=4.*rho(iq)-2.*g(iq)/r(iq)
      gd1=4.*rho(iq1)-2.*g(iq1)/r(iq1)
      ag=(gd+gd1)*hsq+2.*(g(iq)-g(iq1))*hcu
      bg=-(2.*gd+gd1)*hr-3.*(g(iq)-g(iq1))*hsq


      do 1100 il=1,5

      t=.5*hn*(xi(il)+1.)
      rr=r1+t

      gr=g(iq)+t*(gd+t*(bg+t*ag))
      el=ell(iq)+t*(elld+t*(be+ae*t))
      et=(elld+t*(2.*be+t*3.*ae))*rr/el

      grin(iq,il)=gr
      elin(iq,il)=el
      etin(iq,il)=et

      aa=acon(iq)+t*(qacon(1,iq)+t*(qacon(2,iq)+t*qacon(3,iq)))
      cc=ccon(iq)+t*(qccon(1,iq)+t*(qccon(2,iq)+t*qccon(3,iq)))
      ff=fcon(iq)+t*(qfcon(1,iq)+t*(qfcon(2,iq)+t*qfcon(3,iq)))
      ll=lcon(iq)+t*(qlcon(1,iq)+t*(qlcon(2,iq)+t*qlcon(3,iq)))
      nn=ncon(iq)+t*(qncon(1,iq)+t*(qncon(2,iq)+t*qncon(3,iq)))
      rrho=rho(iq)+t*(qrho(1,iq)+t*(qrho(2,iq)+t*qrho(3,iq)))
      rhoin(iq,il)=rrho

      etan=ff/(aa-2.*ll)

      if(iq.le.nic) then ! innner core 
        do i=0,MLL
          do j=1,2*i+1
            da(i,j,iq,il)=0.
            dc(i,j,iq,il)=0.
            dl(i,j,iq,il)=0.
            dn(i,j,iq,il)=0.
            df(i,j,iq,il)=0.
            dr(i,j,iq,il)=0.
          enddo
        enddo
      else if(iq.le.noc) then ! outer core
        do i=0,MLL
          do j=1,2*i+1
            da(i,j,iq,il)=0.
            dc(i,j,iq,il)=0.
            dl(i,j,iq,il)=0.
            dn(i,j,iq,il)=0.
            df(i,j,iq,il)=0.
            dr(i,j,iq,il)=0
          enddo
        enddo
      else if(iq.le.moho) then ! mantle


        rx=(2.*rr-r(noc)-r(moho))/(r(moho)-r(noc))  ! -1 at cmb 1 at moho
        do i=0,MLL
          do j=1,2*i+1

            dvsrel=0.
            dvprel=0.
            drrel=0.
            do k=4,24
              dvsrel = dvsrel + dvsr(i,j,k)*splh(k-4,rx)               
              dvprel = dvprel + dvpr(i,j,k)*splh(k-4,rx) 
              drrel  = drrel  + drhor(i,j,k)*splh(k-4,rx)
            enddo


            da(i,j,iq,il)=2.*aa*dvprel+aa*drrel
            dc(i,j,iq,il)=2.*cc*dvprel+cc*drrel
            dn(i,j,iq,il)=2.*nn*dvsrel+nn*drrel
            dl(i,j,iq,il)=2.*ll*dvsrel+ll*drrel
            df(i,j,iq,il)=etan*(da(i,j,iq,il)-2.*dl(i,j,iq,il))
            dr(i,j,iq,il)=rrho*drrel


          enddo
        enddo
      else                      ! crust
      endif
 1100 continue
 1000 continue

      return
      end
