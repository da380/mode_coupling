
c-----------------------------------------------
c This routine evaluates the harmonic coefficients
c of delta(A), delta(C) etc. using the input heterogeneous
c models. It will need to be changed to accomodate different
c model parameterisations etc.
c

      subroutine intpltnnew(hetmodel)

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

      character*(*) hetmodel
      character*2 ftype
      common/nmdl/lmaxn,nstr,nlp,anmdl(MXMDLL),bnmdl(MXMDLL)
      dimension mask(MXPARM),mask1(MXPARM),lask(0:MLL),lask1(0:MLL)
      dimension dvsr(0:MLL,MM,MXPARM)

c Legendre polynomials

      p1(x)=sqrt(1.5)*x
      p2(x)=sqrt(5./8.)*(3.*x*x-1.)
      p3(x)=sqrt(7./8.)*(5.*x**3-3.*x)
      p4(x)=sqrt(9./128.)*(35.*x**4-3.*x**2+3.)
      p0=sqrt(.5)



      call rdmdl(3)

      open(27,file=hetmodel)
      rewind(27)
      read(27,'(a)') ftype
      close(27)

      if(ftype.eq.'LG') then

ca-------------------------------------------------
ca Read in mantle model in Legendre polynomials
ca-------------------------------------------------

!      write(6,*) 'File type is Legendre polynomials'
      lu=20
      open(lu,file=hetmodel)
      rewind(lu)
      read(lu,'(a)') ftype
      read(lu,*) lmax,kmax,ieven

      do k=1,kmax
         do i=0,lmax
            do j=1,2*i+1
            if(mod(i,2).eq.0) then
               read(lu,*) dvsr(i,j,k)
            else
              if(ieven.eq.1) then
                read(lu,*) dvsr(i,j,k)
              else
                dvsr(i,j,k)=0.0
              endif
            endif
            enddo
         enddo
         if(lmax.lt.MLL) then
         do i=lmax+1,MLL
            do j=1,2*i+1
               dvsr(i,j,k)=0.0
            enddo
         enddo
         endif
      enddo
      close(lu)

      else

ca-------------------------------------------------
ca Read in mantle model in splines (John's format)
ca-------------------------------------------------

!      write(6,*) 'File type is splines (John format)'
c1      write(6,'(a)') 'Reading model: '//hetmodel
      call heread(1,hetmodel,bnmdl,nstr,lmaxn,mask,lask,0)
c1      write(6,'(a,i5,a,i5,a,i5)') 'Model read into array bnmdl. nstr=', nstr
c1     1              ,' lmaxn=',lmaxn
      do i=1,nstr
         mask1(i)=1
      enddo
      do i=0,lmaxn
         lask1(i)=1
      enddo
      call mskmdl(bnmdl,nstr,lmaxn,mask,lask
     1       ,anmdl,nstr,lmaxn,mask1,lask1) ! expanded model now in anmdl
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
            dvsr(i,j,k) = anmdl(leny*(k-1)+ind)  ! actually dvs/vs
         enddo
        enddo
      enddo

      endif
      print *, n
      !APV altered upper limit from n to n-1
      !cf emails with DA-A, 15/Jan/2015
      do 1000 iq=1,n-1  ! do loop over PREM levels iq=1 (centre) to iq=n (surface)
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
c
c
c        f1=2.*cc*alfau*sqrt(rrho/ll)
c        f2=2.*sqrt(ll*rrho)
c        f3=betau*rrho*sqrt(rrho/ll)
        rx=(2.*rr-r(noc)-r(moho))/(r(moho)-r(noc))  ! -1 at cmb 1 at moho
        do i=0,MLL
          do j=1,2*i+1

            if(ftype.eq.'LG') then    ! Legendre polynomials	
            dvsrel=(p0*dvsr(i,j,1)+p1(rx)*dvsr(i,j,2)+p2(rx)*dvsr(i,j,3)+
     1         p3(rx)*dvsr(i,j,4)+p4(rx)*dvsr(i,j,5))
            else                      ! Spline format
            dvsrel=0.
            do k=4,24
              dvsrel=dvsrel+dvsr(i,j,k)*splh(k-4,rx) ! coef of dvs/vs
            enddo
            endif
c
c           figure out adjustments to A,C,L,N,F
c           assuming (because it's built into splh format that 
c           dvp/vp = .5*dvs/vs
c           drho/rho = 0
c


!           turn off the heterogeneity
!            dvsrel = 0.0
            


            dvprel=.5*dvsrel  
            drrel=0.3*dvsrel  
!            turn of the density perturbation
!            drrel=0.0              
!            dvsrel = 0.0


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
