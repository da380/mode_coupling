      program mode_list
c     this program generates a list of those modes
c     from a given mode catelog with frequency lying
c     in given range



      parameter (MMODES=20000)

      parameter (pi = 3.141592653589793238462643383279502884197)

      
      integer lord(MMODES),nord(MMODES),it(MMODES),iwksp(MMODES)

      real w(MMODES)

      real wmin,wmax

      character*1 imt
      character*1 iq(MMODES)
      character*80 file
      character*80 getunx
      character*160 datapath

      common/eig1/n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1
     1        ,u1(222),up1(222),v1(222),vp1(222),ph1(222),php1(222)
      common/eig2/n2,jcom2,l2,om2,q2,cgp2,avert2,ahor2,phis2
     1        ,u2(222),up2(222),v2(222),vp2(222),ph2(222),php2(222)


c     normalization parameters
      common/nond/rn,wn,vn,gn,rhobar

c     common blocks for the PREM model
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(222)
     1            ,rho(222),qrho(3,222),g(222),ell(222),eta(222)
      common/modl2/acon(222),qacon(3,222),ccon(222),qccon(3,222)
     1            ,lcon(222),qlcon(3,222),ncon(222),qncon(3,222)
     2            ,fcon(222),qfcon(3,222)


      integer*2    indf(1000),nrec,maxn,maxl
      common/get/nrec(500,2),maxn(500,2),maxl(2),llu

c=====================================================================c
c=====================================================================c



!      call chekcl('|-lu2:o:1:[/home/seismic/john/dta/foanis05.222]'
!     1        //'|-lu7:o:1:[mdcpl.out]'
!     1         //'|-lu3:o:1:[/home/seismic/john/dta/m1084x2.htm]'
!     1         //'model on unit 3 (rdmdl)'
!     1         //'|-model:o:1:[/home/seismic/john/dta/S20RTS.sph] Model'
!     1         //'|-pc:o:1:[default.pc] startup plotting commands'
!     1         //'|')


      call chekcl('|-p:o:1:[../data]'                    
     1     //'|-lu2:o:1:[foanis05.222]'                                                                       
     1     //'|-lu7:o:1:[mdcpl.out]'                                                                       
     1     //'|-lu3:o:1:[m1084x2.htm] model on unit 3 (rdmdl)'  
     1     //'|-model:o:1:[S20RTS.sph] Model'                   
     1     //'|-ref:o:1:[PREM222.BIN] reference model'                   
     1     //'|-pc:o:1:[default.pc] startup plotting commands'                                           
     1     //'|')

c     location of files
      datapath=getunx('-p',1,lpath)

c     this opens the PREM model file default foanis05.222
      file=datapath(1:lpath)//'/'//getunx('-lu2',1,ll)
      open(2,file=file,status='old')



c     this opens the outout file mdcpl.out
      file=getunx('-lu7',1,ll)
      open(7,file=file)




c     open the PREM catelog
      file=datapath(1:lpath)//'/'//getunx('-ref',1,ll)
      call openfl(1,file,1,0,0,
     1     istat,5364)
!      call openfl(1,'/home/eeyore1/john/dta/PREM.BIN',1,0,0,
!     1    istat,5364)
      call seteig(1)



c     read the reference model
      call modl(2,0,rn,wn,vn,gn,rhobar)


!      open(99,file='mod.out') 
!      do mm = 1,222
!         write(99,*) r(mm),rho(mm),g(mm),ell(mm),eta(mm)
!      end do
!      close(99)




c     get the minimum frequency
      print *, ' minimum frequency (mHz):'
      read(5,*) wmin
c      wmin = 0.1
      wmin = pi*wmin/500.



c     get the maximum frequency
      print *, ' maximum frequency (mHz):'
      read(5,*) wmax
c      wmax = 4.
      wmax = pi*wmax/500.





c     do the spheroidal modes
      ks = 0
      k = 0
      nss = 0
      lmax = 0
      do ind = 1,maxl(1)
         l = ind-1
         do jnd = 1,maxn(ind,2)
            n = jnd -1
            call fetch_dum(n,'s',l,1,ifexst)
            if(ifexst .eq. 1) then
               if(om1 .ge. wmin .and. om1 <= wmax) then
                  ks = ks+1
                  k = k+1
                  lord(k) = l1
                  iq(k) = 's'
                  nord(k) = n1
                  it(k) = 3
                  w(k) = om1                  
                  nss = nss+(2*l1+1)
                  if(l1 > lmax) lmax = l1
!                  open(99,file='eig.out')
!                  do mm = 1,222
!                     write(99,*) r(mm),u1(mm),v1(mm),ph1(mm)
!                     print *,  n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1
!     1                    ,u1(mm),up1(mm),v1(mm),vp1(mm),ph1(mm),php1(mm)
!                  end do
!                  close(99)
               end if
            end if
         end do
      end do




c     do the toroidal modes
      kt = 0
      nts = 0
      do ind = 1,maxl(2)
         l = ind-1
         do jnd = 1,maxn(ind,1)
            n = jnd -1
            call fetch_dum(n,'t',l,1,ifexst)
            if(ifexst .eq. 1) then
               if(om1 .ge. wmin .and. om1 <= wmax) then
                  kt = kt+1
                  k = k+1
                  lord(k) = l1
                  iq(k) = 't'
                  nord(k) = n1
                  it(k) = 2
                  w(k) = om1
                  nts = nts+(2*l1+1)
                  if(l1 > lmax) lmax = l1
               end if
            end if
         end do
      end do

      print *, ' number of spheroidal multiplets = ',ks
      print *, ' number of spheroidal singlets = ',nss
      print *, ' number of toroidal multiplets = ',kt
      print *, ' number of toroidal singlets = ',nts
      print *, ' number of multiplets = ',k
      print *, ' number of singlets = ',nss+nts
      print *, ' maximum degree = ',lmax



c     index the list of modes according to frequency
      call indexx(k,w,iwksp)

c     write out the list of modes
      open(7,file='modes.in',form='formatted')
      write(7,*) k
      do i = 1,k
         write(7,2) nord(iwksp(i)),iq(iwksp(i)),lord(iwksp(i))
      end do
    2 format(i3,1x,a1,1x,i3)      
      close(7)


c     write out all the frequencies to a file for plotting
      open(7,file='modes.in.freq',form='formatted')
      write(7,*) ks,kt,0
      do i = 1,ks
         write(7,*) nord(i),lord(i),w(i)
      end do
      do i = ks+1,k
          write(7,*) nord(i),lord(i),w(i)
      end do
      close(7)

      end
