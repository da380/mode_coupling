program mdcpl

  use nrtype


  ! dimension information stored in a header file
  include 'dimension.h'

  

!  complex(spc), dimension(MDIM,MDIM) :: A0,A1,A2
  complex(dpc), dimension(MDIM,MDIM) :: A0,A1,A2
  complex(dpc), dimension(MDIM) :: ev,w
  complex(dpc), dimension(2*MDIM) :: wk
!  complex(dpc), dimension(-ML:ML,-ML:ML) :: z
  complex(dpc), dimension(-ML:ML,-ML:ML) :: z0
  complex(dpc), dimension(-ML:ML,-ML:ML) :: z1
  complex(dpc), dimension(-ML:ML,-ML:ML) :: z2


  complex(dpc) :: w0,w1,w2,dsw1,dsw2,wsum

  ! w0 is the fiducial frequency used in the calculations

  integer(i4b), dimension(MMODES) :: lord,nord

  real(dp), dimension(MDIM) :: wr,wq 
 
  complex(dpc) :: ctmp
  
  character(len=1), dimension(MMODES) :: iq
  character(len=80) :: file,getunx
  character(len=256) :: afile
  

  ! common blocks for the eigenfunctions U, U', V, V',ph1,ph1'
  ! jcom=1 --> radial mode
  ! jcom=2 --> toroidal mode
  ! jcom=3 --> spheroidal mode
  ! n = radial order
  ! l  = angular order
  ! om = angular frequency in rad/s
  ! cgp = group velocity in km/s
  ! avert, ahor = vertical and horizontal accelerations at the ocean floor
  !               -- used to calculate seismograms.  
  common/eig1/n1,jcom1,l1,om1,q1,cgp1,avert1,ahor1,phis1 & 
       ,u1(222),up1(222),v1(222),vp1(222),ph1(222),php1(222)
  common/eig2/n2,jcom2,l2,om2,q2,cgp2,avert2,ahor2,phis2 & 
       ,u2(222),up2(222),v2(222),vp2(222),ph2(222),php2(222)
!  common/matr/w0,z
  common/matr0/z0,z1,z2




  ! conversion factors returned by subroutine model
  ! rn =     length conversion to m (=6371000.00)
  ! vn =     velocity conversion to m/s
  ! wn =     frequency conversion to s**(-1)
  ! gn =     acceleration conversion to m/s/s
  ! rhobar = density conversion to kg/m**3
  common/nond/rn,wn,vn,gn,rhobar



  
 
!  call chekcl('|-lu2:o:1:[/home/seismic/john/dta/foanis05.222]'                   & 
!       //'|-lu7:o:1:[mdcpl.out]'                                                  &  
!       //'|-lu3:o:1:[/home/seismic/john/dta/m1084x2.htm] model on unit 3 (rdmdl)' & 
!       //'|-model:o:1:[/home/seismic/john/dta/S20RTS.sph] Model'                  & 
!       //'|-pc:o:1:[default.pc] startup plotting commands'                        &
!       //'|')


  
  
  call chekcl('|-lu2:o:1:[/home/david/coupling_standalone/data/foanis05.222]'                   & 
       //'|-lu7:o:1:[mdcpl.out]'                                                                     &  
       //'|-lu3:o:1:[/home/david/coupling_standalone/data/m1084x2.htm] model on unit 3 (rdmdl)' & 
       //'|-model:o:1:[/home/david/coupling_standalone/data/S20RTS.sph] Model'                  & 
       //'|-pc:o:1:[default.pc] startup plotting commands'                                           &
       //'|')



  
  file = getunx('-lu2',1,ll)
  open(2,file=file,status='old')
  file=getunx('-lu7',1,ll)
  open(7,file=file)
  
  file=getunx('-lu3',1,ll)
  open(3,file=file)

  
  ! open the PREM mode catalogue
  
  call openfl(1,'/home/david/coupling_standalone/data/PREM222.BIN',1,0,0,istat,5364)
  call seteig(1)


   
  ! read in the PREM model and set it up in common blocks modl1, modl2
  ! dimensional conversion factors rn, wn, vn, gn, rhobar
  ! are returned. These are conversion factors from a dimensionless
  ! system of units such that
  !   earth's mean density = 1, pi*G = 1, earth radius=1.
  ! The model is stored in dimensionless units.  
  call modl(2,0,rn,wn,vn,gn,rhobar)




!     I guess fus, fup, flp, fls
!     are factors associated with
!     upper mantle S and P, lower mantle S and P
!     They could be automatically set to 1.0
!
!  -> these values are not needed anymore, replaced by 
!     alfau,alfal,betau,betal.
!     fus is still needed in some subroutines, I think it just
!     has to be set to 1.0
!
  fus=1.0
  !Restored by APV
  call intpltnnew(getunx('-model',1,ll))
  
  !! alternate version that reads in separate vs, vp, and rho
  !! models
  !
  !call intpltnnew_sep(
  !!/home/seiraid2/davida/dta/Model_MCM_M2_dvs.sph', & 
  !!'/home/seiraid2/davida/dta/Model_MCM_M2_dvp.sph', & 
  !!'/home/seiraid2/davida/dta/Model_MCM_M2_drho.sph')
       


  write(6,"('Input mins & maxs  eg:0 8')")  
  read(5,*)mins, maxs
  write(6,"('do you want include rotation? yes=1 no=0')")
  read(5,"(i1)")isw1
  write(6,"('do you want include ellip.? yes=1 no=0')")
  read(5,"(i1)")isw2
  write(6,*) " afile = "
  read(5,*) afile





  ! Read in list of modes from input file

  open(25,file='modes.in')
  rewind(25)
  read(25,*) mtot      
  wsum=cmplx(0.,0.)
  nelem=0
  do im=1,mtot
     read(25,2) nord(im),iq(im),lord(im)
     call fetch(nord(im),iq(im),lord(im),1,ifexst)     
     w1=cmplx(om1,.5*om1*q1)
     wsum=wsum+w1
     nelem=nelem+2*lord(im)+1
  end do
  wsum=wsum/mtot
2 format(i3,1x,a1,1x,i3)



  
!  print *, nelem,MDIM
!  stop


  !  get the eigenfunctions from the mode catalogue
  !  and compute the coupling matrices  
  imbeg=0
  do im=1,mtot


     ! Self-coupling


     print *,  'self-coupling for mode: ',nord(im),iq(im),lord(im)

     call fetch(nord(im),iq(im),lord(im),1,ifexst)
     call fetch(nord(im),iq(im),lord(im),2,ifexst)
     if(ifexst.ne.1) then
        write(6,"('mode does not exist in catalogue')")
        write(6,*) 'In self-coupling', nord(im),iq(im),lord(im),ifexst
        stop
     end if


     w1=cmplx(om1,.5*om1*q1)
     w0=om1



     call memdcpl_new(isw1,isw2,fus,mins,maxs)



    do i=1,2*lord(im)+1
        do j=1,2*lord(im)+1
           A0(imbeg+i,imbeg+j)=z0(i-lord(im)-1,j-lord(im)-1)
           A1(imbeg+i,imbeg+j)=z1(i-lord(im)-1,j-lord(im)-1)
           A2(imbeg+i,imbeg+j)=z2(i-lord(im)-1,j-lord(im)-1)
        end do
     end do
        

     

     do i=1,2*lord(im)+1
        A0(imbeg+i,imbeg+i)=A0(imbeg+i,imbeg+i)+w1*w1
     end do
     
     
     if(im == mtot) exit
     

     ! Cross-coupling
     
     jmbeg=imbeg+2*lord(im)+1
     do jm=im+1,mtot
        call fetch(nord(jm),iq(jm),lord(jm),2,ifexst)
        if(ifexst.ne.1) then
           write(6,"('mode does not exist in catalogue')")
           write(6,*) 'In cross coupling', nord(jm),iq(jm),lord(jm),ifexst
           stop
        end if

        w2=cmplx(om2,.5*om2*q2)
        w0=0.5*(om1+om2)

        print *,  'cross-coupling for modes: ',nord(im),iq(im),lord(im)
        print *,  '                          ',nord(jm),iq(jm),lord(jm)

        call memdcpl_new(isw1,isw2,fus,mins,maxs)


        do i=1,2*lord(im)+1
           do j=1,2*lord(jm)+1
              A0(imbeg+i,jmbeg+j)=z0(i-lord(im)-1,j-lord(jm)-1)
              A0(jmbeg+j,imbeg+i)=conjg(A0(imbeg+i,jmbeg+j))
              A1(imbeg+i,jmbeg+j)=z1(i-lord(im)-1,j-lord(jm)-1)
              A1(jmbeg+j,imbeg+i)=conjg(A1(imbeg+i,jmbeg+j))
              A2(imbeg+i,jmbeg+j)=z2(i-lord(im)-1,j-lord(jm)-1)
              A2(jmbeg+j,imbeg+i)=conjg(A2(imbeg+i,jmbeg+j))
           end do
        end do        
        jmbeg=jmbeg+2*lord(jm)+1
     end do
     imbeg=imbeg+2*lord(im)+1
     
  end do



  print *, ' writing out coupling matrices'

  print *, 'nelem = ',nelem

  ! normalize the A2 matrix
  a2 = a2/wn**2


  ! write out matrix parts in binary form
  open(35,file=trim(afile),form='unformatted')
  write(35) a0(1:nelem,1:nelem)
  write(35) a1(1:nelem,1:nelem)
  write(35) a2(1:nelem,1:nelem)
  close(35)



  
end program mdcpl
