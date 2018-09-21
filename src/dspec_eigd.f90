program dspec_eig


  use nrtype
  use module_util
  use module_fourier
  implicit none

  logical(lgt) :: ltmp

!  integer(i4b), parameter :: ML = 75
!  integer(i4b), parameter :: MDIM = 2500
!  integer(i4b), parameter :: MMODES = 500


!  integer(i4b), parameter :: ML     = 150
!  integer(i4b), parameter :: MDIM   = 5000
!  integer(i4b), parameter :: MMODES = 1000

  integer(i4b), parameter :: ML     = 300
  integer(i4b), parameter :: MDIM   = 10000
  integer(i4b), parameter :: MMODES = 2000
 
  integer(i4b), parameter :: maxit = 50
  
  character(len=1), dimension(:), allocatable :: iq
  character(len=256) :: matrix_file
  character(len=256) :: spec_out 
  character(len=256) :: ts_out
  character(len=256) :: source_file
  character(len=256) :: rec_file

  integer(i4b), parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer(i4b) :: mtot,nelem,im,i,j,info,istat, & 
       nindex,ifgot,nw,iw,ir1,nt,i1,i2,mex,qex,nt0,ne, & 
       ntb,md,im,iter,mindex,ib,jb,tbd,it,count,ntime, & 
       acountit,jm

  integer(i4b), dimension(:), allocatable :: nn,ll,ity

  integer(i4b), dimension(:), allocatable :: ipivot,ipivotd

  integer(i4b), dimension(:,:), allocatable :: tbl

  real(sp), dimension(9) :: source
  real(sp), dimension(:,:), allocatable :: station


  real(dp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
       t1,t2,atimeit,timeit1,timeit2,atimeitt
  
  real(dp) :: err
  real(dp), parameter :: tol = 1.e-4_sp

  complex(dpc) :: w,aclr,aclr0,cfac,w0,exct,cexp
  complex(dpc) :: alpha,beta,tmp1,tmp2
  complex(dpc), dimension(:), allocatable :: ww
  complex(dpc), dimension(:,:), allocatable :: a0,a1,a2
  complex(dpc), dimension(:,:), allocatable :: ad,a2i,utmp

  complex(dpc), dimension(:,:), allocatable :: a,atb
  complex(dpc), dimension(:,:), allocatable :: vsd
  complex(dpc), dimension(:), allocatable :: xs,xr,x0

  complex(dpc), dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, & 
                                              sb1,p1,pb1,vt
  integer(i4b) :: saved 
  complex(dpc), dimension(:), allocatable :: xsave
  
  complex(dpc), dimension(:), allocatable :: vs
    complex(dpc), dimension(:,:), allocatable :: vr
  complex(dpc), dimension(2*ML+1) :: vrt,vst

  complex(dpc), dimension(:), allocatable :: acl
  complex(dpc), dimension(:,:), allocatable :: dat,dtmp


  complex(dpc), dimension(:), allocatable :: wwh
  complex(dpc), dimension(:,:), allocatable :: uui,uu
  complex(dpc), dimension(:), allocatable :: xfac

  complex(dpc), dimension(:), allocatable :: work
  real(dp), dimension(:), allocatable :: work2


  real(sp), dimension(12) :: amp
  real(sp), dimension(10) :: aker
  real(sp), dimension(4)  :: ar1,ar2
  real(sp) :: r1,r2
  integer(i4b) :: nord,it,lord

  ! common block for the eigenfunctions parameters
  common/premdata/amp(12),aker(10),ar1(4),ar2(4),r1,r2,nord,it,lord
  
  ! Read in list of modes from input file
  open(io1,file='modes.in')
  rewind(io1)
  read(io1,*) mtot      

  ! allocate some arrays
  allocate(iq(mtot),nn(mtot),ll(mtot),ity(mtot))
  nelem=0
  do im=1,mtot
     read(io1,2) nn(im),iq(im),ll(im)
     nelem=nelem+2*ll(im)+1
     if(iq(im) == 's') then
        ity(im) = 'S'
     else
        ity(im) = 'T'
     end if
  end do
2 format(i3,1x,a1,1x,i3)
  
  close(io1)

  ! allocate everything else
  allocate(ipivot(nelem))
  allocate(ipivotd(2*nelem))
  allocate(tbl(mtot,4))
  allocate(ww(mtot))
  allocate(a0(nelem,nelem),a1(nelem,nelem), &        
       a2(nelem,nelem),a(nelem,nelem),atb(nelem,nelem))
  allocate(vsd(nelem,1))
  allocate(xs(nelem),xr(nelem),x0(nelem))
  allocate(s0(nelem),sb0(nelem),p0(nelem),pb0(nelem), & 
       s1(nelem),sb1(nelem),p1(nelem),pb1(nelem),vt(nelem))
  allocate(xsave(nelem))
  allocate(vr(nelem,nr),vs(nelem))


  allocate(ad(2*nelem,2*nelem),a2i(nelem,nelem),utmp(2*nelem,2*nelem))
  
  allocate(wwh(2*nelem))
  allocate(uui(2*nelem,2*nelem),uu(2*nelem,2*nelem))
  allocate(work(4*nelem))
  allocate(work2(4*nelem))



  ! get parameters for the calculations
  call get_float(' f1 (mhz) = ',f1)
  call get_float(' f2 (mhz) = ',f2)            
  call get_float(' dt (sec) = ',dt)
  call get_float(' tout (hrs) = ',tout)
  call get_float(' df0 (mhz) = ',df0)
  call get_float(' wtb (mhz) = ',wtb)
  call get_float(' t1 (hrs) = ',t1)
  call get_float(' t2 (hrs) = ',t2)  
  call get_string(' spectra file = ',spec_out)
  call get_string(' time-series file = ',ts_out)
  call get_string( ' matrix file = ',matrix_file)
  
  ! read in the coupling matrices

  open(io1,file=trim(matrix_file),form='unformatted')
  read(io1) a0
  read(io1) a1
  read(io1) a2
  close(io1)


  ! read in source parameters
  call get_string(' source file = ',source_file )
  open(99,file=trim(source_file))
  do i = 1,9
     read(99,*) source(i)
  end do
  close(99)


  ! normalize the moment tensor so that the 
  ! output displacements are in cm's
  source(4:9) = source(4:9)*0.95179e-30_sp


  ! receiver location
  call get_string(' receiver file = ',rec_file )
  open(99,file=trim(rec_file))
  read(99,*) nr
  allocate(station(5,nr))
  allocate(rec_string(nr))
  do j = 1,nr
        read(99,*) station(1:5,j),rec_string(j)
  end do
  close(99)

  
  ! open the mode catalog
  call openfl(7,'../data/SPRM1.BIN',1,0,0,istat,4096)
  
  ! set up mode catalog
  call setprem(7,0)
  ! find the depth of the source
  call findd(source(3),ir1)
   
  !  start the loop over the modes      
  nindex = 0
  do im = 1,mtot
     
     ! get the modes
     call getprem(nn(im),ity(im),ll(im),ir1,ifgot)
          
     ! modes complex frequency
     ww(im) = cmplx(amp(1),0.5*amp(1)*amp(2))   


      if(im ==1) then
         fmin = amp(1)/twopi_d
         fmax = fmin
      else
         if(amp(1)/twopi_d < fmin) fmin = amp(1)/twopi_d 
         if(amp(1)/twopi_d > fmax) fmax = amp(1)/twopi_d 
      end if
   
     ! compute the source vector
     call csveca(vst,source)
          
     do i=1,2*ll(im)+1
        vs(nindex+i)=vst(i)
     enddo

     ! begin loop over the recievers
     do ir = 1,nr
          
        ! compute the receiver vector 
        call crvecd(vrt,station(:,ir))
        
        do i=1,2*ll(im)+1
           vr(nindex+i,ir)=vrt(i)
        enddo

     end do
     ! end loop over the receivers
          
     nindex = nindex+2*ll(im)+1
     
  end do
  ! end loop over the modes


  ! invert the kinetic energy matrix
  a2i = 0.0_sp
  do i = 1,nelem
     a2i(i,i) = 1.0_sp
  end do
  call zgetrf(nelem,nelem,a2,nelem,ipivot,info)
  call zgetrs('N',nelem,nelem,a2,nelem,ipivot,a2i,nelem,info)



  ! put together the double sized matrix
  ad = 0.0_dp
  do i = 1,nelem
     ad(i,i+nelem) = 1.0_sp
  end do
  ad(nelem+1:2*nelem,1:nelem) = -matmul(a2i,a0)
  ad(nelem+1:2*nelem,nelem+1:2*nelem) = -matmul(a2i,a1)

  
  ! perform the eigenvalue decomposition of double sized matrix
  call zgeev('v','v',2*nelem,ad,2*nelem,wwh,uui,2*nelem,uu,2*nelem, & 
       work,4*nelem,work2,info)

  ! invert the  matrix of right-eigenvectors
  utmp = uu
  call zgetrf(2*nelem,2*nelem,utmp,2*nelem,ipivotd,info)
  uui = 0.0_sp
  do i = 1,2*nelem
     uui(i,i) = 1.0_dp
  end do
  call zgetrs('N',2*nelem,2*nelem,utmp,2*nelem,ipivotd,uui,2*nelem,info)


  ! get the normalization factors
!  do j = 1,nelem
!     xfac(j) = 0.0_sp
!     do i = 1,nelem
!        xfac(j) = xfac(j) + uui(i,j)*uu(j,i)
!     end do
!     print *, xfac(j)
!  end do

  ! get the eigenvalues
  wwh = sqrt(wwh)
  where(real(wwh) < 0.0_sp) wwh = -real(wwh) +ii*abs(imag(wwh))


  ! print out the eigenvectors
  open(99,file='eigs.out')  
  nindex = 0
  do im = 1,mtot
     do i = 1,2*ll(im)+1
        write(99,*) 1000.0_dp*real(wwh(nindex+i))/twopi
     end do
     nindex = nindex+2*ll(im)+1
  enddo
  close(99)


  ! transform the source vector
  vs = matmul(uui(1:nelem,1:nelem),vs)



  ! begin loop over the recievers
  do ir = 1,nr

     ! transform the receiver vector  
     vr(:,ir) = matmul(transpose(uu(1:nelem,1:nelem)),vr(:,ir))
  end do
  

  ! get parameters for the calculations
  call get_float(' f1 (mhz) = ',f1)
  call get_float(' f2 (mhz) = ',f2)            
  call get_float(' dt (sec) = ',dt)
  call get_float(' tout (hrs) = ',tout)
  call get_float(' df0 (mhz) = ',df0)
  call get_float(' t1 (hrs) = ',t1)
  call get_float(' t2 (hrs) = ',t2)
  
  call get_string(' spectra file = ',spec_out)
  call get_string(' time-series file = ',ts_out)
  
  f1 = f1/1000.0_sp
  f2 = f2/1000.0_sp
  df0 = df0/1000.0_sp
  tout = tout*3600.0_sp
  t1 = t1*3600.0_sp
  t2 = t2*3600.0_sp
  if(t2 > tout) t2 = tout



 
  ! work out frequency step size etc
  mex = 5
  qex = 4 
  call fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
  print *, ' f1 = ',f1*1000.0_sp
  print *, ' f2 = ',f2*1000.0_sp
  print *, ' df = ',df*1000.0_sp
  print *, ' nt = ',nt

  wr1 = twopi_d*f1
  wr2 = twopi_d*f2
  dwr = twopi_d*df
  wi = ep


 ! initialize the output arrays
  allocate(acl(nt))  


  acl = 0.0_sp
  do i = 1,nelem

     exct = vs(i)*vr(i)
     cexp = exp(ii*wwh(i)*dt)

     acl(1) = acl(1) + exct

     do it = 2,nt
        
        exct = exct*cexp

       acl(it) = acl(it) + exct
        
     end do

  end do

  acl = real(acl)

  allocate(dat(1,nt))
  dat(1,:) = acl(:)


  ! parameters for time-domain filter
  fac = 0.5_sp
  t11 = t1
  t22 = t2
  t12 = t11+fac*(t22-t11)
  t21 = t22-fac*(t22-t11)
  
  ! write out the time-series
  open(io1,file=trim(ts_out))
  do i = 1,nt
     t = (i-1)*dt
     filt =  hann(t,t11,t12,t21,t22)
     if(t > tout) exit
!     write(io1,*) t,filt*real(dat(1,i))
     write(io1,*) t,real(dat(1,i))
  end do
  close(io1)
  
  ! filter the time series
  do i = 1,nt
     t = (i-1)*dt
     filt = hann(t,t11,t12,t21,t22)
     dat(1,i) = filt*dat(1,i)
  end do
    
  ! pad the time series if necessary
  nt0 = floor(1.0_sp/df0*dt)
  if(nt0 > nt) then
     ne = log(real(nt0))/log(2.0_sp)+1
     nt0 = 2**ne
     allocate(dtmp(1,nt))
     dtmp = dat
     deallocate(dat)
     allocate(dat(1,nt0))
     dat(1,1:nt) = dtmp(1,:)
     dat(1,nt+1:nt0) = 0.0_sp
     nt = nt0
     df = 1.0_dp/(nt*dt)
     i1 = max(floor(f1/df),1)
     i2 = floor(f2/df)+1
     deallocate(dtmp)
  end if
  
  ! do the Fourier transform
  call fourrow_dp(dat,-1)
  dat = dat*dt



  ! parameters for frequency-domain filter
  fac = 0.05_sp
  f11 = f1
  f22 = f2
  f12 = f11+fac*(f22-f11)
  f21 = f22-fac*(f22-f11)

    do i = 1,nt/2+1
     f = (i-1)*df
     filt = hann(f,f11,f12,f21,f22)
!     dat(1,i) = filt*dat(1,i)
  end do

  ! converts displacements into cm (according to Arwen's program)
  acl(:)  = dat(1,:)
  
  ! write out the spectra
  open(io1,file=trim(spec_out)) 
  do i = i1,i2
     f = (i-1)*df*1000.0_sp
     write(io1,'(3e20.12e3)') f,real(acl(i)),imag(acl(i))
!     write(io1,'(2e20.12e3)') f,abs(acl(i))
  end do
  
  close(io1)
  


contains


  subroutine fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    use nrtype
    implicit none
    real(sp), intent(inout) :: f1,f2
    real(sp), intent(in) :: dt,tout
    integer(i4b), intent(in) :: mex,qex
    real(sp), intent(out) :: df,ep
    integer(i4b), intent(out) :: nt,i1,i2
    
    integer(i4b) :: ne
    real(sp) :: fn
    
    ! check that the Nyquist frequency for the given
    ! value of dt lies above f2
    fn = 0.5_sp/dt
    if(fn < f2) stop ' f2 is greater than the Nyquist frequency for the time step'
    
    ep = mex/tout
    
    df = ep/(twopi_d*qex)
    
    nt = 1.0_sp/(df*dt)
    
    ne = log(real(nt))/log(2.0_sp)+1
    nt = 2**ne
    
    df = 1.0_sp/(nt*dt)
    
    i1 = max(floor(f1/df),1)
    f1 = (i1-1)*df          
    
    i2 = floor(f2/df)+2
    f2 = (i2-1)*df

    return
  end subroutine fcal

  function hann(t,t11,t12,t21,t22)
    use nrtype
    implicit none
    logical(lgt) :: ltmp
    real(sp) :: hann
    real(sp), intent(in) :: t,t11,t12,t21,t22

    ltmp = (t11 == 0.0_sp .and. t12 == 0.0_sp & 
         .and. t21 == 0.0_sp .and. t22 == 0.0_sp)
    if(.not.ltmp) then
       if(t < t11) then
          hann = 0.0_sp
       else if(t >= t11 .and. t < t12) then
          hann = pi_d*(t-t11)/(t12-t11)
          hann = 0.5_sp*(1.0_dp-cos(hann))
       else if(t >= t12 .and. t < t21) then
          hann = 1.0_dp
       else if(t >= t21 .and. t < t22) then
          hann = pi_d*(t22-t)/(t22-t21)
          hann = 0.5_sp*(1.0_dp-cos(hann))
       else if(t >= t22) then
          hann = 0.0_sp
       end if       
    end if
    return
  end function hann



  

end program dspec_eig
