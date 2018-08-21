

program dspec5

  use nrtype
  use module_util
  use module_fourier
  !use SUNPERF
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
 
  integer(i4b), parameter :: maxit = 1000
  
  character(len=1), dimension(:), allocatable :: iq
  character(len=256) :: spec_out 
  character(len=256) :: ts_out
  character(len=256) :: source_file
  character(len=256) :: rec_file
  character(len=256) :: time_file

  integer(i4b), parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer(i4b) :: mtot,nelem,im,i,j,info,istat, & 
       nindex,ifgot,nw,iw,ir1,nt,i1,i2,mex,qex,nt0,ne, & 
       ntb,md,iter,mindex,ib,jb,tbd,it,count,ntime, & 
       acountit,dimtb,im2,ifcor

  integer(i4b), dimension(:), allocatable :: nn,ll,ity

  integer(i4b), dimension(:), allocatable :: ipivot

  integer(i4b), dimension(:,:), allocatable :: tbl

  real(sp), dimension(9) :: source
  real(sp), dimension(2) :: station

  real(sp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
       t1,t2,atimeit,timeit1,timeit2,atimeitt
  
  real(sp) :: err
  real(sp), parameter :: tol = 1.e-4_sp

  complex(spc) :: w,aclr,aclr0,cfac,w0
  complex(spc) :: alpha,beta,tmp1,tmp2
  complex(spc), dimension(:), allocatable :: ww
  complex(spc), dimension(:,:), allocatable :: a0,a1,a2
  complex(dpc), dimension(:,:), allocatable :: a0d,a1d,a2d

  complex(spc), dimension(:,:), allocatable :: a,atb
  complex(spc), dimension(:,:), allocatable :: vsd
  complex(spc), dimension(:), allocatable :: xs,xr,x0

  complex(spc), dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, & 
                                              sb1,p1,pb1,vt
  integer(i4b) :: saved 
  complex(spc), dimension(:), allocatable :: xsave
  
  complex(spc), dimension(:), allocatable :: vr,vs
  complex(spc), dimension(2*ML+1) :: vrt,vst

  complex(spc), dimension(:), allocatable :: acl
  complex(dpc), dimension(:,:), allocatable :: dat,dtmp

  real(sp), dimension(12) :: amp
  real(sp), dimension(10) :: aker
  real(sp), dimension(4)  :: ar1,ar2
  real(sp) :: r1,r2
  integer(i4b) :: nord,lord

  ! common block for the eigenfunctions parameters
  common/premdata/amp,aker,ar1,ar2,r1,r2,nord,it,lord
  print *, bigg
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
        ity(im) = transfer('S   ',1)
     else
        ity(im) = transfer('T   ',1)
     end if
  end do
2 format(i3,1x,a1,1x,i3)
  
  close(io1)

  ! allocate everything else
  allocate(ipivot(nelem))
  allocate(tbl(mtot,4))
  allocate(ww(mtot))
  allocate(a0(nelem,nelem),a1(nelem,nelem), &        
       a2(nelem,nelem),a(nelem,nelem),atb(nelem,nelem))
  allocate(vsd(nelem,1))
  allocate(xs(nelem),xr(nelem),x0(nelem))
  allocate(s0(nelem),sb0(nelem),p0(nelem),pb0(nelem), & 
       s1(nelem),sb1(nelem),p1(nelem),pb1(nelem),vt(nelem))
  allocate(xsave(nelem))
  allocate(vr(nelem),vs(nelem))
 
  ! read in the coupling matrices 
  open(io1,file='matrix_parts_sp.bin',form='unformatted')
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


  call get_string(' receiver file = ',rec_file )
  open(99,file=trim(rec_file))
  do i = 1,2
     read(99,*) station(i)
  end do
  close(99)

  
  ! open the mode catalog
  call openfl(7,'/home/eeyore_data3/john/dta/SPRM1.BIN', & 
       1 ,0,0,istat,4096)
  
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
          
     ! compute the receiver vector 
     call crveca(vrt,station)
     
     do i=1,2*ll(im)+1
        vr(nindex+i)=vrt(i)
     enddo
     
     nindex = nindex+2*ll(im)+1
     
  end do
  ! end loop over the modes

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
  
  f1 = f1/1000.0_sp
  f2 = f2/1000.0_sp
  df0 = df0/1000.0_sp
  tout = tout*3600.0_sp
  t1 = t1*3600.0_sp
  t2 = t2*3600.0_sp
  if(t2 > tout) t2 = tout
  wtb = wtb*pi_d/500.0_sp
  
  ! work out frequency step size etc
  mex = 5
  qex = 2 
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
  acl(:)  = 0.0_sp

  ntime = 0
  atime = 0
  atimeitt =0
  acount = 0
  saved = 0


  ! switch whether to use the use the corilis approximation
  ifcor = 0


  ! implement the kinetic and coriolis approximations
  if(ifcor == 1) then
     nindex = 0
     do im = 1,mtot
        do i = 1,2*ll(im)+1
           mindex = 0
           do im2 = 1,mtot
              w0 = 0.5_sp*real(ww(im)+ww(im2))
              do j = 1,2*ll(im2)+1
                 a0(nindex+i,mindex+j) = a0(nindex+i,mindex+j)      & 
                      +w0*a1(nindex+i,mindex+j)  & 
                      +w0*w0*a2(nindex+i,mindex+j)  
                 a1(nindex+i,mindex+j) = 0.0_sp
                 a2(nindex+i,mindex+j) = 0.0_sp
              end do
              mindex = mindex + 2*ll(im2)+1
           end do
        end do
        nindex = nindex + 2*ll(im)+1
     end do
  end if

  call cpu_time(time1)




  ! begin loop over frequency
  do iw = i1,i2
     


     wr =(iw-1)*dwr
     w = wr-ii*wi
     
!     print *, ' frequencies to go = ',i2-iw+1

     !----------------------------!
     ! build the coupling matrix  !
     !----------------------------!
     do i = 1,nelem
        do j = 1,nelem
           a(i,j) = a0(i,j)+w*a1(i,j)+w*w*a2(i,j)
        end do
        a(i,i) = a(i,i)-w*w
     end do


     
     

     
     !-----------------------------------!
     ! build the preconditioning matrix  !
     !-----------------------------------!
     
     ! work out which modes are in the target block
     ntb = 0
     nindex = 0
     mindex = 0
     dimtb = 0
     do im = 1,mtot
        if(check_mode(w,ww(im),wtb)) then
           ntb = ntb+1
           tbl(ntb,1) = im         ! mode number for the block
           tbl(ntb,2) = nindex     ! starting index of the sub-block in the full array
           tbl(ntb,3) = mindex     ! starting index of the sub-block in the reduced array
           tbl(ntb,4) = 2*ll(im)+1 ! dimension of the sub-block
           dimtb = dimtb + tbl(ntb,4) 
           mindex = mindex+2*ll(im)+1
        end if
        nindex = nindex+2*ll(im)+1
     end do

     



     

     ! deal with the target block contribution!
     if(ntb /= 0) then
        


        ! dimension of the target system
        tbd = 0
        do ib = 1,ntb
           tbd = tbd + tbl(ib,4)
        end do
        
        ! build the coupling matrix for the 
        ! target block
        do ib = 1,ntb
           do jb = 1,ntb
              atb(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4), & 
                  tbl(jb,3)+1:tbl(jb,3)+tbl(jb,4)) = & 
                  a(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4), & 
                  tbl(jb,2)+1:tbl(jb,2)+tbl(jb,4))    
           end do
        end do
        
        ! perform an LU decomposition of the target-block
        call getrf(tbd,tbd,atb,nelem,ipivot,info)
        
        ! assemble the source vector
        ! note the assumed step-function time-dependence of source
        do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) & 
                = vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))/(ii*w)
        end do

        ! solve the linear system for target block 
        call getrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
        do ib = 1,ntb
           xs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = & 
                vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
        end do
        
     end if
     
     ! deal with the rest of the matrix
     if(ntb /= mtot) then
        
        nindex = 0
        do im = 1,mtot
           md = 2*ll(im)+1
           if(.not.check_mode(w,ww(im),wtb)) then
              cfac = 1.0_sp/(ww(im)**2-w**2)
              do i = 1,md
                 xs(nindex+i) = cfac*vs(nindex+i)/(ii*w)
              end do
           end if
           nindex = nindex+md
        end do

     end if
    
     ! set the initial guess for the solution    
     x0 = xs


    
     ! set the intial values for the vectors
     s0 = matmul(a,x0)
     if(ntb /= 0) then
        do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
        end do
        call getrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
        do ib = 1,ntb
           s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = & 
                vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
        end do
     end if
     if(ntb /= mtot) then
        nindex = 0
        do im = 1,mtot
           md = 2*ll(im)+1
           ltmp = check_mode(w,ww(im),wtb)
           if(.not.ltmp) then
              do i = 1,md
                 cfac = 1.0_sp/(ww(im)**2-w**2)
                 s0(nindex+i) = cfac*s0(nindex+i)
              end do
           end if
           nindex = nindex+md
        end do             
     end if
     s0 = xs-s0
     sb0 = s0
     p0  = s0
     pb0 = sb0
     xr  = x0
     

     

     ! check the initial error
     tmp1 = my_dot(sb0,s0)
     err = sqrt(abs(tmp1/my_dot(xs,xs)))
     if(err < tol) goto 1000

     
     count = 0   
     atimeit = 0 
     it = 0
     ! start the iteration loop 


     
     do it = 1,maxit

        
        call cpu_time(timeit1)

        ! iterate the vectors
        vt = matmul(a,p0)     
        if(ntb /= 0) then
           do ib = 1,ntb
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                   vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
           end do
           call getrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
           do ib = 1,ntb
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = & 
                   vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
           end do
        end if
        if(ntb /= mtot) then
           nindex = 0
           do im = 1,mtot
              md = 2*ll(im)+1
              if(.not.check_mode(w,ww(im),wtb)) then
                 do i = 1,md
                    cfac = 1.0_sp/(ww(im)**2-w**2)
                    vt(nindex+i) = cfac*vt(nindex+i)
                 end do
              end if
              nindex = nindex+md
           end do
        end if       
        alpha = tmp1/my_dot(pb0,vt)
        s1 = s0-alpha*vt
               
        if(ntb /= 0) then
           do ib = 1,ntb
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                   pb0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
           end do
           call getrs('T',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
           do ib = 1,ntb
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = & 
                   vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
           end do
        end if
        if(ntb /= mtot) then
           nindex = 0
           do im = 1,mtot
              md = 2*ll(im)+1
              if(.not.check_mode(w,ww(im),wtb)) then
                 do i = 1,md
                    cfac = 1.0_sp/(ww(im)**2-w**2)
                    vt(nindex+i) = cfac*pb0(nindex+i)
                 end do
              end if
              nindex = nindex+md
           end do
        end if
        vt = matmul(transpose(a),vt)
        
        sb1  = sb0-alpha*vt
        tmp2 = my_dot(sb1,s1)
        beta = tmp2/tmp1
        p1   = s1+beta*p0
        pb1  = sb1+beta*pb0

        ! update the solution
        xr = xr+alpha*p0
        
        ! estimate the error
        vt = matmul(a,xr)
        if(ntb /= 0) then
           do ib = 1,ntb
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                   vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
           end do
           call getrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
           do ib = 1,ntb
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = & 
                   vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
           end do
        end if
        if(ntb /= mtot) then
           nindex = 0
           do im = 1,mtot
              md = 2*ll(im)+1
              if(.not.check_mode(w,ww(im),wtb)) then
                 do i = 1,md
                    cfac = 1.0_sp/(ww(im)**2-w**2)
                    vt(nindex+i) = cfac*vt(nindex+i)
                 end do
              end if
              nindex = nindex+md
           end do
        end if
        count = it
        vt = vt-xs
        
        err = sqrt(abs(my_dot(vt,vt)/my_dot(xs,xs))) 
        if(err < tol) exit
                 
        ! update the vectors
        s0 = s1
        sb0 = sb1
        p0 = p1
        pb0 = pb1
        tmp1 = tmp2
        
!        if(it == maxit) stop 'no convergence'        
        
        if(it == maxit) goto 1000
        
        call cpu_time(timeit2)

        atimeit = (count*atimeit+(timeit2-timeit1))/(count+1)
     

   
     end do
     
     
1000 continue
     



     
     ! form the product with the receiver vector
     aclr = my_dot(vr,xr)   
     

          
     acl(iw) = aclr

     ! convert displacement to acceleration
     acl(iw) = -w*w*acl(iw)
     

     
     
  end do
  ! end loop over frequency


  
 
  




  
  ! write out list of modes and frequencies
  open(io1,file='modes.out')  
  do im = 1,mtot
     if(iq(im) == 't') iq(im) = 'T'
     if(iq(im) == 's') iq(im) = 'S'
     write(io1,3) nn(im),iq(im),ll(im),real(ww(im))
  end do
3 format(i3,1x,a1,1x,i3,f16.8)
  
  close(io1)

  ! parameters for frequency-domain filter
!  fac = 0.05_sp
  fac = 0.01_sp
  f11 = f1
  f22 = f2
  f12 = f11+fac*(f22-f11)
  f21 = f22-fac*(f22-f11)
  
  allocate(dat(1,nt))  
  dat(1,:) = acl(:)

  ! filter the spectra
  do i = 1,nt/2+1
     f = (i-1)*df     
     filt = hann(f,f11,f12,f21,f22)
     dat(1,i) = filt*dat(1,i)
  end do
    
  ! do the negative frequencies
  j = 0
  do i = nt/2+2,nt
     j = j+1
     dat(1,i) = conjg(dat(1,nt/2-j+1))
  end do  
  
  ! do the inverse Fourier transform
  call fourrow_dp(dat,1)
  dat = dat/(dt*nt)  
      
  ! undo the exponential decay on the time-series
  do i = 1,nt
     t = (i-1)*dt
     if(t > tout) exit
     dat(1,i) = dat(1,i)*exp(ep*t)
  end do

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
  
  ! filter the spectra again
  do i = 1,nt/2+1
     f = (i-1)*df
     filt = hann(f,f11,f12,f21,f22)
     dat(1,i) = filt*dat(1,i)
  end do

  ! converts displacements into cm (according to Arwen's program)
  acl(:)  = dat(1,:)!*0.95179e-30_sp
  
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

  function check_mode(w,wm,wtb)
    use nrtype
    implicit none
    logical(lgt) :: check_mode
    complex(spc) :: w,wm
    real(sp) :: wtb

!    if(abs(real(ww(im)-w)) <= wtb) then
    if(abs(real(wm-w)) <= wtb) then
       check_mode = .true.
    else
       check_mode = .false.
    end if

    return
  end function check_mode




end program dspec5
 
