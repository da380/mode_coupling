program dspec_pro

  use nrtype
  use module_util
  use module_fourier
  implicit none

  character(len=1), dimension(:), allocatable :: iq
  character(len=256), parameter :: spec_in = 'dspec_cal.out.'
  character(len=256), parameter :: job_pref = 'dspec_job.'
  character(len=256) :: spec_out
  character(len=256) :: spec_file
  character(len=256) :: time_out
  character(len=256) :: time_file
  character(len=10) :: string
  character(len=10), dimension(:), allocatable :: rec_string

  integer(i4b), parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer(i4b) :: mtot,nelem,im,i,j,info,istat, & 
       nindex,ifgot,nw,iw,ir1,nt,i1,i2,mex,qex,nt0,ne, & 
       ntb,md,iter,mindex,ib,jb,tbd,it,count,ntime,ijob,njob, & 
       i11,i22,nr,ir

  integer(i4b), dimension(:), allocatable :: nn,ll,ity

  integer(i4b), dimension(:), allocatable :: ipivot

  integer(i4b), dimension(:,:), allocatable :: tbl

  real(sp), dimension(9) :: source
  real(sp), dimension(:,:), allocatable :: station

  real(dp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
       t1,t2,samp

  complex(dpc), dimension(:), allocatable :: ww
  complex(dpc), dimension(:,:), allocatable :: acl
  complex(dpc), dimension(:,:), allocatable :: dat,dtmp
  
  
  real(sp), dimension(12) :: amp
  real(sp), dimension(10) :: aker
  real(sp), dimension(4)  :: ar1,ar2
  real(sp) :: r1,r2
  integer(i4b) :: nord,lord

  ! common block for the eigenfunctions parameters
  common/premdata/amp,aker,ar1,ar2,r1,r2,nord,it,lord




  ! get the total number of jobs
  call get_integer(' njob = ',njob)
  call get_string(' spectra file = ',spec_out)
  call get_string(' time series file = ',time_out)
  
  ! read in the first job to get parameters
  ijob = 1
  write(string,'(i10)') ijob
  j = floor(log10(real(ijob)))
  open(io1,file=trim(job_pref)//string(10-j:10), & 
       form='unformatted')

  
  read(io1) mtot,nt
  allocate(iq(mtot),nn(mtot),ll(mtot),ity(mtot))
  read(io1) iq(1:mtot),nn(1:mtot),ll(1:mtot),ity(1:mtot)
  read(io1) source(1:9)
  read(io1) nr
  allocate(station(5,nr))
  allocate(rec_string(nr))
  read(io1) station(1:5,1:nr),rec_string(1:nr)
  read(io1) f11,f22,dt,tout,df0,wtb,t1,t2,df,i11,i22,ep,nelem
  close(io1)

  ! get the minimun frequency
  f1 = f11
  i1 = i11
  f2 = f22
  i2 = i22

  ! allocate the data arrays
  allocate(acl(nt,nr))
  acl(:,:) = 0.0_dp



  ! read in the first block of the spectra  
  open(io1,file=trim(spec_in)//string(10-j:10), & 
       form='unformatted')
  read(io1) acl(i11:i22,1:nr)
  close(io1)


  ! loop over the rest of the jobs
  do ijob = 2,njob

     write(string,'(i10)') ijob
     j = floor(log10(real(ijob)))
     open(io1,file=trim(job_pref)//string(10-j:10), & 
          form='unformatted')          
     read(io1) mtot,nt
     read(io1) iq(1:mtot),nn(1:mtot),ll(1:mtot),ity(1:mtot)
     read(io1) source(1:9)
     read(io1) nr
     read(io1) station(1:5,1:nr),rec_string(1:nr)
     read(io1) f11,f22,dt,tout,df0,wtb,t1,t2,df,i11,i22,ep,nelem
     close(io1)


     ! get the maximum frequency
     if(ijob == njob) then
        f2 = f22
        i2 = i22
     end if

     open(io1,file=trim(spec_in)//string(10-j:10), & 
          form='unformatted')
     read(io1) acl(i11:i22,1:nr)
     close(io1)
    
  end do


  ! open the mode catalog
  call openfl(7,'/home/david/coupling_standalone/data/SPRM1.BIN', & 
       1 ,0,0,istat,4096)
  
  ! set up mode catalog
  call setprem(7,0)
  ! find the depth of the source
  call findd(source(3),ir1)

  allocate(ww(mtot))
  
  !  start the loop over the modes      
  nindex = 0
  do im = 1,mtot
     
     ! get the modes
     call getprem(nn(im),ity(im),ll(im),ir1,ifgot)
          
     ! modes complex frequency
     ww(im) = cmplx(amp(1),0.5*amp(1)*amp(2))   
  
     
  end do
  ! end loop over the modes





  ! write out list of modes and frequencies
  open(io1,file='modes.out')  
  do im = 1,mtot
     if(iq(im) == 't') iq(im) = 'T'
     if(iq(im) == 's') iq(im) = 'S'
!     write(io1,3) nn(im),iq(im),ll(im),real(ww(im))
     write(io1,*) real(ww(im))
  end do
3 format(i3,1x,a1,1x,i3,f16.8)  
  close(io1)



  ! parameters for frequency-domain filter
  fac = 0.1_sp
  f11 = f1
  f22 = f2
  f12 = f11+fac*(f22-f11)
  f21 = f22-fac*(f22-f11)
  
  ! put accelerations in a temporary array
  allocate(dat(nr,nt))  
  do ir = 1,nr
     dat(ir,:) = acl(:,ir)
  end do


 
  ! filter the spectra
  do i = 1,nt/2+1
     f = (i-1)*df     
     filt = hann(f,f11,f12,f21,f22)
     dat(:,i) = filt*dat(:,i)
  end do
    
  ! do the negative frequencies
  j = 0
  do i = nt/2+2,nt
     j = j+1
     dat(:,i) = conjg(dat(:,nt/2-j+1))
  end do  
  
  ! do the inverse Fourier transform
  call fourrow_dp(dat,1)
  dat = dat/(dt*nt)  
      
  ! undo the exponential decay on the time-series
  do i = 1,nt
     t = (i-1)*dt
     if(t > tout) exit
     dat(:,i) = dat(:,i)*exp(ep*t)
  end do


  
  ! parameters for time-domain filter
  fac = 0.5_sp
  t11 = t1
  t22 = t2
  t12 = t11+fac*(t22-t11)
  t21 = t22-fac*(t22-t11)

  ! begin loop over recievers
  do ir = 1,nr
  

     time_file = trim(time_out)//'.'//trim(rec_string(ir))

     ! write out the time series
     open(io1,file=trim(time_file)) 
     

     do i = 1,nt
        t = (i-1)*dt
        if(t > tout) exit
        filt = hann(t,t11,t12,t21,t22)
        write(io1,'(2e25.12e3)') t,real(dat(ir,i))
     end do
     close(io1)

     
  end do
  ! end loop over recievers

  
  
  ! filter the time series
  do i = 1,nt
     t = (i-1)*dt
     filt = hann(t,t11,t12,t21,t22)
     dat(:,i) = filt*dat(:,i)
  end do
    
  ! pad the time series if necessary
  nt0 = floor(1.0_dp/df0*dt)
  if(nt0 > nt) then
     ne = log(real(nt0))/log(2.0_dp)+1
     nt0 = 2**ne
     allocate(dtmp(nr,nt))
     dtmp = dat
     deallocate(dat)
     allocate(dat(nr,nt0))
     dat(:,1:nt) = dtmp(:,:)
     dat(:,nt+1:nt0) = 0.0_dp
     nt = nt0
     df = 1.0_dp/(nt*dt)
     i1 = max(floor(f1/df),1)
     i2 = floor(f2/df)+1
     deallocate(dtmp)
  end if


  print *, 'padded time =',nt*dt/3600.0_dp
  
  ! do the Fourier transform
  call fourrow_dp(dat,-1)
  dat = dat*dt


  
  ! filter the spectra again
!  do i = 1,nt/2+1
!     f = (i-1)*df
!     filt = hann(f,f11,f12,f21,f22)
!     dat(:,i) = filt*dat(:,i)
!  end do



  ! begin loop over recievers
  do ir = 1,nr
     


     spec_file = trim(spec_out)//'.'//trim(rec_string(ir))
     
     ! write out the spectra
     open(io1,file=trim(spec_file)) 
     do i = i1,i2
        f = (i-1)*df*1000.0_sp
        samp = abs(dat(ir,i))
        write(io1,'(4e25.12e3)') f,real(dat(ir,i)),imag(dat(ir,i)),samp
     end do

     close(io1)
     
  end do
  ! end loop over recievers

  


contains 

  
  function hann(t,t11,t12,t21,t22)
    use nrtype
    implicit none
    logical(lgt) :: ltmp
    real(dp) :: hann
    real(dp), intent(in) :: t,t11,t12,t21,t22

    ltmp = (t11 == 0.0_dp .and. t12 == 0.0_dp & 
         .and. t21 == 0.0_dp .and. t22 == 0.0_dp)
    if(.not.ltmp) then
       if(t < t11) then
          hann = 0.0_dp
       else if(t >= t11 .and. t < t12) then
          hann = pi_d*(t-t11)/(t12-t11)
          hann = 0.5_dp*(1.0_dp-cos(hann))
       else if(t >= t12 .and. t < t21) then
          hann = 1.0_dp
       else if(t >= t21 .and. t < t22) then
          hann = pi_d*(t22-t)/(t22-t21)
          hann = 0.5_dp*(1.0_dp-cos(hann))
       else if(t >= t22) then
          hann = 0.0_dp
       end if       
    end if
    return
  end function hann


end program dspec_pro
