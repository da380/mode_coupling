program dspec_cal

  use nrtype
  use module_util
  use module_fourier 
  implicit none

  logical(lgt) :: ltmp

  ! dimension information in header file
  include 'dimension.h'

  integer(i4b), parameter :: maxit = 10000

  character(len=256) :: datapath  
  
  character(len=1), dimension(:), allocatable :: iq
  character(len=256), parameter :: spec_out = 'dspec_cal.out.'
  character(len=10) :: string
  character(len=256), parameter :: job_pref = 'dspec_job.'

  character(len=256) :: matrix_file,file
  character(len=80) :: getunx,string2
  character(len=10), dimension(:), allocatable :: rec_string

  integer(i4b), parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer(i4b) :: mtot,nelem,im,i,j,info,istat, & 
       nindex,ifgot,nw,iw,ir1,nt,i1,i2,mex,qex,nt0,ne, & 
       ntb,md,iter,mindex,ib,jb,tbd,count,ntime,ijob, &
       ifcor,im2,nr,ir,nbyts,lmod,lpath

  integer(i4b), dimension(:), allocatable :: nn,ll,ity

  integer(i4b), dimension(:), allocatable :: ipivot

  integer(i4b), dimension(:,:), allocatable :: tbl

  real(sp), dimension(9) :: source
  real(sp), dimension(:,:), allocatable :: station

!  real(sp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
!       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
!       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
!       t1,t2

  real(dp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
       t1,t2
  
!  real(sp) :: err
  real(dp) :: err



  real(dp), parameter :: tol = 1.e-6_dp


  complex(dpc) :: w,aclr,aclr0,cfac,w0
  complex(dpc) :: alpha,beta,tmp1,tmp2
  complex(dpc), dimension(:), allocatable :: ww
  complex(dpc), dimension(:,:), allocatable :: a0,a1,a2

  complex(dpc), dimension(:,:), allocatable :: a,atb
  complex(dpc), dimension(:,:), allocatable :: vsd
  complex(dpc), dimension(:), allocatable :: xs,xr,x0


  complex(dpc), dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, & 
                                              sb1,p1,pb1,vt

  integer(i4b) :: saved 
  complex(dpc), dimension(:), allocatable :: xsave
  
  complex(dpc), dimension(:), allocatable :: vs
  complex(dpc), dimension(:,:), allocatable :: vr
  complex(spc), dimension(2*ML+1) :: vrt,vst

  complex(dpc), dimension(:,:), allocatable :: acl

  real(sp), dimension(12) :: amp
  real(sp), dimension(10) :: aker
  real(sp), dimension(4)  :: ar1,ar2
  real(sp) :: r1,r2
  integer(i4b) :: nord,it,lord

  ! common block for the eigenfunctions parameters
  common/premdata/amp,aker,ar1,ar2,r1,r2,nord,it,lord

  ! get the job id number
  !call get_command_argument(1,string)
  !read(string,*) ijob
  call chekcl('|-p:o:1:[../data]'                   & 
       //'|-j:r:1:job ID'                 &  
       //'|-m:o:1:[SPRM1.BIN] mode catalog'  &                   
       //'|')

      string2 = getunx('-j',1,nbyts)
      read(string2,*) ijob

  ! open the job input file
  write(string,'(i10)') ijob
  j = floor(log10(real(ijob)))
  open(io1,file=trim(job_pref)//string(10-j:10), & 
       form='unformatted')



  ! read data to the job file
  read(io1) mtot,nt
  allocate(iq(mtot),nn(mtot),ll(mtot),ity(mtot))
  read(io1) iq(1:mtot),nn(1:mtot),ll(1:mtot),ity(1:mtot)
  read(io1) source(1:9)
  read(io1) nr
  allocate(station(5,nr))
  allocate(rec_string(nr))
  read(io1) station(1:5,1:nr),rec_string(1:nr)
  read(io1) f1,f2,dt,tout,df0,wtb,t1,t2,df,i1,i2,ep,nelem
  read(io1) matrix_file
  close(io1)


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
  allocate(vr(nelem,nr),vs(nelem))
 
  ! read in the coupling matrices 
  open(io1,file=trim(matrix_file),form='unformatted')
  read(io1) a0
  read(io1) a1
  read(io1) a2
  close(io1)

! location of files
  datapath=getunx('-p',1,lpath)
  
  ! open the mode catalog
  file = datapath(1:lpath)//'/'//getunx('-m',1,lmod)

  call openfl(7,file, & 
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
     ww(im) = dcmplx(amp(1),0.5*amp(1)*amp(2))   

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


  f1 = (i1-1)*df
  f2 = (i2-1)*df
  

  wr1 = twopi_d*f1
  wr2 = twopi_d*f2
  dwr = twopi_d*df
  wi = ep


  nt = i2-i1+1

  ! initialize the output arrays
  allocate(acl(nt,nr))  
  acl(:,:)  = 0.0_sp

  ntime = 0
  atime = 0
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
                 a1(nindex+i,mindex+j) = 0.0_dp
                 a2(nindex+i,mindex+j) = 0.0_dp
              end do
              mindex = mindex + 2*ll(im2)+1
           end do
        end do
        nindex = nindex + 2*ll(im)+1
     end do
  end if
  
  


  ! begin loop over frequency
  do iw = 1,nt
     
     call cpu_time(time1)

     wr = wr1+(iw-1)*dwr
     w = wr-ii*wi
     
     print *, ' frequencies to go = ',nt-iw+1

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
     do im = 1,mtot
        if(check_mode(w,ww(im),wtb)) then
           ntb = ntb+1
           tbl(ntb,1) = im         
           tbl(ntb,2) = nindex     
           tbl(ntb,3) = mindex     
           tbl(ntb,4) = 2*ll(im)+1 
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
        call zgetrf(tbd,tbd,atb,nelem,ipivot,info)
        
        ! assemble the source vector
        do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) & 
                = vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))/(ii*w)
        end do

        ! solve the linear system for target block 
        call zgetrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
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
              do i = 1,md
                 cfac = 1.0_sp/(ww(im)**2-w**2)
                 xs(nindex+i) = cfac*vs(nindex+i)/(ii*w)
              end do
           end if
           nindex = nindex+md
        end do

     end if
    
     ! set the initial guess for the solution    
!     if(saved == 0) then       
        x0 = xs
!     else
!        x0 = xsave
!     end if          
    
     ! set the intial values for the vectors
     s0 = matmul(a,x0)
     if(ntb /= 0) then
        do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
        end do
        call zgetrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
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
     if(err < tol) then
        print *, 'no iterations',' err = ',err/tol
        goto 1000
     end if
     
     count = 0   
     ! start the iteration loop 
     
     do it = 1,maxit
        

        ! iterate the vectors
        vt = matmul(a,p0)     
        if(ntb /= 0) then
           do ib = 1,ntb
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = & 
                   vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) 
           end do
           call zgetrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
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
           call zgetrs('T',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
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
           call zgetrs('N',tbd,1,atb,nelem,ipivot,vsd,nelem,info)
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
        
        print *, 'iteration = ',it,' err = ',err/tol

        if(err < tol) exit

                 
        ! update the vectors
        s0 = s1
        sb0 = sb1
        p0 = p1
        pb0 = pb1
        tmp1 = tmp2
        
        if(it == maxit) stop 'no convergence'        
!        if(it == maxit) goto 1000
        
     end do
      
     
1000 continue
    
     ! begin loop over the recievers
     do ir = 1,nr
        ! form the product with the receiver vector
        aclr = my_dot(vr(:,ir),xr)   
     
        acl(iw,ir) = aclr

     end do
        
     ! output is the acceleration spectrum
     acl(iw,:) = -w*w*acl(iw,:)

     ! output is displacement spectrum for Dirac source
!     acl(iw,:) = ii*w*acl(iw,:)
     
     xsave = xr
     saved = 1
     
     call cpu_time(time2)
     atime = (ntime*atime+(time2-time1))/(ntime+1)
     acount = (ntime*acount+count)/(ntime+1)
     ntime = ntime+1
     
  end do
  ! end loop over frequency
  
  print *, '========================='
  print *, 'atime = ',atime
  print *, 'acount = ',acount


  ! write out the results as a binary file

  write(string,'(i10)') ijob
  j = floor(log10(real(ijob)))
  open(io1,file=trim(spec_out)//string(10-j:10), & 
       form='unformatted')
  write(io1) acl

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



  function check_mode(w,wm,wtb)
    use nrtype
    implicit none
    logical(lgt) :: check_mode
    complex(dpc) :: w,wm
    real(dp) :: wtb

!    if(abs(real(ww(im)-w)) <= wtb) then
    if(abs(real(wm-w)) <= wtb) then
       check_mode = .true.
    else
       check_mode = .false.
    end if

    return
  end function check_mode


end program dspec_cal
