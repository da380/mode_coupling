program QnmS
!program modified from dspec_pre

  use nrtype
  use nrutil
  use module_util
  use module_fourier
  !use SUNPERF
  implicit none

  logical(lgt) :: ltmp,exists

  ! dimension information in header file
  include 'dimension.h'

  integer(i4b), dimension(:,:), allocatable :: indj
 
  integer(i4b), parameter :: maxit = 100

  character(len=256) :: matrix_file
  
  character(len=1), dimension(:), allocatable :: iq
  character(len=256), parameter :: job_pref = 'dspec_job.'
  character(len=256) :: job_out,outdir,sysstring
  character(len=10) :: string
  character(len=256) :: source_file,cmtsolution,tempstr
  character(len=256) :: rec_file
  character(len=80) :: getunx
  character(len=10), dimension(:), allocatable :: rec_string



  integer(i4b), parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer(i4b) :: mtot,nelem,im,i,j,info,istat, & 
       nindex,ifgot,nw,iw,ir1,nt,i1,i2,mex,qex,nt0,ne, & 
       ntb,md,iter,mindex,ib,jb,tbd,it,count,ntime, & 
       njob,ncal,ntmp,i11,i22,nr,ir,nbyts,loutdir,ierror,lmat

  integer(i4b), dimension(:), allocatable :: nn,ll,ity

  real(sp), dimension(9) :: source
  real(sp), dimension(:,:), allocatable :: station

  real(dp) :: wr1,wr2,wr,wi,dwr,f1,f2,df,dt,tout,ep, & 
       f11,f12,f21,f22,t,t11,t12,t21,t22,fac,filt, & 
       f,df0,wtb,qq,time1,time2,atime,acount,fmin,fmax, &
       t1,t2
  

  complex(spc) :: w,aclr,aclr0,cfac
  complex(spc) :: alpha,beta,tmp1,tmp2
  complex(spc), dimension(:), allocatable :: ww

  real(sp), dimension(12) :: amp
  real(sp), dimension(10) :: aker
  real(sp), dimension(4)  :: ar1,ar2
  real(sp) :: r1,r2
  integer(i4b) :: nord,lord
  
  ! cmt source
  integer(i4b) :: iyr,imo,ida,iho,imi,jda,julday,ityphdur,itypso,ierr
  real(sp) :: hdur,epla,eplo,dep,fsec,tcmt
  real(sp), dimension(6) :: xm

  ! common block for the eigenfunctions parameters
  common/premdata/amp,aker,ar1,ar2,r1,r2,nord,it,lord
  
  call chekcl('|   -f1:o:1:minimum frequency for the output spectra (in mHz)[0.1] '  & 
       //'|   -f2:o:1:maximum frequency for the output spectra (in mHz) [2]'  & 
       //'|   -dt:o:1:time step for the time-series (in seconds) [20]'  & 
       //'|   -tout:o:1:length of time series (in hours) [128]'  & 
       //'|   -df0:o:1:fractional width of a filter in the frequency domain mhz [0.01]'  & 
       //'|   -wtb:o:1:a parameter used in the IDSM [0.05]'  & 
       //'|   -t1:o:1:start time for a cosine bell window (in hours) [0.]'  & 
       //'|   -t2:o:1:end time for a cosine bell window (in hours) [tout]'  & 
       //'|   -cm:o:1:coupling matrix [matrix_parts.bin]'  & 
       //'|   -C:o:1:CMT source file [CMTSOLUTION]'  & 
       //'|   -L:o:1:List with records [DATA/RECORDHEADERS]'     &  
       //'|   -D:o:1:Directory for output files [SYNT]'    &
       //'|   -s:o:1:Scratch directory for temporary files [TEMP]'    &
       //'|   -njob:o:1:number of processors [4]'  &                   
       //'|   -v:o:0,1:Verbosity level [0]'  &                   
       //'|')

  tempstr = getunx('-njob',1,nbyts)
  read (tempstr,*) njob
  tempstr = getunx('-f1',1,nbyts)
  read (tempstr,*) f1
  tempstr = getunx('-f2',1,nbyts)
  read (tempstr,*) f2
  tempstr = getunx('-dt',1,nbyts)
  read (tempstr,*) dt
  tempstr = getunx('-tout',1,nbyts)
  read (tempstr,*) tout
  tempstr = getunx('-df0',1,nbyts)
  read (tempstr,*) df0
  tempstr = getunx('-wtb',1,nbyts)
  read (tempstr,*) wtb
  tempstr = getunx('-t1',1,nbyts)
  read (tempstr,*) t1  
  tempstr = getunx('-t2',1,nbyts)
  if (tempstr(1:lnblnk(tempstr)) == 'tout') then
    t2=tout
  else
    read (tempstr,*) t2
  endif

  matrix_file=getunx('-cm',1,lmat)
  inquire(file=matrix_file,exist=exists)
    if(.not.exists) then
          write(6,"('input file does not exist:',a)") matrix_file(1:lmat)
          call exit(1)
    endif
    

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


  ! read in source parameters
  !call get_string(' source file = ',source_file )
  !open(99,file=trim(source_file))
  !do i = 1,9
  !   read(99,*) source(i)
  !end do
  !close(99)

!---- get the filename for moment tensor
!
  cmtsolution=getunx('-C',1,nbyts)
  inquire(file=cmtsolution,exist=exists)
  if(.not.exists) then
    write(6,"('moment tensor file does not exist:',a)") cmtsolution(1:nbyts)
    call exit(1)
  endif!

!---- read the CMT information
!
  open(1,file=cmtsolution)
  call rdsource(1,iyr,imo,ida,iho,imi,fsec,tcmt, & 
          ityphdur,hdur,epla,eplo,dep,itypso,xm,ierr)
  jda=julday(iyr,imo,ida)
  write(6,"('time shift:',f12.4)") tcmt
  if(ityphdur.eq.1) then
    write(6,"('boxcar time function')")
  else
    write(6,"('triangular time function')")
  endif
  write(6,"('half duration:',f9.4)") hdur
  write(6,"('latitude:  ',f12.4)") epla
  write(6,"('longitude: ',f12.4)") eplo
  write(6,"('depth:     ',f12.4)") dep
  if(itypso.eq.0) then
        write(6,"('Mrr:    ',e15.6)") xm(1)
        write(6,"('Mtt:    ',e15.6)") xm(2)
        write(6,"('Mpp:    ',e15.6)") xm(3)
        write(6,"('Mrt:    ',e15.6)") xm(4)
        write(6,"('Mrp:    ',e15.6)") xm(5)
        write(6,"('Mtp:    ',e15.6)") xm(6)
  else if(itypso.eq.1) then
        write(6,"('CFr:    ',e15.6)") xm(1)
        write(6,"('CFt:    ',e15.6)") xm(2)
        write(6,"('CFp:    ',e15.6)") xm(3)
        write(6,"('... CF moment tensor not tested yet')") 
        call exit(1)
  endif
  close(1)

  ! move values to source for rest of the code
  source(1)=epla
  source(2)=eplo
  source(3)=dep
  do j=4,9
     source(j)=xm(j-3)
  end do
  ! normalize the moment tensor so that the 
  ! output displacements are in cm's
  source(4:9) = source(4:9)*0.95179e-30_sp


  ! receiver location
  !call get_string(' receiver file = ',rec_file )
  rec_file=getunx('-L',1,nbyts)
  inquire(file=rec_file,exist=exists)
  if(.not.exists) then
    write(6,"('receiver file does not exist:',a)") rec_file(1:nbyts)
    call exit(1)
  endif!
  
  open(99,file=trim(rec_file))
  read(99,*) nr
  allocate(station(5,nr))
  allocate(rec_string(nr))
  do j = 1,nr
        read(99,*) station(1:5,j),rec_string(j)
  end do
  close(99)

!
!---- get the name of the output directory
!
  outdir=getunx('-D',1,loutdir)
  if(loutdir.gt.0) then
    inquire(file=outdir,exist=exists)
    if(.not.exists) then
          sysstring='mkdir -p '//outdir(1:loutdir)
          ierror=system(sysstring)
          if(ierror.ne.0) then
            write(6,"('unable to create output directory',a)") outdir(1:loutdir)
            call exit(1)
          endif
    endif
  endif
  
  ! get parameters for the calculations
!   call get_float(' f1 (mhz) = ',f1)
!   call get_float(' f2 (mhz) = ',f2)            
!   call get_float(' dt (sec) = ',dt)
!   call get_float(' tout (hrs) = ',tout)
!   call get_float(' df0 (mhz) = ',df0)
!   call get_float(' wtb (mhz) = ',wtb)
!   call get_float(' t1 (hrs) = ',t1)
!   call get_float(' t2 (hrs) = ',t2)
!   call get_integer(' number of jobs = ',njob)
!   call get_string( ' matrix file = ',matrix_file)
  
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
  qex = 4
  call fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)


  print *, 'tout = ',tout/3600.0_dp
  print *, 'dt = ',dt
  print *, 'nt = ',nt
  print *, 'nt*dt = ',nt*dt/3600.0_dp


  ! divide the frequencies between the processors
  allocate(indj(njob,2))


  call jobs(i1,i2,njob,indj)

  
  ! number of frequencies to do 
  ncal = i2-i1+1



  do i = 1,njob

     ! open the binary job file     
     write(string,'(i10)') i
     j = floor(log10(real(i)))
     open(io1,file=trim(job_pref)//string(10-j:10), & 
          form='unformatted')


     
     i11 = indj(i,1)
     i22 = indj(i,2)



     ! write data to the job file
     write(io1) mtot,nt
     write(io1) iq(1:mtot),nn(1:mtot),ll(1:mtot),ity(1:mtot)
     write(io1) source(1:9)
     write(io1) nr
     write(io1) station(1:5,1:nr),rec_string(1:nr)
     write(io1) f1,f2,dt,tout,df0,wtb,t1,t2,df,i11,i22,ep,nelem
     write(io1) matrix_file
     close(io1)

     

  end do

  

  

contains


  subroutine fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    use nrtype
    implicit none
    real(dp), intent(inout) :: f1,f2
    real(dp), intent(in) :: dt,tout
    integer(i4b), intent(in) :: mex,qex
    real(dp), intent(out) :: df,ep
    integer(i4b), intent(out) :: nt,i1,i2
    
    integer(i4b) :: ne
    real(dp) :: fn
    
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
    
    i1 = max(floor(f1/df),2)
    f1 = (i1-1)*df          
    
    i2 = floor(f2/df)+2
    f2 = (i2-1)*df

    return
  end subroutine fcal


  subroutine jobs(i1,i2,njob,indj)
    use nrtype
    use nrutil
    implicit none
    integer(i4b), intent(in) :: i1,i2
    integer(i4b), intent(in) :: njob
    integer(i4b), dimension(njob,2) :: indj


    integer(i4b) :: count,ijob,ncal,mn,mx,imn,imx
    integer(i4b), dimension(njob) :: ncalj
    real(dp) :: r,f,fac


    fac = 0.9_dp
    ncal = i2-i1+1

    ncalj = 0
    count = 0
    ijob = 1
    do      
       if(count == ncal) exit
       call random_number(r)
       f = fac*real(ijob)/real(njob+1)
       if(r >= f) then
          ncalj(ijob) = ncalj(ijob)+1
          count = count+1
       end if
       if(ijob < njob) then
          ijob = ijob+1
       else
          ijob = 1
       end if
    end do

    ! check for null-jobs
    do 
       
       imn = iminloc(ncalj)
       mn = ncalj(imn)

       if(mn > 1) exit

       imx = imaxloc(ncalj)
       mx = ncalj(imx)

       if(mx == 3) stop ' too many jobs'

       ncalj(imn) = 2
       ncalj(imx) = ncalj(imx)-2

    end do


    
    indj(1,1) = i1
    do ijob = 1,njob
       indj(ijob,2) = indj(ijob,1)+ncalj(ijob)-1
       if(ijob < njob) then
          indj(ijob+1,1) = indj(ijob,2)+1
       end if
    end do


    
    
    
    
    return
  end subroutine jobs




end program QnmS
