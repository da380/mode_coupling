module module_fourier



  contains

    subroutine fourrow_dp(data,isign)
      use nrtype; use nrutil, only : assert,swap
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: isign
      integer(i4b) :: n,i,istep,j,m,mmax,n2
      real(dp) :: theta
      complex(dpc), dimension(size(data,1)) :: temp
      complex(dpc) :: w,wp
      complex(dpc) :: ws
      n=size(data,2)
      call assert(iand(n,n-1) == 0, 'n must b a power of 2 in fourrow_dp')
      n2=n/2
      j=n2
      do i=1,n-2
         if(j > i) call swap(data(:,j+1),data(:,i+1))
         m=n2
         do
            if(m < 2 .or. j < m) exit
            j=j-m
            m=m/2
         end do
         j=j+m
      end do
      mmax=1
      do 
         if(n <= mmax) exit
         istep=2*mmax
         theta=pi_d/(isign*mmax)
         wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
         w=cmplx(1.0_dp,0.0_dp,kind=dpc)
         do m=1,mmax
            ws=w
            do i=m,n,istep
               j=i+mmax
               temp=ws*data(:,j)
               data(:,j)=data(:,i)-temp
               data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
         end do
         mmax=istep
      end do
    end subroutine fourrow_dp



    subroutine high_pass(data,i1,i2)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: i1,i2
      integer(i4b) :: i
      real(dp) :: r1

      if(i1 >= i2) stop 'bad input to low_pass'

      do i=1,min(i2,size(data,2)/2+1)
         if(i < i1) then
            data(:,i)=0.0_dp
         else
            r1=pi_d*real(i-i1)/real(i2-i1)
            r1=0.5_dp*(1.0_dp-cos(r1))
            data(:,i)=r1*data(:,i)
         end if
      end do

      return
    end subroutine high_pass



    subroutine low_pass(data,i1,i2)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: i1,i2
      integer(i4b) :: i
      real(dp) :: r1

      if(i1 > i2) stop 'bad input to low_pass'

      do i=i1,size(data,2)/2+1
         if(i > i2) then
            data(:,i)= 0.0_dp
         else
            r1=pi_d*real(i2-i)/real(i2-i1)
            r1=0.5_dp*(1.0_dp-cos(r1))
            data(:,i)=r1*data(:,i)
         end if
      end do


         return
    end subroutine low_pass

    subroutine time_chop(data,dt,ts,t1,t2,t3,t4)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt,ts,t1,t2,t3,t4
      integer(i4b) :: i,n,i1,i2
      real(dp) :: arg,t

      n=size(data,2)

      do i=1,n
         t=ts+(i-1)*dt
         if(t < t1) then
            data(:,i)=0.0_dp
         else if(t >= t1 .and. t < t2) then
            arg=pi_d*(t-t1)/(t2-t1)
            arg=0.5_dp*(1.0_dp-cos(arg))
            data(:,i)=arg*data(:,i)
         else if(t > t3 .and. t <= t4) then
            arg=pi_d*(t4-t)/(t4-t3)
            arg=0.5_dp*(1.0_dp-cos(arg))
            data(:,i)=arg*data(:,i)
         else if(t > t4) then
            data(:,i)=0.0_dp
         end if
         
      end do
      
      return
    end subroutine time_chop


    subroutine convolve(data,dt,tau)
      use nrtype; use module_util
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), intent(in) :: tau
      real(dp) :: w,dw,sinc

      n=size(data,2)
      dw=1.0_dp/(dt*n)

      do i=2,n/2+1
         w=twopi_d*(i-1)*dw
         sinc=sin(w*tau)/(w*tau)
         data(:,i)=data(:,i)*sinc
      end do
         
      return
    end subroutine convolve


  subroutine rick_convolve(data,dt,tau)
      use nrtype; use module_util
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), intent(in) :: tau
      real(dp) :: w,dw,sinc,t
      complex(dpc), dimension(:,:), allocatable :: ra

      n=size(data,2)
      dw=1.0_dp/(dt*n)
      allocate(ra(1,n))
      do i=1,n
         t=(i-1)*dt
         ra(1,i)=rick(t,tau)
      end do
      call fourrow_dp(ra,1)

      do i=2,n/2+1
         data(:,i)=data(:,i)*ra(1,i)
      end do
         
      return
    end subroutine rick_convolve


    subroutine accel(data,dt)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), parameter :: tau=10.0_dp
      real(dp) :: w,dw,sinc

      n=size(data,2)
      dw=1.0_dp/(dt*n)

      do i=2,n/2+1
         w=twopi_d*(i-1)*dw
         data(:,i)=-w**2*data(:,i)
      end do
         
      return
    end subroutine accel


    subroutine filter(data,dt,ts,w1,w2,w3,w4,acc,box,ric)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt,ts,w1,w2,w3,w4
      logical(lgt), intent(in), optional :: acc
      real(dp), intent(in), optional :: box,ric
      integer(i4b) :: n,i1,i2,i3,i4,i
      real(dp) :: dw,w,tau,t1,t2,t3,t4
      real(dp), parameter :: frac=0.8
      complex(dpc) :: filt


      n=size(data,2)
      dw=1.0_dp/(n*dt)
      i1=floor(w1/dw)+1
      i2=floor(w2/dw)+1
      i3=floor(w3/dw)+1
      i4=floor(w4/dw)+1
      t4=ts+dt*(size(data,2)-1)
      t3=frac*t4
      t1=ts
      t2=(1.0_dp-frac)*t4
      call time_chop(data,dt,ts,t1,t2,t3,t4)
      call fourrow_dp(data,1)
      if(i1 /= i2) then
         call high_pass(data,i1,i2)
      end if
      if(i3 /= i4) then
         call low_pass(data,i3,i4)
      end if
      if(acc) then
         call accel(data,dt)
      end if
      if(present(box)) then
         if(box > 0.0_dp) then
            tau=box
            call convolve(data,dt,tau)
         end if
      end if
      if(present(ric)) then
         if(ric > 0.0_dp) then
            tau=ric
            call rick_convolve(data,dt,tau)
         end if
      end if
      data(:,n/2+2:n)=conjg(data(:,n/2:2:-1))
      call fourrow_dp(data,-1)
      data=data/n

      return
    end subroutine filter



    subroutine get_power_spectra(data,dt,power,freq)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(in) :: data
      real(dp), intent(in) :: dt
      real(dp), dimension(:), allocatable, intent(out) :: freq
      real(dp), dimension(:,:), allocatable, intent(out) :: power
      integer(i4b) :: i,n
      real(dp) :: nfreq
      complex(dpc), dimension(:,:), allocatable :: data_tmp
  
      n=size(data,2)
      allocate(data_tmp(size(data,1),n))
      data_tmp=data
      call fourrow_dp(data_tmp,1)
      allocate(power(size(data,1),n/2+1),freq(n/2+1))
      power(:,1:n/2+1)=data_tmp(:,1:n/2+1)*conjg(data_tmp(:,1:n/2+1))
      do i=1,n/2+1
         freq(i)=(i-1)/(n*dt)
      end do
      deallocate(data_tmp)
      return
    end subroutine get_power_spectra

    
    function  rick(t,tau)
      use nrtype
      implicit none
      real(dp) :: rick
      real(dp), intent(in)  :: t,tau
      rick=-(twopi_d/tau)**2*(t-0.5_dp*tau) & 
           *exp(-2.0*pi_d**2*(t/tau-0.5_dp)**2)
      return
    end function rick
    
    



end module module_fourier
