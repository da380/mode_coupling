module module_util

  interface my_reallocate
     module procedure my_reallocate_iv,my_reallocate_ia,my_reallocate_rv, &
          my_reallocate_dv,my_reallocate_cv,my_reallocate_zv
  end interface

  interface poly_interp
     module procedure poly_interp_r, poly_interp_z
  end interface


  interface bilinear
     module procedure bilinear_d
  end interface

  
  interface get_float
     module procedure get_float_sp, get_float_dp
  end interface

  interface my_dot
     module procedure my_dot_c,my_dot_z
  end interface



contains
    

    subroutine my_reallocate_iv(a,n)
      use nrtype
      implicit none
      integer(i4b), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      integer(i4b), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_iv


    subroutine my_reallocate_ia(a,n,m)
      use nrtype
      implicit none
      integer(i4b), dimension(:,:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n,m
      integer(i4b), dimension(:,:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store,m_tmp,m_store
      n_tmp=size(a,1); m_tmp=size(a,2)
      n_store=min(n,n_tmp); m_store=min(m,m_tmp)      
      allocate(a_tmp(n_store,m_store))
      a_tmp(1:n_store,1:m_store)=a(1:n_store,1:m_store)
      deallocate(a); allocate(a(n,m))
      a(1:n_store,1:m_store)=a_tmp(1:n_store,1:m_store)
      deallocate(a_tmp)
      return
    end subroutine my_reallocate_ia



    subroutine my_reallocate_rv(a,n)
      use nrtype
      implicit none
      real(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_rv


    subroutine my_reallocate_dv(a,n)
      use nrtype
      implicit none
      real(dp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      real(dp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_dv
      

    subroutine my_reallocate_cv(a,n)
      use nrtype
      implicit none
      complex(sp), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(sp), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_cv


    subroutine my_reallocate_zv(a,n)
      use nrtype
      implicit none
      complex(dpc), dimension(:), allocatable, intent(inout) :: a
      integer(i4b), intent(in) :: n
      complex(dpc), dimension(:), allocatable :: a_tmp
      integer(i4b) :: n_tmp,n_store
      n_tmp=size(a)
      n_store=min(n,n_tmp)     
      allocate(a_tmp(n_store)); a_tmp(1:n_store)=a(1:n_store)
      deallocate(a); allocate(a(n))
      a(1:n_store)=a_tmp(1:n_store); deallocate(a_tmp)
      return
    end subroutine my_reallocate_zv


    subroutine poly_interp_r(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      real(dp), intent(in) :: f1,fp1,f2,fp2
      real(dp), intent(out) :: f,fp
      real(dp) :: c2,c3,delta,deltai,deltai2

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_r


    subroutine poly_interp_z(x1,f1,fp1,x2,f2,fp2,x,f,fp)
      use nrtype
      implicit none
      real(dp), intent(in) :: x,x1,x2
      complex(dpc), intent(in) :: f1,fp1,f2,fp2
      complex(dpc), intent(out) :: f,fp
      real(dp) :: delta,deltai,deltai2
      complex(dpc) :: c2,c3

      if(x == x1) then
        f=f1
        fp=fp1
        return
      else if(x == x2) then
         f=f2
         fp=fp2
         return
      end if

      delta=x2-x1

      if(delta == 0) stop 'poly_interp: delta = 0'
      
      deltai=1.0_dp/delta
      deltai2=deltai*deltai

      c2=3.0_dp*deltai2*(f2-f1)-deltai*(fp2+2.0_dp*fp1)
      c3=deltai2*(fp2+fp1)-2.0_dp*deltai*deltai2*(f2-f1)

      delta=x-x1
      deltai=delta*delta
      deltai2=deltai*delta

      f=f1+fp1*delta+c2*deltai+c3*deltai2
      fp=fp1+2.0_dp*c2*delta+3.0_dp*c3*deltai


      return
    end subroutine poly_interp_z


    function bilinear_d(x1,x2,y1,y2,f11,f12,f21,f22,x,y)
      use nrtype
      implicit none
      real(dp) :: bilinear_d
      real(dp), intent(in) :: x1,x2,y1,y2,f11,f12,f21,f22,x,y
      

      bilinear_d = (f11*(x2-x)*(y2-y) & 
                   +f12*(x2-x)*(y-y1) & 
                   +f21*(x-x1)*(y2-y) &
                   +f22*(x-x1)*(y-y1))/((x2-x1)*(y2-y1))         

      return
    end function bilinear_d


    subroutine get_string(tag,string)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      character(len=*), intent(out) :: string

      write(6,'("'//tag//'")',advance='no')
      read(5,*) string
      string=adjustl(string)
      return
    end subroutine get_string

    
    subroutine get_integer(tag,int)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      integer(i4b), intent(out) :: int
      character(len=256) :: string,form,sls
      integer(i4b) :: sl,ios
      
      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         write(sls,*) sl
         sls=adjustl(sls)
         sl=len_trim(sls)
         form='(i'//sls(1:sl)//')'
         read(string,form,iostat=ios) int
         if(ios == 0) exit
         print *, ' input must be an integer'
      end do
            
      return
    end subroutine get_integer


    subroutine get_float_dp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_dp

    subroutine get_float_sp(tag,fpn)
      use nrtype
      implicit none
      character(len=*), intent(in) :: tag
      real(sp), intent(out) :: fpn
      character(len=256) :: string,form,sls,ips
      integer(i4b) :: sl,ios,ip

      do 
         write(6,'("'//tag//'")',advance='no')
         read(5,*) string
         string=adjustl(string)
         sl=len_trim(string)
         ip=index(string,'.')
         if(ip == 0) ip=sl
         write(sls,*) sl
         write(ips,*) sl-ip
         sls=adjustl(sls)
         sl=len_trim(sls)
         ips=adjustl(ips)
         ip=len_trim(ips)
         form='(f'//sls(1:sl)//'.'//ips(1:ip)//')'
         read(string,form,iostat=ios) fpn
         if(ios == 0) exit
         print *, ' input must be a floating point number'
      end do


      return
    end subroutine get_float_sp



    function string_length(string,sign)
      use nrtype
      implicit none
      integer(i4b) :: string_length
      integer(i4b), intent(in) :: sign
      character (len=*), intent(in) :: string
      select case(sign)
      case(1)
         string_length=0
         do 
            if(string(string_length+sign:string_length+sign) == ' ') exit
            string_length=string_length+sign
         end do
      case(-1)
         string_length=len(string)
         do 
            if(string(string_length:string_length) == ' ') exit
            string_length=string_length+sign
         end do
     case default
        stop 'bad input of sign to string_length'
     end select
     return
   end function string_length


    subroutine bound(xx,x,i1,i2,err)
      ! Given a monotonically increasing array, xx, and a value
      ! x, returns the indices i1 and i2 such that
      ! xx(i1) <= x <= xx(i2). The routine checks that x is in 
      ! range, and if not returns err=.true.
      use nrtype
      implicit none
      real(dp), dimension(:), intent(in) :: xx
      real(dp), intent(in) :: x
      integer(i4b), intent(out) :: i1,i2
      logical(lgt), intent(out) :: err
      integer(i4b) :: i,n
      real(dp), parameter :: tol=1.0e-5_dp
      real(dp) :: rtol
      
      err=.false.
      n=size(xx)
      rtol=tol*xx(n)
      if(x > xx(n) .and. x < xx(1)) then
         err=.true.
         return
      end if

      if(abs(x-xx(n)) < rtol) then
         i2=n
         i1=n-1
         return
      end if
      if(abs(x-xx(1)) < rtol) then
         i2=2
         i1=1
         return
      end if
      do i=n-1,1,-1
         if(x <= xx(i+1) .and. x >= xx(i)) then
            i2=i+1
            i1=i
            return
         end if
      end do
      err=.true.
      return
    end subroutine bound



    function my_dot_c(x,y)
      use nrtype
      implicit none
      complex(spc) :: my_dot_c
     complex(spc), dimension(:), intent(in) :: x,y
      integer(i4b) :: i,n
      
      my_dot_c = 0.0_sp
     n = size(x,1)
      do i = 1,n
         my_dot_c = my_dot_c + x(i)*y(i)
      end do
      
      return
      
    end function my_dot_c


    function my_dot_z(x,y)
      use nrtype
      implicit none
      complex(dpc) :: my_dot_z
      complex(dpc), dimension(:), intent(in) :: x,y
      integer(i4b) :: i,n
      
      my_dot_z = 0.0_dp
      n = size(x,1)
      do i = 1,n
         my_dot_z = my_dot_z + x(i)*y(i)
      end do
     return
      
    end function my_dot_z



!    function my_dot(x,y)
!      use nrtype
!      implicit none
!      complex(dpc) :: my_dot
!      complex(dpc), dimension(:), intent(in) :: x,y
 !     integer(i4b) :: i,n
     
!      my_dot = 0.0_dp
!      n = size(x,1)
!      do i = 1,n
!         my_dot = my_dot + x(i)*y(i)
!      end do
!     return
      
!   end function my_dot



    subroutine string_cat_int(stri,i,stro)
      use nrtype
      implicit none
      
      character(len=*), intent(in) :: stri
      integer(i4b), intent(in) :: i
      character(len=*), intent(out) :: stro

      character(len=10) :: stmp
      integer(i4b) :: j

      write(stmp,'(i10)') i
      j = floor(log10(real(i)))
      stro = trim(stri)//stmp(10-j:10)
      
      return
    end subroutine string_cat_int


end module module_util
