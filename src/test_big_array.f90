      program test_big_array
        integer,parameter::NDIM=10000
      complex*8,dimension(NDIM,NDIM)::big
      
      fac=1.3
      big=3.

      big=big/fac

      print *,"Done"
      end program
