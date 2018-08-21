      program test_big_array
      parameter NDIM=10000
      complex*8 big(NDIM,NDIM)
      
      fac=1.3
      big=3.

      big=big/fac

      print *,"Done"
      end program
