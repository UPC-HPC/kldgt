      subroutine agc_wrapper(pinput, pscalar, ntrc, nsamp, gatelen)

      include 'msdginterp.fi'
 
      call bsdg_agc_gxcon(ntrc,nsamp,pinput,pscalar,gatelen)

      end 
