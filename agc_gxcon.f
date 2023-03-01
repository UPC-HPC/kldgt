c<<
c NAME
c     agc_gxcon - GX phase of agc with constant size gates.
c
c SYNOPSIS
      subroutine bsdg_agc_gxcon(count,length,pinput,pscalar,gate)
      include 'msdginterp.fi'

      integer*8 base
      integer count
      integer   status
      integer i,j, k, ii
      integer gate, applypt, scgate, scapplypt
      integer trid, stride
      integer skip, length, domain
      real    normvalue, threshold, deadgate
      real    pwr_sngl
      real*8  pwr_dbl
      logical l_double

      integer ist,lst
      logical l_save, threshopt, l_wholetr
      integer sclamp

      integer*8 ix, ix_scr_all, ix_scr2, ix_wtr
      integer*8 ix_scr, ln_scr
      integer*8 ix_lim, ix_thrlim
      integer ln_wtr
      integer*8 ix_rmed, ix_imed
      integer ln_med
      integer*8 ix_wrk
      integer*8 ix_gainlimit, cur_scr
      integer ntime, cur_gate, cur_applypt
      real    threshval,temp
      real    weight1,weight2,scalarf
      real    gainlimit

      integer*8 mftime, mflength, mflist
      integer mfnpairs
      real median
      integer*8 tempwork,tempscl
      integer indexbegin, indexend

      real    trimalpha, trimgatelen
      real    meanval
      integer*8 meanlist, trimlist,ixbase
      integer ibeg, iend,nin

      real pinput(*)
      real ptrace(length)
      real pscalar(*)
      real pscalartr(length)
c.......................................................................

c.....Get and zero scratch space for the scaletrace

      l_double = .true.
      normvalue = 1.0
      applypt = gate/2
c.......................................................................
c.....For each trace.
c.......................................................................

      do i = 0,count-1
         
         do j = 1,length
            pscalartr(j) = 0.0
            ptrace(j) = pinput(i*length+j)
         end do

c........Find first and last non-zero sample.

         lst = ML_imut(length,ptrace,-1)
         if (lst .ne. 0 ) then 
            ist = ML_imut(length,ptrace,1)
         endif

c........Calculate scalers for live input samples.

         if (lst.ne. 0) then
            
            cur_gate = gate
            cur_applypt = applypt
                        
            call bsdg_agc_con_scal_amp(length, ist, lst, cur_gate, cur_applypt,
     *           normvalue, l_double, ptrace, pscalartr )
            
         endif

         do j = 1,length
            pscalar(i*length+j)=pscalartr(j)
         end do

      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c<<
c NAME
c     agc_con_scal_amp - agc scalar with constant size gates var precision.
c
c SYNOPSIS
      subroutine bsdg_agc_con_scal_amp(length,ist,lst,gate,applypt,normvalue,
     *   l_dbl,vecin,sclout)
      include 'msdginterp.fi'    

c
c ARTGUMENTS
c     myid - Id of my structure.
c
c     length   I   vector length
c     ist      I   first live sample limit 1
c     lst      I   last live sample  limit length
c     gate     I   gate length in samples
c     applypt  I 
c     l_dbl    I logical true for double precision 
c     vecin    I   input vector
c     sclout   O   work aray length

c
      integer length
      integer ist
      integer lst
      integer gate
      integer applypt
      real    normvalue
      logical l_dbl
      real    vecin(1:*)
      real    sclout(1:*)
c DESCRIPTION
c     The scalar calculation used in agc algorithms.
c
c REVISION
c 
c 16-Jul-1998	APD	v20	Gnats PR 
c                               make this a subroutine for use in other places
c
c>>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer k
      real    scalval, deadgate
      real    rin, rout
      real    pwr_sngl
      real*8  pwr_dbl
      integer nnrm 
      integer ibeg, iend
      integer iobeg, ioend, ioapp

      real AGC_EPS

      deadgate = 0.0

      AGC_EPS = 1.e-15

      if ( l_dbl ) then 
c
c     Double precision loop for amp scaling   
c     
c Form first scaler sum.
            pwr_dbl = 0.0
            nnrm = 0

c Calculate pwr from live to live+gate.
            ibeg = 1 + ist -1
            iend = 1 + min (ist + gate -1, lst-1)
            do k=ibeg,iend
               if(vecin(k).gt.AGC_EPS .or. vecin(k).lt.-AGC_EPS) then
                 pwr_dbl = pwr_dbl+abs(vecin(k))
                 nnrm = nnrm +1 
               endif
            enddo

c Store first scaler in scaletrace, from begining of gate to applypt.
            iobeg = 1+ist -1
            ioend = 1+min(ist-1+applypt,lst-1)
            if (pwr_dbl .gt. AGC_EPS) then
               scalval = nnrm*normvalue/pwr_dbl
               do k=iobeg,ioend
                  sclout(k) = scalval
               enddo 
            else
               do k=iobeg,ioend
                  sclout(k) = deadgate
               enddo
            endif

c Now do middle portion of trace.
            ibeg =  1 + ist -1
            iend =  1 + lst-gate-2
            ioapp = ioend +1
            do k= ibeg,iend
               rout = vecin(k)
               rin  = vecin(k +gate +1)
               if(rout.gt.AGC_EPS .or. rout.lt.-AGC_EPS) nnrm = nnrm - 1
               if(rin .gt.AGC_EPS .or. rin .lt.-AGC_EPS) nnrm = nnrm + 1
               pwr_dbl = pwr_dbl-abs(rout)+abs(rin)
               if (pwr_dbl .gt. AGC_EPS) then
                  sclout(ioapp) = nnrm*normvalue/pwr_dbl
               else
                  sclout(ioapp) = deadgate
               endif
               ioapp = ioapp + 1
            enddo

c Now do end of trace if apply pt not at end of gate.
            iobeg = ioapp
            ioend = 1 +length -1
            if (pwr_dbl .gt. AGC_EPS ) then
               scalval = nnrm*normvalue/pwr_dbl
               do k=iobeg,ioend
                  sclout(k) = scalval
	       enddo
            else
               do k = iobeg,ioend
                  sclout(k) = deadgate
	       enddo
            endif
      else 
c
c     Single precision loop for amp scaling   
c     
c Form first scaler sum.
            pwr_sngl = 0.0
            nnrm = 0

c Calculate pwr from live to live+gate.
            ibeg = 1 + ist -1
            iend = 1 + min (ist + gate -1, lst-1)

            do k=ibeg,iend
               if(vecin(k).gt.AGC_EPS .or. vecin(k).lt.-AGC_EPS) then
                 pwr_sngl = pwr_sngl+abs(vecin(k))
                 nnrm = nnrm +1 
               endif
            enddo

c Store first scaler in scaletrace, from begining of gate to applypt.
            iobeg = 1+ist -1
            ioend = 1+min(ist-1+applypt,lst-1)

            if (pwr_sngl .gt. AGC_EPS) then
               scalval = nnrm*normvalue/pwr_sngl 
               do k=iobeg,ioend
                  sclout(k) = scalval
               enddo 
            else
               do k=iobeg,ioend
                  sclout(k) = deadgate
               enddo
            endif

c Now do middle portion of trace.
            ibeg =  1 + ist -1
            iend =  1 + lst-gate-2
            ioapp = ioend +1

            do k= ibeg,iend
               rout = vecin(k)
               rin  = vecin(k +gate +1)
               if(rout.gt.AGC_EPS .or. rout.lt.-AGC_EPS) nnrm = nnrm - 1
               if(rin .gt.AGC_EPS .or. rout.lt.-AGC_EPS) nnrm = nnrm + 1
               pwr_sngl = pwr_sngl-abs(rout)+abs(rin)
               if (pwr_sngl .gt. AGC_EPS) then
                 sclout(ioapp) = nnrm*normvalue/pwr_sngl
               else
                 sclout(ioapp) = deadgate
               endif
               ioapp = ioapp + 1
            enddo

c Now do end of trace if apply pt not at end of gate.
            iobeg = ioapp
            ioend = 1 +length -1
            if (pwr_sngl .gt. AGC_EPS ) then
               scalval = nnrm*normvalue/pwr_sngl
               do k=iobeg,ioend
                  sclout(k) = scalval
	       enddo
            else
               do k = iobeg,ioend
                  sclout(k) = deadgate
	       enddo
            endif
       endif 

       return 
       end







