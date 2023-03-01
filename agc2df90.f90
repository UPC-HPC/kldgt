subroutine agc2df90(winagc2d,taptagc2d,sclmaxagc2d,insamppad,isi,intin,scalaragc2d,datapad)

  real, dimension(insamppad,intin) , intent(INOUT) :: datapad
  real, dimension(insamppad), intent(OUT) :: scalaragc2d

  real :: winagc2d,taptagc2d,rsi
  integer :: insamppad,isi

  integer :: iagc2d,idnnct,idntapt,idnwinmoveup
  integer :: idnnwin,ifsamp,imsamp,ilsamp

  real, allocatable :: dnwval(:)
  real, allocatable :: dnposin(:),dnposou(:)

  iagc2d = 1

  rsi = real(isi)/1000.0

! Variables if we need to normalise the data....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (iagc2d == 1) then
   idnnct  = winagc2d /rsi
   idntapt = taptagc2d/rsi

   idnwinmoveup = idnnct-idntapt
   idnnwin  = 0

   do
      idnnwin = idnnwin + 1
      ifsamp = 1+(idnnwin - 1)*idnwinmoveup
      ilsamp = ifsamp + idnnct - 1
      imsamp = (ifsamp+ilsamp)/2
      if (imsamp + idnwinmoveup > insamppad) exit
   enddo
   
   if (idnnwin<5) then
      write(*,*) 'ERROR : The interpolation code used for option AGC2D '
      write(*,*) 'requires at least 5 points for interpolation.'
      write(*,*) 'Please code a smaller value for the time window length WINT' 
      call abort
   end if

   allocate(dnwval(idnnwin),stat=istat)!;write (9,*) "Alloc wval",inwin,istat
      
   allocate(dnposin(idnnwin),stat=istat)!;write (9,*) "Alloc posin",inwin,istat

   allocate(dnposou(insamppad),stat=istat)!;write (9,*) "Alloc posou",insamppad,istat

   ! Initisliase values for sparseness normalisation if required.
   do i = 1, insamppad
      dnposou(i) = real(i) * rsi
   enddo
   do i = 1, idnnwin
      ifsamp = 1+(i - 1)*idnwinmoveup
      ilsamp = ifsamp + idnnct - 1
      dnposin(i) = rsi*0.5*real(ifsamp+ilsamp)
   enddo
endif

! Balance the input data if required
if (iagc2d == 1) then
   call agc2d_balance(isi,               & ! IN Sample interval; microseconds
                          insamppad,         & ! IN Number of samples
                          intin,             & ! IN Number of ps
                          sclmaxagc2d,        & ! IN Maximum scalar
                          idnnct,            & ! IN Window size for weight balancing
                          idnwinmoveup,      & ! IN Window moveup
                          idnnwin,           & ! IN Number of windows for normalisation
                          dnposin,           & ! IN Input positions
                          dnposou,           & ! IN Output positions
                          dnwval,            & ! SC Amplitude values for each window.
                          scalaragc2d,         & ! SC Interpolated amplitude values
                          datapad,           & ! IN Taup domain
                          datapad)             ! OU Sparseness weights

endif

!!!!scalaragc2d(1:insamppad)=1.0

deallocate(dnwval)
deallocate(dnposin)
deallocate(dnposou)

end subroutine agc2df90

subroutine agc2d_balance(isi, insamppad, ipfold, &
     rscalmax, iwinsize, iwinmoveup, inwin, posin, posou, wval, wvalint, taup, sparse)

implicit none

integer :: isi, insamppad, ipfold
real    :: rscalmax
integer :: iwinsize, iwinmoveup, inwin
real    :: posin(inwin), posou(insamppad), wval(inwin), wvalint(insamppad)
real    :: taup(insamppad,ipfold),sparse(insamppad,ipfold)
!!!!!!!!!!!!!!!!!!!!
integer :: iwin, ifsamp, ilsamp, insamp, ip, isamp, inval, ifsampint, ilsampint
real    :: rtmp
real    :: rsi
!!!!!!!!!!!!!!!!!!!!

rsi = real(isi)/1000.0

ifsampint = ceiling(posin(1)    /rsi)
ilsampint = ceiling(posin(inwin)/rsi)

! If necessary, apply some balancing to the sparseness weights
! Calculate value for each window
DO IWIN = 1, INWIN
   IFSAMP = 1+(IWIN-1)*IWINMOVEUP
   ILSAMP = IFSAMP + IWINSIZE - 1
   ILSAMP = MIN(ILSAMP,INSAMPPAD)
   INSAMP = ILSAMP-IFSAMP + 1
   
   wval(iwin) = 0.
   
   DO IP = 1, IPFOLD
      DO ISAMP = IFSAMP, ILSAMP
         WVAL(IWIN) = WVAL(IWIN) + ABS(SPARSE(ISAMP,IP))
      ENDDO
   ENDDO
   WVAL(IWIN) = WVAL(IWIN)/REAL(IPFOLD*INSAMP)
ENDDO

! Do the cubic spline interpolation for the scaling trace
inval = ilsampint-ifsampint+1
CALL AGC2D_DTIAKIM(POSIN,WVAL,INWIN,POSOU(ifsampint),WVALINT(ifsampint),inval)

! Set the values outside the cubic spline range....
wvalint(1:ifsampint)         = WVAL(1)
wvalint(ilsampint:insamppad) = WVAL(INWIN)
! Smooth the values and invert to find the scalars
rtmp = 0.
DO ISAMP = 1, INSAMPPAD
   ! Calculate the maximum
   wvalint(isamp) = max(wvalint(isamp),0.)
   rtmp = max(rtmp,wvalint(isamp))
enddo
! Now scale the data with limit
rtmp = rtmp / rscalmax
DO ISAMP = 1, INSAMPPAD
   wvalint(isamp) = max(rtmp,wvalint(isamp))
   if (wvalint(isamp) > 0.) then
!      wvalint(isamp) = 1./wvalint(isamp)
      SPARSE(ISAMP,:) = SPARSE(ISAMP,:)/wvalint(isamp)
   ENDIF
enddo

end subroutine agc2d_balance


      SUBROUTINE AGC2D_DTIAKIM (X,Y,N,Z,B,M)
      IMPLICIT NONE

      INTEGER I,IP,J,JU,M,N,NI

      REAL DX,DZ,HT,HTT,P2,P3,S4,S5,T1,T2,TT,XX,YA,YB,YC,YD,Z1,Z2,Z21,Z3,Z32,Z4,Z43,Z5,Z54
      REAL B(M),V(5),X(N),Y(N),Z(M)

      DATA Z1,Z2,Z3,Z4,Z5 / 5*0 /

      DX = 0.;DZ = 0.

      IP = -1
      IF (DX .EQ. 0.) THEN
         V(1)  = X(N-1) - X(N-2)
         V(2)  = X(N) - X(N-1)
         V(3)  = X(3) - X(2)
         V(4)  = X(2) - X(1)
      ELSE
         DO 90 I=1,5
            V(I) = DX
   90    CONTINUE
      ENDIF
      Z1 = (Y(N-1)-Y(N-2)) / V(1)
      Z2 = (Y(N)-Y(N-1)) / V(2)
      Z3 = Z2 * 2. - Z1
      Z4 = Z3 * 2. - Z2
      YC = V(1) * Z3
      YD = V(2) * Z4
      Z1 = (Y(3)-Y(2)) / V(3)
      Z2 = (Y(2)-Y(1)) / V(4)
      Z3 = Z2 * 2. - Z1
      Z4 = Z3 * 2. - Z2
      YB = V(3) * Z3
      YA = V(4) * Z4
      NI = 1
      DO  100 J=1,M
         IF (DZ .NE. 0.) THEN
            TT = Z(1) + (J-1) * DZ
         ELSE
            TT = Z(J)
         ENDIF
         DO  120 I=NI,N
            IF (DX .NE. 0.) THEN
               XX = X(1) + (I-1) * DX
            ELSE
               XX = X(I)
            ENDIF
            HTT = TT - XX
            IF (HTT .LT. 0.)  GOTO 130
  120    CONTINUE
         IF (HTT .GT. 0) THEN
            M = J - 1
         ELSE
            IF (HTT .EQ. 0.) THEN
               M = J
               B(J) = Y(N)
            ENDIF
         ENDIF
         GOTO 150
  130    CONTINUE
         NI = I
         I = I - 1
         IF (DX .NE. 0.) THEN
            HT = HTT + DX
         ELSE
            HT = TT - X(I)
         ENDIF
         IF (I .NE. IP) THEN
            IP = I
            IF (DX .EQ. 0.) THEN
               IF (I .GT. 2) THEN
                  IF (I .LT. (N-2)) THEN
                     DO 110 JU=1,5
                        V(JU) = X(I-2+JU) - X(I-3+JU)
  110                CONTINUE
                  ELSE
                     DO  111 JU=1,3
                        V(JU) = X(I-2+JU) - X(I-3+JU)
  111                CONTINUE
                     V(5) = V(3)
                     V(4) = V(2)
                     IF (I .EQ. (N-2)) V(4) = X(N) - X(N-1)
                  ENDIF
               ELSE
                  DO 112 JU=3,5
                     V(JU) = X(I-2+JU) - X(I-3+JU)
  112             CONTINUE
                  V(1) = V(3)
                  V(2) = V(4)
                  IF (I .EQ. 2) V(2) = X(2) - X(1)
               ENDIF
            ENDIF
            IF (I .GT. 1) THEN
               YA = YB
               IF (I .GT. 2) YA = Y(I-1) - Y(I-2)
               YB = Y(I) - Y(I-1)
            ENDIF
            Z1 = YA / V(1)
            Z2 = YB / V(2)
            Z3 = (Y(I+1) - Y(I)) / V(3)
            IF (I .LT. (N-2)) THEN
               Z4 = Y(I+2) - Y(I+1)
               Z5 = Y(I+3) - Y(I+2)
            ELSE
               IF (I .EQ. (N-2)) THEN
                  Z4 = Y(N) - Y(N-1)
                  Z5 = YC
               ELSE
                  IF (I .EQ. (N-1)) THEN
                     Z4 = YC
                     Z5 = YD
                  ENDIF
               ENDIF
            ENDIF
            Z4  = Z4 / V(4)
            Z5  = Z5 / V(5)
            Z54 = ABS(Z5-Z4)
            Z43 = ABS(Z4-Z3)
            Z32 = ABS(Z3-Z2)
            Z21 = ABS(Z2-Z1)
            S4  = Z43 + Z21
            S5  = Z54 + Z32
            IF (S4 .EQ. 0.) THEN
               T1 = (Z2 + Z3) / 2.
            ELSE
               T1 = (Z43 * Z2 + Z21 * Z3) / S4
            ENDIF
            IF (S5 .EQ. 0.) THEN
               T2 = (Z3 + Z4) / 2.
            ELSE
               T2 = (Z54 * Z3 + Z32 * Z4) / S5
            ENDIF
         ENDIF
         P2      = ( Z3*3. - T1*2. - T2 )/ V(3)
         P3      = ( T1 + T2 - Z3*2. ) / V(3)**2
         B(J) = ((P3 * HT + P2) * HT + T1) * HT + Y(I)
  100 CONTINUE
  150 CONTINUE

      END SUBROUTINE AGC2D_DTIAKIM
