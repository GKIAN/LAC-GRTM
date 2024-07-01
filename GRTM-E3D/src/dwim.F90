!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 17 Aug 2021 10:24:52 AM CST
!-------------------------------------------------------------------------------

#define CAL_OMEGA(f) 2.0_MK * pi * f - feps * (0.0_MK, 1.0_MK)

#ifdef FFTW

#ifdef DOUBLETYPE
#ifdef __GFORTRAN__
#define FFTW_(a) dfftw_/**/a
#else
#define FFTW_(a) dfftw_ ## a
#endif
#else
#ifdef __GFORTRAN__
#define FFTW_(a) sfftw_/**/a
#else
#define FFTW_(a) sfftw_ ## a
#endif
#endif

#endif

module dwimMod

  !$ use omp_lib
  use comm
  use math
  use paraMod, only: ckwlt, tRec, faLim, kfCri, kLim, dt => dtRec, dkLim, &
    & sxyz, rxyz, fxyz, mxyz, svInty, swType, swTime, swFreq, swFile, wkTime, &
    & isForce, lSrc, lRec, z, beta, alpha
  use grtcMod
  implicit none
  private

#ifdef FFTW
      integer(kind = DP) :: prc, pcr
      include 'fftw3.f'
#endif

  integer, public :: ntRec

  integer :: nt, nf, ivMin
  real(kind = MK) :: df, dk
  real(kind = MK) :: fMax, feps, dmin

  integer :: nig
  real(kind = MK) :: r, theta
  real(kind = MK) :: ISH1, IPS1, IPS0
  real(kind = MK) :: RSH2, RSH1, RPS2, RPS1, RPS01, RPS02

  integer :: fiEx = -1
  complex(kind = MK), allocatable :: Swf(:)
  real(kind = MK), allocatable :: Swt(:)
  real(kind = MK), allocatable :: wabs(:)

  integer, parameter :: nIntg = 30, npt = 10
  logical :: suga(nIntg, 2)
  integer :: ipts(nIntg, 2)
  real(kind = MK) :: vpts(nIntg, 2, npt), suIntg(nIntg, 2, 3)
  complex(kind = MK) :: uIntg(nIntg), duIntg(nIntg)
  !$OMP THREADPRIVATE(suga, ipts, vpts, suIntg, uIntg, duIntg)

  complex(kind = MK), allocatable :: iUr(:), iUt(:), iUz(:)
  complex(kind = MK), allocatable :: iTrz(:), iTtz(:), iTzz(:)
  complex(kind = MK), allocatable :: iTrr(:), iTrt(:), iTtt(:)
  real(kind = MK), allocatable :: u(:)
  real(kind = MK), allocatable, public :: t(:)
  real(kind = MK), allocatable, public :: ur(:), ut(:), ux(:), uy(:), uz(:)
  real(kind = MK), allocatable, public :: trz(:), ttz(:), txz(:), tyz(:), tzz(:)
  real(kind = MK), allocatable, public :: trr(:), trt(:), ttt(:), &
    & txx(:), txy(:), tyy(:)

  public dwimInitialize, dwimRun, dwimTransform, dwimFinalize

  contains

    subroutine dwimInitialize()
      real(kind = MK) :: L, vMin, vMax, aLim
      complex(kind = MK) :: omg
      integer :: m, fileID, iros, iloc(1)
      integer :: i, j, ios
      real(kind = MK) :: rTemp
#ifdef FFTW
      complex(kind = MK), allocatable :: Af(:)
#endif

      ntRec = int(tRec / dt) + 1
      m = ceiling( log10( real(ntRec, kind = MK) ) / log10(2.0_MK) )

      nt = 2 ** m
      nf = nt / 2 + 1
      fMax = 1.0_MK / (2.0_MK * dt)
      df = fMax / (nf - 1)
      feps = pi / ((nt - 1) * dt)

      ! azimuth angle: between the line from epicenter to station with the north
      ! under the right-hand rule: x -> North, y -> East, z -> Down
      r = sqrt( (rxyz(1) - sxyz(1)) ** 2 + (rxyz(2) - sxyz(2)) ** 2 )
      if(r == 0.0_MK) then
        theta = 0.0_MK
      else
        if(rxyz(2) >= sxyz(2)) then
          theta = acos( (rxyz(1) - sxyz(1)) / r )
        else
          theta = - acos( (rxyz(1) - sxyz(1)) / r )
        end if
      end if

      vMin = minval(beta)
      vMax = maxval(alpha)
      !=> space period length: it should be so long that wave won't arrive at
      ! any receiver in the time window?
      L = vMax * nt * dt + r + vMax / ((vMin + vMax) / 2.0_MK) &
        & * sqrt(r * r + (rxyz(3) - sxyz(3)) ** 2) + 100.0D+3
      dk = min(dkLim, 2.0_MK * pi / L)

      allocate(Swf(nt)); Swf = (0.0_MK, 0.0_MK)
      allocate(Swt(nt)); Swt = 0.0_MK
      allocate(wabs(nf)); wabs = 0.0_MK

#ifdef FFTW
      allocate(Af(nt)); Af = (0.0_MK, 0.0_MK)
      call FFTW_(plan_dft_r2c_1d) (prc, nt, Swt, Swf, FFTW_ESTIMATE)
      call FFTW_(plan_dft_c2r_1d) (pcr, nt, Swf, Swt, FFTW_ESTIMATE)
#endif

      if(trim(adjustl(swType)) /= 'File') then
        do i = 1, nf
          omg = CAL_OMEGA( df * (i - 1) )
          Swf(i) = mathWavelet(omg, swType, swTime, swFreq)
        end do
        if(ckwlt) then
#ifdef FFTW
          Af = Swf
          call FFTW_(execute_dft_c2r) (pcr, Af, Swt)
          Swt = Swt * df
#else
          Swf(nt:nt/2 + 2:-1) = conjg(Swf(2:nt/2))
          Swt = real( fft(nt, Swf, -1) )
          Swt = Swt * df
#endif
          open(newunit = fileID, file = swFile, status = 'replace')
            write(fileID, '(10X, A, 18X, A)') 'time (s)', 'amplitude'
            do i = 1, nt
              write(fileID, '(2(2X, ES25.17E3))') (i - 1) * dt, Swt(i)
            end do
          close(fileID)
        end if
      else
        i = 0
        open(newunit = fileID, file = swFile, status = 'old')
          read(fileID, *, iostat = ios)
          do while(ios == 0 .and. i < nt)
            i = i + 1
            read(fileID, *, iostat = ios) rTemp, Swt(i)
          end do
        close(fileID)
        if(i /= nt) then
          call commErrorExcept(FAIL2CHECK, 'Insufficient data length in ' &
            & // 'the source wavelet file <' // trim(adjustl(swFile)) // '>.')
        else
          write(*, '(A, I0, A)') 'Read the first ', nt, ' data points from ' &
            & // 'the source wavelet file <' // trim(adjustl(swFile)) // '>.'
#ifdef FFTW
          call FFTW_(execute_dft_r2c) (prc, Swt, Swf)
          Swf = Swf * dt
#else
          Swf = fft(nt, cmplx(Swt, kind = MK), 1)
          Swf = Swf * dt
#endif
        end if
        do i = 1, nf
          omg = CAL_OMEGA( df * (i - 1) )
        end do
      end if

      fiEx = nf
      wabs = abs(Swf(:nf))
      iloc = maxloc(wabs)
      aLim = wabs(iloc(1)) * faLim
      do i = nf, iloc(1), -1
        if(wabs(i) >= aLim) exit
      end do
      fiEx = min(i + 1, nf)

      if(lRec < lSrc) then
        iros = min(lRec + 1, lSrc - 1)
        iloc = minloc( beta(iros:lSrc) )
        ivMin = iloc(1) + iros - 1
        dmin = min( sxyz(3) - z(lSrc - 1), z(lSrc - 1) - rxyz(3) )
        do j = lRec, lSrc - 2
          dmin = min( dmin, z(j + 1) - z(j) )
        end do
      else if(lRec == lSrc) then
        ivMin = lSrc
        dmin = abs(sxyz(3) - rxyz(3))
      else if(lRec > lSrc) then
        iros = max(lSrc + 1, lRec - 1)
        iloc = minloc( beta(lSrc:iros) )
        ivMin = iloc(1) + lSrc - 1
        dmin = min( z(lSrc) - sxyz(3), rxyz(3) - z(lSrc) )
        do j = lSrc + 1, lRec - 1
          dmin = min( dmin, z(j) - z(j - 1) )
        end do
      end if

      if(isForce) then
        nig = 15
        ISH1 = fxyz(2) * cos(theta) - fxyz(1) * sin(theta)
        IPS1 = fxyz(1) * cos(theta) + fxyz(2) * sin(theta)
        IPS0 = fxyz(3)
      else
        nig = 30
        RSH2  = ( mxyz(2, 2) - mxyz(1, 1) ) / 2.0_MK * sin(2.0_MK * theta) &
          & + mxyz(1, 2) * cos(2.0_MK * theta)
        RSH1  = mxyz(1, 3) * sin(theta) - mxyz(2, 3) * cos(theta)
        RPS2  = - ( mxyz(2, 2) - mxyz(1, 1) ) / 2.0_MK * cos(2.0_MK * theta) &
          & + mxyz(1, 2) * sin(2.0_MK * theta)
        RPS1  = mxyz(1, 3) * cos(theta) + mxyz(2, 3) * sin(theta)
        RPS01 = ( mxyz(1, 1) + mxyz(2, 2) ) / 2.0_MK
        RPS02 = mxyz(3, 3)
      end if

      suga = .true.
      vpts = 0.0_MK

      allocate(iUr(nt)); iUr = (0.0_MK, 0.0_MK)
      allocate(iUt(nt)); iUt = (0.0_MK, 0.0_MK)
      allocate(iUz(nt)); iUz = (0.0_MK, 0.0_MK)
      allocate(iTrz(nt)); iTrz = (0.0_MK, 0.0_MK)
      allocate(iTtz(nt)); iTtz = (0.0_MK, 0.0_MK)
      allocate(iTzz(nt)); iTzz = (0.0_MK, 0.0_MK)
      allocate(iTrr(nt)); iTrr = (0.0_MK, 0.0_MK)
      allocate(iTrt(nt)); iTrt = (0.0_MK, 0.0_MK)
      allocate(iTtt(nt)); iTtt = (0.0_MK, 0.0_MK)
      allocate(u(nt)); u = 0.0_MK

      allocate(t(ntRec)); t = [ ( (i - 1) * dt, i = 1, ntRec ) ]
      allocate(ur(ntRec)); ur = 0.0_MK
      allocate(ut(ntRec)); ut = 0.0_MK
      allocate(ux(ntRec)); ux = 0.0_MK
      allocate(uy(ntRec)); uy = 0.0_MK
      allocate(uz(ntRec)); uz = 0.0_MK
      allocate(trz(ntRec)); trz = 0.0_MK
      allocate(ttz(ntRec)); ttz = 0.0_MK
      allocate(txz(ntRec)); txz = 0.0_MK
      allocate(tyz(ntRec)); tyz = 0.0_MK
      allocate(tzz(ntRec)); tzz = 0.0_MK
      allocate(trr(ntRec)); trr = 0.0_MK
      allocate(trt(ntRec)); trt = 0.0_MK
      allocate(ttt(ntRec)); ttt = 0.0_MK
      allocate(txx(ntRec)); txx = 0.0_MK
      allocate(txy(ntRec)); txy = 0.0_MK
      allocate(tyy(ntRec)); tyy = 0.0_MK

#ifdef FFTW
      if(allocated(Af)) deallocate(Af)
#endif

#ifdef DEBUG
      write(*, '(A, 2(1X, G0))') 'dwimInit: nt =', m, nt
      write(*, '(A, 2(1X, G0))') 'dwimInit: df =', fMax, df
      write(*, '(A, 2(1X, G0))') 'dwimInit: dk =', L, dk
#endif
    end subroutine dwimInitialize

    subroutine dwimRun()
      real(kind = MK) :: k, kCri
      complex(kind = MK) :: omg, Sw2p
      integer :: i, j, jj, ni
      logical :: isMain = .true.

      ni = 0
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, jj, k, omg, kCri, Sw2p, &
      !$OMP   & isMain) COPYIN(suga, vpts)
      !$ isMain = (omp_get_thread_num() == 0)

      !$OMP DO SCHEDULE(GUIDED, 4)
      do i = 1, fiEx
#ifdef PROGBAR
        if(isMain) &
          & call commProgressBar(ni, fiEx, 'run k-integration by frequency')
#else
        if(mod(i, 10) == 0) write(*, '(A, I0, A, I0, A)') &
          & 'Run k-integration for frequency point ', i, '/', fiEx, ' ...'
#endif
        omg = CAL_OMEGA( df * (i - 1) )
        call grtcSetMedia(omg)
        Sw2p = svInty * Swf(i) * exp( - omg * wkTime * (0.0_MK, 1.0_MK) )
        Sw2p = Sw2p / ( 2.0_MK * pi )
        kCri = kfCri * kmax(i, omg)

        ! when under the critical k-value, integrate by trapezoidal rule
        k = dk
        suIntg = 0.0_MK
        do while(.true.)
          duIntg = integrand(k, omg) * dk
          suIntg(:, 1, 2) = suIntg(:, 1, 2) + real(duIntg)
          suIntg(:, 2, 2) = suIntg(:, 2, 2) + imag(duIntg)
          if(k > kCri) exit
          suIntg(:, :, 1) = suIntg(:, :, 2)
          k = k + dk
        end do

        ! when beyond the critical k-value, integrate by peak-trough averaging
        ipts = 0
        do while(.true.)
          if(k > kLim) then
            call commErrorExcept(NOFATALERR, &
              & 'Integrate to the limited k-value. Result NOT reliable.')
            exit
          end if

          k = k + dk
          duIntg = integrand(k, omg) * dk
          suIntg(:, 1, 3) = suIntg(:, 1, 2) + real(duIntg)
          suIntg(:, 2, 3) = suIntg(:, 2, 2) + imag(duIntg)

          do j = 1, 2
            do jj = 1, nig
              suga(jj, j) = ptamGather(suIntg(jj, j, :), ipts(jj, j), &
                & vpts(jj, j, :))
            end do
          end do
          if(all(suga)) exit

          suIntg(:, :, 1) = suIntg(:, :, 2)
          suIntg(:, :, 2) = suIntg(:, :, 3)
        end do

        do j = 1, 2
          do jj = 1, nig
            call ptamReduce(ipts(jj, j), vpts(jj, j, :))
          end do
        end do
        uIntg = cmplx(vpts(:, 1, 1), vpts(:, 2, 1), kind = MK)

        ! assemble the spectra of displacements
        if(isForce) then
          iUr (i) =   Sw2p * ( uIntg(1 ) + uIntg(2 ) )
          iUt (i) =   Sw2p *   uIntg(3 )
          iUz (i) = - Sw2p * ( uIntg(4 ) + uIntg(5 ) )

          iTrz(i) =   Sw2p * ( uIntg(6 ) + uIntg(7 ) )
          iTtz(i) =   Sw2p *   uIntg(8 )
          iTzz(i) = - Sw2p * ( uIntg(9 ) + uIntg(10) )

          iTrr(i) =   Sw2p * ( uIntg(11) + uIntg(12) )
          iTrt(i) =   Sw2p *   uIntg(13)
          iTtt(i) = - Sw2p * ( uIntg(14) + uIntg(15) )

        else
          iUr (i) =   Sw2p * ( uIntg(1 ) + uIntg(2 ) + uIntg(3 ) + uIntg(4 ) )
          iUt (i) =   Sw2p * ( uIntg(5 ) + uIntg(6 ) )
          iUz (i) = - Sw2p * ( uIntg(7 ) + uIntg(8 ) + uIntg(9 ) + uIntg(10) )

          iTrz(i) =   Sw2p * ( uIntg(11) + uIntg(12) + uIntg(13) + uIntg(14) )
          iTtz(i) =   Sw2p * ( uIntg(15) + uIntg(16) )
          iTzz(i) = - Sw2p * ( uIntg(17) + uIntg(18) + uIntg(19) + uIntg(20) )

          iTrr(i) =   Sw2p * ( uIntg(21) + uIntg(22) + uIntg(23) + uIntg(24) )
          iTrt(i) =   Sw2p * ( uIntg(25) + uIntg(26) )
          iTtt(i) = - Sw2p * ( uIntg(27) + uIntg(28) + uIntg(29) + uIntg(30) )
        end if

#ifndef STRAIN
        iTrr(i) = iTrr(i) + eta * iTzz(i)
        iTtt(i) = iTtt(i) + eta * iTzz(i)
#endif
        ni = ni + 1
      end do
      !$OMP END DO

#ifdef PROGBAR
      if(isMain) &
        & call commProgressBar(fiEx, fiEx, 'run k-integration by frequency')
#endif

      !$OMP END PARALLEL

#ifdef FFTW
      call FFTW_(execute_dft_c2r) (pcr, iUr , u); ur  = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iUt , u); ut  = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iUz , u); uz  = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTrz, u); trz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTtz, u); ttz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTzz, u); tzz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTrr, u); trr = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTrt, u); trt = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTtt, u); ttt = u(1:ntRec) * df * exp(feps * t)
#else
      !=> For the DFT, we have $ X_{N - m} = X_{-m} $. And specially for purely
      ! real input, the output is Hermitian-symmetric, that is, the negative-
      ! frequency terms are just the conjugates of the corresponding positive-
      ! frequency terms.
      iUr (nt:nt/2 + 2:-1) = conjg( iUr (2:nt/2) )
      iUt (nt:nt/2 + 2:-1) = conjg( iUt (2:nt/2) )
      iUz (nt:nt/2 + 2:-1) = conjg( iUz (2:nt/2) )
      iTrz(nt:nt/2 + 2:-1) = conjg( iTrz(2:nt/2) )
      iTtz(nt:nt/2 + 2:-1) = conjg( iTtz(2:nt/2) )
      iTzz(nt:nt/2 + 2:-1) = conjg( iTzz(2:nt/2) )
      iTrr(nt:nt/2 + 2:-1) = conjg( iTrr(2:nt/2) )
      iTrt(nt:nt/2 + 2:-1) = conjg( iTrt(2:nt/2) )
      iTtt(nt:nt/2 + 2:-1) = conjg( iTtt(2:nt/2) )
      u = real( fft(nt, iUr , -1) ); ur  = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iUt , -1) ); ut  = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iUz , -1) ); uz  = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTrz, -1) ); trz = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTtz, -1) ); ttz = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTzz, -1) ); tzz = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTrr, -1) ); trr = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTrt, -1) ); trt = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTtt, -1) ); ttt = u(1:ntRec) * df * exp(feps * t)
#endif
    end subroutine dwimRun

    subroutine dwimTransform()
      ux  = ur  * cos(theta) - ut  * sin(theta)
      uy  = ut  * cos(theta) + ur  * sin(theta)
      txz = trz * cos(theta) - ttz * sin(theta)
      tyz = ttz * cos(theta) + trz * sin(theta)
      txx = trr * cos(theta) * cos(theta) + ttt * sin(theta) * sin(theta) &
          & - 2.0_MK * trt * cos(theta) * sin(theta)
      tyy = ttt * cos(theta) * cos(theta) + trr * sin(theta) * sin(theta) &
          & + 2.0_MK * trt * cos(theta) * sin(theta)
      txy = trt * ( cos(theta) * cos(theta) - sin(theta) * sin(theta) ) &
          & + (trr - ttt) * cos(theta) * sin(theta)
    end subroutine dwimTransform

    subroutine dwimFinalize()
      if(allocated(wabs)) deallocate(wabs)

      if(allocated(iUr)) deallocate(iUr)
      if(allocated(iUt)) deallocate(iUt)
      if(allocated(iUz)) deallocate(iUz)
      if(allocated(iTrz)) deallocate(iTrz)
      if(allocated(iTtz)) deallocate(iTtz)
      if(allocated(iTzz)) deallocate(iTzz)
      if(allocated(iTrr)) deallocate(iTrr)
      if(allocated(iTrt)) deallocate(iTrt)
      if(allocated(iTtt)) deallocate(iTtt)
      if(allocated(u)) deallocate(u)

      if(allocated(t)) deallocate(t)
      if(allocated(ur)) deallocate(ur)
      if(allocated(ut)) deallocate(ut)
      if(allocated(ux)) deallocate(ux)
      if(allocated(uy)) deallocate(uy)
      if(allocated(uz)) deallocate(uz)
      if(allocated(trz)) deallocate(trz)
      if(allocated(ttz)) deallocate(ttz)
      if(allocated(txz)) deallocate(txz)
      if(allocated(tyz)) deallocate(tyz)
      if(allocated(tzz)) deallocate(tzz)
      if(allocated(trr)) deallocate(trr)
      if(allocated(trt)) deallocate(trt)
      if(allocated(ttt)) deallocate(ttt)
      if(allocated(txx)) deallocate(txx)
      if(allocated(txy)) deallocate(txy)
      if(allocated(tyy)) deallocate(tyy)

#ifdef FFTW
      call FFTW_(destroy_plan) (prc)
      call FFTW_(destroy_plan) (pcr)
#endif
    end subroutine dwimFinalize

    real(kind = MK) function kmax(fi, omg)
      integer, intent(in) :: fi
      complex(kind = MK), value :: omg
      real(kind = MK) :: ka

      if(fi <= 3) omg = CAL_OMEGA(df * 3)
      ka = abs(omg/gbeta(ivMin))

      if(dmin < 0.2D+3) then
        if(fi < 10) then
          kmax = min(ka + 2.0D-3, 20 * ka)
        else
          kmax = min(ka + 2.0D-3, 2  * ka)
        end if
      else
        kmax = sqrt( (3.0_MK/dmin) ** 2 + ka * ka )
      end if
      kmax = max(kmax, ka + 0.5D-3)
    end function kmax

    function integrand(k, omg) result(intg)
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg
      complex(kind = MK) :: intg(nIntg)
      real(kind = MK) :: kr, jn(0:4), ej(0:6)

      call grtcSetEigen(k, omg)
#ifdef SH
      call grtcCoefficientSH()
      call grtcKernelCoeffiSH()
      if(isForce) then
        call grtcKernelForceSH()
      else
        call grtcKernelMomentSH(k)
      end if
#endif
#ifdef PS
      call grtcCoefficientPS()
      call grtcKernelCoeffiPS()
      if(isForce) then
        call grtcKernelForcePS(k, omg)
      else
        call grtcKernelMomentPS(k, omg)
      end if
#endif

      kr = k * r
      jn(0) = bessel_j0(kr)
      jn(1) = bessel_j1(kr)
      jn(2) = bessel_jn(2, kr)
      jn(3) = bessel_jn(3, kr)
      ej(0) = - jn(1)
      ej(1) = ( jn(0) - jn(2) ) / 2.0_MK
      ej(2) = ( jn(0) + jn(2) ) / 2.0_MK
      ej(3) = ( jn(1) - jn(3) ) / 2.0_MK
      ej(4) = ( jn(1) + jn(3) ) / 2.0_MK

      if(isForce) then
        !ref.: eqs. (4-7a, b and c)
        intg(1 ) = IPS1 * k * ( gSH1    * jn(0)  + (gPS1(1) - gSH1) * ej(1) )
        intg(2 ) = IPS0 * k *   gPS0(1) * ej(0)
        intg(3 ) = ISH1 * k * ( gPS1(1) * jn(0)  + (gSH1 - gPS1(1)) * ej(1) )
        intg(4 ) = IPS1 * k *   gPS1(2) * jn(1)
        intg(5 ) = IPS0 * k *   gPS0(2) * jn(0)

#ifndef STRAIN
        intg(6 ) = IPS1 * k * ( tSH1    * jn(0)  + (tPS1(1) - tSH1) * ej(1) )
        intg(7 ) = IPS0 * k *   tPS0(1) * ej(0)
        intg(8 ) = ISH1 * k * ( tPS1(1) * jn(0)  + (tSH1 - tPS1(1)) * ej(1) )
        intg(9 ) = IPS1 * k *   tPS1(2) * jn(1)
        intg(10) = IPS0 * k *   tPS0(2) * jn(0)

        intg(11) = - IPS1 * k * ( gPS1(1) * k * ( muj * ej(3) + kap * jn(1) ) &
                              & + gSH1    * k *   muj * ej(4) )
        intg(12) =   IPS0 * k *   gPS0(1) * k * ( muj * jn(2) - kap * jn(0) )
        intg(13) = - ISH1 * k * ( gPS1(1) * k *   muj * ej(4) &
                              & + gSH1    * k *   muj * ej(3) )
        intg(14) = - IPS1 * k * ( gPS1(1) * k * ( muj * ej(3) - kap * jn(1) ) &
                              & + gSH1    * k *   muj * ej(4) )
        intg(15) =   IPS0 * k *   gPS0(1) * k * ( muj * jn(2) + kap * jn(0) )
#else
        intg(6 ) =   IPS1 * k * ( tSH1    * ej(2) &
                            & + ( tPS1(1) - gPS1(2) * k ) * ej(1) ) / 2.0_MK
        intg(7 ) = - IPS0 * k * ( tPS0(1) - gPS0(2) * k ) * jn(1)   / 2.0_MK
        intg(8 ) =   ISH1 * k * ( tSH1    * ej(1) &
                            & + ( tPS1(1) - gPS1(2) * k ) * ej(2) ) / 2.0_MK
        intg(9 ) =   IPS1 * k *   tPS1(2) * jn(1)
        intg(10) =   IPS0 * k *   tPS0(2) * jn(0)

        intg(11) = - IPS1 * k * ( gPS1(1) * k * jn(1) &
                            & - ( gPS1(1) - gSH1 ) * k * ej(4) / 2.0_MK )
        intg(12) = - IPS0 * k *   gPS0(1) * k * ej(1)
        intg(13) = - ISH1 * k * ( gSH1    * k * jn(1) / 2.0_MK &
                            & + ( gPS1(1) - gSH1 ) * k * ej(4) / 2.0_MK )
        intg(14) =   IPS1 * k * ( gPS1(1) - gSH1 ) * k * ej(4) / 2.0_MK
        intg(15) =   IPS0 * k *   gPS0(1) * k * ej(2)
#endif

      else
        jn(4) = bessel_jn(4, kr)
        ej(5) = ( jn(0) - jn(4) ) / 2.0_MK
        ej(6) = ( jn(0) + jn(4) ) / 2.0_MK

        !ref.: eqs. (4-9a, b and c)
        intg(1 ) =   RPS2  * k * ( qSH2    * jn(1)  + (qPS2(1) - qSH2) * ej(3) )
        intg(2 ) = - RPS1  * k * ( qSH1    * jn(0)  + (qPS1(1) - qSH1) * ej(1) )
        intg(3 ) =   RPS02 * k *   qPS0(1) * ej(0)
        intg(4 ) = - RPS01 * k *   qPS2(1) * ej(0)
        intg(5 ) =   RSH2  * k * ( qPS2(1) * jn(1)  + (qSH2 - qPS2(1)) * ej(3) )
        intg(6 ) =   RSH1  * k * ( qPS1(1) * jn(0)  + (qSH1 - qPS1(1)) * ej(1) )
        intg(7 ) =   RPS2  * k *   qPS2(2) * jn(2)
        intg(8 ) = - RPS1  * k *   qPS1(2) * jn(1)
        intg(9 ) =   RPS02 * k *   qPS0(2) * jn(0)
        intg(10) = - RPS01 * k *   qPS2(2) * jn(0)

#ifndef STRAIN
        intg(11) =   RPS2  * k * ( bSH2    * jn(1)  + (bPS2(1) - bSH2) * ej(3) )
        intg(12) = - RPS1  * k * ( bSH1    * jn(0)  + (bPS1(1) - bSH1) * ej(1) )
        intg(13) =   RPS02 * k *   bPS0(1) * ej(0)
        intg(14) = - RPS01 * k *   bPS2(1) * ej(0)
        intg(15) =   RSH2  * k * ( bPS2(1) * jn(1)  + (bSH2 - bPS2(1)) * ej(3) )
        intg(16) =   RSH1  * k * ( bPS1(1) * jn(0)  + (bSH1 - bPS1(1)) * ej(1) )
        intg(17) =   RPS2  * k *   bPS2(2) * jn(2)
        intg(18) = - RPS1  * k *   bPS1(2) * jn(1)
        intg(19) =   RPS02 * k *   bPS0(2) * jn(0)
        intg(20) = - RPS01 * k *   bPS2(2) * jn(0)

        intg(21) =   RPS2  * k * ( qPS2(1) * k * ( muj * ej(6) - kap * jn(2) ) &
                               & + qSH2    * k *   muj * ej(5) )
        intg(22) =   RPS1  * k * ( qPS1(1) * k * ( muj * ej(3) + kap * jn(1) ) &
                               & + qSH1    * k *   muj * ej(4) )
        intg(23) =   RPS02 * k *   qPS0(1) * k * ( muj * jn(2) - kap * jn(0) )
        intg(24) = - RPS01 * k *   qPS2(1) * k * ( muj * jn(2) - kap * jn(0) )
        intg(25) =   RSH2  * k * ( qPS2(1) * k *   muj * ej(5) &
                               & + qSH2    * k *   muj * ej(6) )
        intg(26) = - RSH1  * k * ( qPS1(1) * k *   muj * ej(4) &
                               & + qSH1    * k *   muj * ej(3) )
        intg(27) =   RPS2  * k * ( qPS2(1) * k * ( muj * ej(6) + kap * jn(2) ) &
                               & + qSH2    * k *   muj * ej(5) )
        intg(28) =   RPS1  * k * ( qPS1(1) * k * ( muj * ej(3) - kap * jn(1) ) &
                               & + qSH1    * k *   muj * ej(4) )
        intg(29) =   RPS02 * k *   qPS0(1) * k * ( muj * jn(2) + kap * jn(0) )
        intg(30) = - RPS01 * k *   qPS2(1) * k * ( muj * jn(2) + kap * jn(0) )
#else
        intg(11) =   RPS2  * k * ( bSH2    * ej(4) &
                             & + ( bPS2(1) - qPS2(2) * k ) * ej(3) ) / 2.0_MK
        intg(12) = - RPS1  * k * ( bSH1    * ej(2) &
                             & + ( bPS1(1) - qPS1(2) * k ) * ej(1) ) / 2.0_MK
        intg(13) = - RPS02 * k * ( bPS0(1) - qPS0(2) * k ) * jn(1) / 2.0_MK
        intg(14) =   RPS01 * k * ( bPS2(1) - qPS2(2) * k ) * jn(1) / 2.0_MK
        intg(15) =   RSH2  * k * ( bSH2    * ej(3) &
                             & + ( bPS2(1) - qPS2(2) * k ) * ej(4) ) / 2.0_MK
        intg(16) =   RSH1  * k * ( bSH1    * ej(1) &
                             & + ( bPS1(1) - qPS1(2) * k ) * ej(2) ) / 2.0_MK
        intg(17) =   RPS2  * k *   bPS2(2) * jn(2)
        intg(18) = - RPS1  * k *   bPS1(2) * jn(1)
        intg(19) =   RPS02 * k *   bPS0(2) * jn(0)
        intg(20) = - RPS01 * k *   bPS2(2) * jn(0)

        intg(21) =   RPS2  * k * ( qPS2(1) * k * ej(1) &
                             & - ( qPS2(1) - qSH2 ) * k * ej(5) / 2.0_MK )
        intg(22) =   RPS1  * k * ( qPS1(1) * k * jn(1) &
                             & - ( qPS1(1) - qSH1 ) * k * ej(4) / 2.0_MK )
        intg(23) = - RPS02 * k *   qPS0(1) * k * ej(1)
        intg(24) =   RPS01 * k *   qPS2(1) * k * ej(1)
        intg(25) =   RSH2  * k * ( qPS2(1) * k * jn(0) &
                             & - ( qPS2(1) - qSH2 ) * k * ej(6) ) / 2.0_MK
        intg(26) = - RSH1  * k * ( qSH1    * k * jn(1) &
                             & + ( qPS1(1) - qSH1 ) * k * ej(4) ) / 2.0_MK
        intg(27) =   RPS2  * k * ( qPS2(1) * k * ej(2) &
                             & - ( qPS2(1) - qSH2 ) * k * ej(5) / 2.0_MK )
        intg(28) = - RPS1  * k * ( qPS1(1) - qSH1 ) * k * ej(4) / 2.0_MK
        intg(29) =   RPS02 * k *   qPS0(1) * k * ej(2)
        intg(30) = - RPS01 * k *   qPS2(1) * k * ej(2)
#endif
      end if
    end function integrand

    logical function ptamGather(S, ipt, vpt) result(gotAll)
      real(kind = MK), intent(in) :: S(:)
      integer, intent(inout) :: ipt
      real(kind = MK), intent(inout) :: vpt(:)
      real(kind = MK) :: toput, a(3)
      logical :: ifput
      ifput = .false.
      if(S(1) == S(2) .and. S(2) == S(3)) then
        ifput = .true.
        toput = S(1)
      else if((S(2) - S(1)) * (S(2) - S(3)) > 0.0_MK .or. S(2) == S(3)) then
        ifput = .true.
        a(1) = 2.0_MK * S(3) - 4.0_MK * S(2) + 2.0_MK * S(1)
        if(a(1) == 0.0_MK) then
          !=> if there is only a tiny difference between S(1) with S(2) when
          ! S(2) is equal to S(3)
          toput = S(2)
        else
          a(2) = 4.0_MK * S(2) - S(3) - 3.0_MK * S(1)
          a(3) = S(1)
          toput = a(3) - a(2) * a(2) / (4.0_MK * a(1))
        end if
      end if
      if(ifput) then
        ipt = ipt + 1
        if(ipt > npt) then
          vpt(1:npt - 1) = vpt(2:npt)
          ipt = npt
        end if
        vpt(ipt) = toput
      end if
      gotAll = (ipt >= npt)
    end function ptamGather

    subroutine ptamReduce(ipt, vpt)
      integer, intent(in) :: ipt
      real(kind = MK), intent(inout) :: vpt(:)
      integer :: i, j
      do i = ipt - 1, 1, -1
        do j = 1, i, +1
          vpt(j) = (vpt(j) + vpt(j + 1)) / 2.0_MK
        end do
      end do
    end subroutine ptamReduce

end module dwimMod

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
