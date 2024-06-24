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

  integer, parameter :: nIntg = 20, npt = 10
  logical :: suga(nIntg, 2)
  integer :: ipts(nIntg, 2)
  real(kind = MK) :: vpts(nIntg, 2, npt), suIntg(nIntg, 2, 3)
  complex(kind = MK) :: uIntg(nIntg), duIntg(nIntg)
  !$OMP THREADPRIVATE(suga, ipts, vpts, suIntg, uIntg, duIntg)

  complex(kind = MK), allocatable :: iUr(:), iUt(:), iUz(:)
  complex(kind = MK), allocatable :: iTr(:), iTt(:), iTz(:)
  real(kind = MK), allocatable :: u(:)
  real(kind = MK), allocatable, public :: t(:)
  real(kind = MK), allocatable, public :: ur(:), ut(:), ux(:), uy(:), uz(:)
  real(kind = MK), allocatable, public :: tr(:), tt(:), tx(:), ty(:), tz(:)

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
        nig = 10
        ISH1 = fxyz(2) * cos(theta) - fxyz(1) * sin(theta)
        IPS1 = fxyz(1) * cos(theta) + fxyz(2) * sin(theta)
        IPS0 = fxyz(3)
      else
        nig = 20
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
      allocate(iTr(nt)); iTr = (0.0_MK, 0.0_MK)
      allocate(iTt(nt)); iTt = (0.0_MK, 0.0_MK)
      allocate(iTz(nt)); iTz = (0.0_MK, 0.0_MK)
      allocate(u(nt)); u = 0.0_MK

      allocate(t(ntRec)); t = [ ( (i - 1) * dt, i = 1, ntRec ) ]
      allocate(ur(ntRec)); ur = 0.0_MK
      allocate(ut(ntRec)); ut = 0.0_MK
      allocate(ux(ntRec)); ux = 0.0_MK
      allocate(uy(ntRec)); uy = 0.0_MK
      allocate(uz(ntRec)); uz = 0.0_MK
      allocate(tr(ntRec)); tr = 0.0_MK
      allocate(tt(ntRec)); tt = 0.0_MK
      allocate(tx(ntRec)); tx = 0.0_MK
      allocate(ty(ntRec)); ty = 0.0_MK
      allocate(tz(ntRec)); tz = 0.0_MK

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
          iUr(i) =   Sw2p * ( uIntg(1) + uIntg(2) )
          iUt(i) =   Sw2p *   uIntg(3)
          iUz(i) = - Sw2p * ( uIntg(4) + uIntg(5) )

          iTr(i) =   Sw2p * ( uIntg(6) + uIntg(7) )
          iTt(i) =   Sw2p *   uIntg(8)
          iTz(i) = - Sw2p * ( uIntg(9) + uIntg(10) )

        else
          iUr(i) =   Sw2p * ( uIntg(1 ) + uIntg(2 ) + uIntg(3 ) + uIntg(4 ) )
          iUt(i) =   Sw2p * ( uIntg(5 ) + uIntg(6 ) )
          iUz(i) = - Sw2p * ( uIntg(7 ) + uIntg(8 ) + uIntg(9 ) + uIntg(10) )

          iTr(i) =   Sw2p * ( uIntg(11) + uIntg(12) + uIntg(13) + uIntg(14) )
          iTt(i) =   Sw2p * ( uIntg(15) + uIntg(16) )
          iTz(i) = - Sw2p * ( uIntg(17) + uIntg(18) + uIntg(19) + uIntg(20) )
        end if

        ni = ni + 1
      end do
      !$OMP END DO

#ifdef PROGBAR
      if(isMain) &
        & call commProgressBar(fiEx, fiEx, 'run k-integration by frequency')
#endif

      !$OMP END PARALLEL

#ifdef FFTW
      call FFTW_(execute_dft_c2r) (pcr, iUr, u); ur = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iUt, u); ut = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iUz, u); uz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTr, u); tr = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTt, u); tt = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (pcr, iTz, u); tz = u(1:ntRec) * df * exp(feps * t)
#else
      !=> For the DFT, we have $ X_{N - m} = X_{-m} $. And specially for purely
      ! real input, the output is Hermitian-symmetric, that is, the negative-
      ! frequency terms are just the conjugates of the corresponding positive-
      ! frequency terms.
      iUr(nt:nt/2 + 2:-1) = conjg( iUr(2:nt/2) )
      iUt(nt:nt/2 + 2:-1) = conjg( iUt(2:nt/2) )
      iUz(nt:nt/2 + 2:-1) = conjg( iUz(2:nt/2) )
      iTr(nt:nt/2 + 2:-1) = conjg( iTr(2:nt/2) )
      iTt(nt:nt/2 + 2:-1) = conjg( iTt(2:nt/2) )
      iTz(nt:nt/2 + 2:-1) = conjg( iTz(2:nt/2) )
      u = real( fft(nt, iUr, -1) ); ur = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iUt, -1) ); ut = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iUz, -1) ); uz = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTr, -1) ); tr = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTt, -1) ); tt = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTz, -1) ); tz = u(1:ntRec) * df * exp(feps * t)
#endif
    end subroutine dwimRun

    subroutine dwimTransform()
      ux = ur * cos(theta) - ut * sin(theta)
      uy = ut * cos(theta) + ur * sin(theta)
      tx = tr * cos(theta) - tt * sin(theta)
      ty = tt * cos(theta) + tr * sin(theta)
    end subroutine dwimTransform

    subroutine dwimFinalize()
      if(allocated(wabs)) deallocate(wabs)

      if(allocated(iUr)) deallocate(iUr)
      if(allocated(iUt)) deallocate(iUt)
      if(allocated(iUz)) deallocate(iUz)
      if(allocated(iTr)) deallocate(iTr)
      if(allocated(iTt)) deallocate(iTt)
      if(allocated(iTz)) deallocate(iTz)
      if(allocated(u)) deallocate(u)

      if(allocated(t)) deallocate(t)
      if(allocated(ur)) deallocate(ur)
      if(allocated(ut)) deallocate(ut)
      if(allocated(ux)) deallocate(ux)
      if(allocated(uy)) deallocate(uy)
      if(allocated(uz)) deallocate(uz)
      if(allocated(tr)) deallocate(tr)
      if(allocated(tt)) deallocate(tt)
      if(allocated(tx)) deallocate(tx)
      if(allocated(ty)) deallocate(ty)
      if(allocated(tz)) deallocate(tz)

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
      real(kind = MK) :: kr, jn(0:3), dj(0:2)

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
      dj(0) = - jn(1)
      dj(1) = ( jn(0) - jn(2) ) / 2.0_MK

      if(isForce) then
        !ref.: eqs. (4-7a, b and c)
        intg(1 ) = IPS1 * ( gSH1    * k * jn(0)  + (gPS1(1) - gSH1) * k * dj(1) )
        intg(2 ) = IPS0 *   gPS0(1) * k * dj(0)
        intg(3 ) = ISH1 * ( gPS1(1) * k * jn(0)  + (gSH1 - gPS1(1)) * k * dj(1) )
        intg(4 ) = IPS1 *   gPS1(2) * k * jn(1)
        intg(5 ) = IPS0 *   gPS0(2) * k * jn(0)

        intg(6 ) = IPS1 * ( tSH1    * k * jn(0)  + (tPS1(1) - tSH1) * k * dj(1) )
        intg(7 ) = IPS0 *   tPS0(1) * k * dj(0)
        intg(8 ) = ISH1 * ( tPS1(1) * k * jn(0)  + (tSH1 - tPS1(1)) * k * dj(1) )
        intg(9 ) = IPS1 *   tPS1(2) * k * jn(1)
        intg(10) = IPS0 *   tPS0(2) * k * jn(0)

      else
        jn(3) = bessel_jn(3, kr)
        dj(2) = ( jn(1) - jn(3) ) / 2.0_MK

        !ref.: eqs. (4-9a, b and c)
        intg(1 ) =   RPS2  * ( qSH2    * k * jn(1)  + (qPS2(1) - qSH2) * k * dj(2) )
        intg(2 ) = - RPS1  * ( qSH1    * k * jn(0)  + (qPS1(1) - qSH1) * k * dj(1) )
        intg(3 ) = - RPS01 *   qPS2(1) * k * dj(0)
        intg(4 ) =   RPS02 *   qPS0(1) * k * dj(0)
        intg(5 ) =   RSH2  * ( qPS2(1) * k * jn(1)  + (qSH2 - qPS2(1)) * k * dj(2) )
        intg(6 ) =   RSH1  * ( qPS1(1) * k * jn(0)  + (qSH1 - qPS1(1)) * k * dj(1) )
        intg(7 ) =   RPS2  *   qPS2(2) * k * jn(2)
        intg(8 ) = - RPS1  *   qPS1(2) * k * jn(1)
        intg(9 ) =   RPS02 *   qPS0(2) * k * jn(0)
        intg(10) = - RPS01 *   qPS2(2) * k * jn(0)

        intg(11) =   RPS2  * ( bSH2    * k * jn(1)  + (bPS2(1) - bSH2) * k * dj(2) )
        intg(12) = - RPS1  * ( bSH1    * k * jn(0)  + (bPS1(1) - bSH1) * k * dj(1) )
        intg(13) = - RPS01 *   bPS2(1) * k * dj(0)
        intg(14) =   RPS02 *   bPS0(1) * k * dj(0)
        intg(15) =   RSH2  * ( bPS2(1) * k * jn(1)  + (bSH2 - bPS2(1)) * k * dj(2) )
        intg(16) =   RSH1  * ( bPS1(1) * k * jn(0)  + (bSH1 - bPS1(1)) * k * dj(1) )
        intg(17) =   RPS2  *   bPS2(2) * k * jn(2)
        intg(18) = - RPS1  *   bPS1(2) * k * jn(1)
        intg(19) =   RPS02 *   bPS0(2) * k * jn(0)
        intg(20) = - RPS01 *   bPS2(2) * k * jn(0)
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
