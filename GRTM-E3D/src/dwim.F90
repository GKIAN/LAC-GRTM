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
  use paraMod, only: tRec, faLim, kfCri, kLim, dt => dtRec, dkLim, sxyz, rxyz, &
    & fxyz, mxyz, svInty, swType, swTime, swFreq, srTime, isForce, lSrc, lRec, &
    & z, beta, alpha
  use grtcMod
  implicit none
  private

  integer, public :: ntRec

  integer :: nt, nf, ivMin
  real(kind = MK) :: df, dk
  real(kind = MK) :: fMax, feps, dmin

  integer :: nig
  real(kind = MK) :: r, theta
  real(kind = MK) :: ISH1, IPS1, IPS0
  real(kind = MK) :: RSH2, RSH1, RPS2, RPS1, RPS01, RPS02

  integer :: fiEx = -1
  real(kind = MK), allocatable :: wabs(:)

  integer, parameter :: nIntg = 16, npt = 10
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
      real(kind = MK) :: L, vMin, vMax, faMax
      complex(kind = MK) :: omg
      integer :: m, iros, iloc(1), fi1, fi2
      integer :: i, j

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

      fi1 = 1
      fi2 = nf
      allocate(wabs(fi1:fi2 + 5)); wabs = 0.0_MK
      do i = fi1, fi2 + 5
        omg = CAL_OMEGA( df * (i - 1) )
        wabs(i) = abs(mathWavelet(omg, swType, swTime, swFreq, srTime))
      end do
      iloc = maxloc(wabs)
      j = iloc(1)
      faMax = wabs(j)
      do i = j, fi2
        if(sum(wabs(i:i + 5)) / 6 < faLim * faMax) then
          fiEx = i
          exit
        end if
      end do

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
        nig = 16
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

#ifdef FFTW
      integer(kind = 8) :: p
      include 'fftw3.f'
#endif

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
        Sw2p = svInty * mathWavelet(omg, swType, swTime, swFreq, srTime)
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
              & 'Integrate to the limit k-value. Result NOT reliable.')
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
          iUr(i) =   Sw2p * ( uIntg(1 ) + uIntg(2 ) + uIntg(3 ) )
          iUt(i) =   Sw2p * ( uIntg(4 ) + uIntg(5 ) )
          iUz(i) = - Sw2p * ( uIntg(6 ) + uIntg(7 ) + uIntg(8 ) )

          iTr(i) =   Sw2p * ( uIntg(9 ) + uIntg(10) + uIntg(11) )
          iTt(i) =   Sw2p * ( uIntg(12) + uIntg(13) )
          iTz(i) = - Sw2p * ( uIntg(14) + uIntg(15) + uIntg(16) )
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
      call FFTW_(plan_dft_c2r_1d) (p, nt, iUr, u, FFTW_ESTIMATE)

      call FFTW_(execute_dft_c2r) (p, iUr, u); ur = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iUt, u); ut = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iUz, u); uz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iTr, u); tr = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iTt, u); tt = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iTz, u); tz = u(1:ntRec) * df * exp(feps * t)

      call FFTW_(destroy_plan) (p)

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
      real(kind = MK) :: kr, j0, j1, j2, dj0, dj1, dj2

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
      j0 = bessel_j0(kr)
      j1 = bessel_j1(kr)
      j2 = bessel_jn(2, kr)
      dj0 = - j1
      dj1 = (j0 - j2) / 2.0_MK

      if(isForce) then
        !ref.: eqs. (4-7a, b and c)
        intg(1 ) = IPS1 * ( gSH1    * k * j0  + (gPS1(1) - gSH1) * k * dj1 )
        intg(2 ) = IPS0 *   gPS0(1) * k * dj0
        intg(3 ) = ISH1 * ( gPS1(1) * k * j0  + (gSH1 - gPS1(1)) * k * dj1 )
        intg(4 ) = IPS1 *   gPS1(2) * k * j1
        intg(5 ) = IPS0 *   gPS0(2) * k * j0

        intg(6 ) = IPS1 * ( tSH1    * k * j0  + (tPS1(1) - tSH1) * k * dj1 )
        intg(7 ) = IPS0 *   tPS0(1) * k * dj0
        intg(8 ) = ISH1 * ( tPS1(1) * k * j0  + (tSH1 - tPS1(1)) * k * dj1 )
        intg(9 ) = IPS1 *   tPS1(2) * k * j1
        intg(10) = IPS0 *   tPS0(2) * k * j0

      else
        dj2 = ( j1 - bessel_jn(3, kr) ) / 2.0_MK

        !ref.: eqs. (4-9a, b and c)
        intg(1 ) =   RPS2  * ( qSH2    * k * j1  + (qPS2(1) - qSH2) * k * dj2 )
        intg(2 ) = - RPS1  * ( qSH1    * k * j0  + (qPS1(1) - qSH1) * k * dj1 )
        intg(3 ) = - RPS01 *   qPS2(1) * k * dj0 + RPS02  * qPS0(1) * k * dj0
        intg(4 ) =   RSH2  * ( qPS2(1) * k * j1  + (qSH2 - qPS2(1)) * k * dj2 )
        intg(5 ) =   RSH1  * ( qPS1(1) * k * j0  + (qSH1 - qPS1(1)) * k * dj1 )
        intg(6 ) =   RPS2  *   qPS2(2) * k * j2
        intg(7 ) = - RPS1  *   qPS1(2) * k * j1
        intg(8 ) =   RPS02 *   qPS0(2) * k * j0  - RPS01  * qPS2(2) * k * j0

        intg(9 ) =   RPS2  * ( bSH2    * k * j1  + (bPS2(1) - bSH2) * k * dj2 )
        intg(10) = - RPS1  * ( bSH1    * k * j0  + (bPS1(1) - bSH1) * k * dj1 )
        intg(11) = - RPS01 *   bPS2(1) * k * dj0 + RPS02  * bPS0(1) * k * dj0
        intg(12) =   RSH2  * ( bPS2(1) * k * j1  + (bSH2 - bPS2(1)) * k * dj2 )
        intg(13) =   RSH1  * ( bPS1(1) * k * j0  + (bSH1 - bPS1(1)) * k * dj1 )
        intg(14) =   RPS2  *   bPS2(2) * k * j2
        intg(15) = - RPS1  *   bPS1(2) * k * j1
        intg(16) =   RPS02 *   bPS0(2) * k * j0  - RPS01  * bPS2(2) * k * j0
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
      do i = ipt - 1, 1, - 1
        do j = 1, i, 1
          vpt(j) = (vpt(j) + vpt(j + 1)) / 2.0_MK
        end do
      end do
    end subroutine ptamReduce

end module dwimMod

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
