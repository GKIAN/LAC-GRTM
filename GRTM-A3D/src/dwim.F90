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
    & svInty, swType, swTime, swFreq, srTime, lSrc, lRec, z, alpha
  use grtcMod
  implicit none
  private

  integer, public :: ntRec

  integer :: nt, nf, ivMin
  real(kind = MK) :: df, dk
  real(kind = MK) :: fMax, feps, dmin
  real(kind = MK) :: r

  integer :: fiEx = -1
  real(kind = MK), allocatable :: wabs(:)

  integer, parameter :: nIntg = 3, npt = 10
  logical :: suga(nIntg, 2)
  integer :: ipts(nIntg, 2)
  real(kind = MK) :: vpts(nIntg, 2, npt), suIntg(nIntg, 2, 3)
  complex(kind = MK) :: uIntg(nIntg), duIntg(nIntg)
  !$OMP THREADPRIVATE(suga, ipts, vpts, suIntg, uIntg, duIntg)

  complex(kind = MK), allocatable :: iUr(:), iUz(:)
  complex(kind = MK), allocatable :: iTz(:)
  real(kind = MK), allocatable :: u(:)
  real(kind = MK), allocatable, public :: t(:)
  real(kind = MK), allocatable, public :: ur(:), ux(:), uy(:), uz(:)
  real(kind = MK), allocatable, public :: tz(:)

  public dwimInitialize, dwimRun, dwimTransform, dwimFinalize

  contains

    subroutine dwimInitialize()
      real(kind = MK) :: L, vMin, vMax, aLim
      complex(kind = MK) :: omg
      integer :: m, iros, iloc(1)
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

      vMin = minval(alpha)
      vMax = maxval(alpha)
      !=> space period length: it should be so long that wave won't arrive at
      ! any receiver in the time window?
      L = vMax * nt * dt + r + vMax / ((vMin + vMax) / 2.0_MK) &
        & * sqrt(r * r + (rxyz(3) - sxyz(3)) ** 2) + 100.0D+3
      dk = min(dkLim, 2.0_MK * pi / L)

      fiEx = nf
      allocate(wabs(nf)); wabs = 0.0_MK
      do i = 1, nf
        omg = CAL_OMEGA( df * (i - 1) )
        wabs(i) = abs(mathWavelet(omg, swType, swTime, swFreq, srTime))
      end do
      iloc = maxloc(wabs)
      j = iloc(1)
      aLim = wabs(j) * faLim
      do i = nf, j, -1
        if(wabs(i) >= aLim) exit
      end do
      fiEx = min(i + 1, nf)

      if(lRec < lSrc) then
        iros = min(lRec + 1, lSrc - 1)
        iloc = minloc( alpha(iros:lSrc) )
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
        iloc = minloc( alpha(lSrc:iros) )
        ivMin = iloc(1) + lSrc - 1
        dmin = min( z(lSrc) - sxyz(3), rxyz(3) - z(lSrc) )
        do j = lSrc + 1, lRec - 1
          dmin = min( dmin, z(j) - z(j - 1) )
        end do
      end if

      allocate(iUr(nt)); iUr = (0.0_MK, 0.0_MK)
      allocate(iUz(nt)); iUz = (0.0_MK, 0.0_MK)
      allocate(iTz(nt)); iTz = (0.0_MK, 0.0_MK)
      allocate(u(nt)); u = 0.0_MK

      allocate(t(ntRec)); t = [ ( (i - 1) * dt, i = 1, ntRec ) ]
      allocate(ur(ntRec)); ur = 0.0_MK
      allocate(ux(ntRec)); ux = 0.0_MK
      allocate(uy(ntRec)); uy = 0.0_MK
      allocate(uz(ntRec)); uz = 0.0_MK
      allocate(tz(ntRec)); tz = 0.0_MK

#ifdef DEBUG
      write(*, '(A, 2(1X, G0))') 'dwimInit: nt =', m, nt
      write(*, '(A, 2(1X, G0))') 'dwimInit: df =', fMax, df
      write(*, '(A, 2(1X, G0))') 'dwimInit: dk =', L, dk
#endif
    end subroutine dwimInitialize

    subroutine dwimRun()
      real(kind = MK) :: k, kCri
      complex(kind = MK) :: omg, Sw
      integer :: i, j, jj, ni
      logical :: isMain = .true.

#ifdef FFTW
      integer(kind = 8) :: p
      include 'fftw3.f'
#endif

      ni = 0
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, jj, k, omg, kCri, Sw, &
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
        Sw = svInty * mathWavelet(omg, swType, swTime, swFreq, srTime)
        kCri = kfCri * kmax(i, omg)

        ! when under the critical k-value, integrate by trapezoidal rule
        k = dk
        suIntg = 0.0_MK
        do while(.true.)
          duIntg = integrand(k, omg, Sw) * dk
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
          duIntg = integrand(k, omg, Sw) * dk
          suIntg(:, 1, 3) = suIntg(:, 1, 2) + real(duIntg)
          suIntg(:, 2, 3) = suIntg(:, 2, 2) + imag(duIntg)

          do j = 1, 2
            do jj = 1, nIntg
              suga(jj, j) = ptamGather(suIntg(jj, j, :), ipts(jj, j), &
                & vpts(jj, j, :))
            end do
          end do
          if(all(suga)) exit

          suIntg(:, :, 1) = suIntg(:, :, 2)
          suIntg(:, :, 2) = suIntg(:, :, 3)
        end do

        do j = 1, 2
          do jj = 1, nIntg
            call ptamReduce(ipts(jj, j), vpts(jj, j, :))
          end do
        end do
        uIntg = cmplx(vpts(:, 1, 1), vpts(:, 2, 1), kind = MK)

        ! assemble the spectra of displacements
        iUr(i) =   uIntg(1) / (2.0_MK * pi)
        iUz(i) = - uIntg(2) / (2.0_MK * pi)
        iTz(i) = - uIntg(3) / (2.0_MK * pi)

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
      call FFTW_(execute_dft_c2r) (p, iUz, u); uz = u(1:ntRec) * df * exp(feps * t)
      call FFTW_(execute_dft_c2r) (p, iTz, u); tz = u(1:ntRec) * df * exp(feps * t)

      call FFTW_(destroy_plan) (p)

#else
      !=> For the DFT, we have $ X_{N - m} = X_{-m} $. And specially for purely
      ! real input, the output is Hermitian-symmetric, that is, the negative-
      ! frequency terms are just the conjugates of the corresponding positive-
      ! frequency terms.
      iUr(nt:nt/2 + 2:-1) = conjg( iUr(2:nt/2) )
      iUz(nt:nt/2 + 2:-1) = conjg( iUz(2:nt/2) )
      iTz(nt:nt/2 + 2:-1) = conjg( iTz(2:nt/2) )

      u = real( fft(nt, iUr, -1) ); ur = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iUz, -1) ); uz = u(1:ntRec) * df * exp(feps * t)
      u = real( fft(nt, iTz, -1) ); tz = u(1:ntRec) * df * exp(feps * t)

#endif
    end subroutine dwimRun

    subroutine dwimTransform()
      ux = ur * cos(pi / 4.0_MK)
      uy = ur * sin(pi / 4.0_MK)
    end subroutine dwimTransform

    subroutine dwimFinalize()
      if(allocated(wabs)) deallocate(wabs)

      if(allocated(iUr)) deallocate(iUr)
      if(allocated(iUz)) deallocate(iUz)
      if(allocated(iTz)) deallocate(iTz)
      if(allocated(u)) deallocate(u)

      if(allocated(t)) deallocate(t)
      if(allocated(ur)) deallocate(ur)
      if(allocated(ux)) deallocate(ux)
      if(allocated(uy)) deallocate(uy)
      if(allocated(uz)) deallocate(uz)
      if(allocated(tz)) deallocate(tz)
    end subroutine dwimFinalize

    real(kind = MK) function kmax(fi, omg)
      integer, intent(in) :: fi
      complex(kind = MK), value :: omg
      real(kind = MK) :: ka

      if(fi <= 3) omg = CAL_OMEGA(df * 3)
      ka = abs(omg/galpha(ivMin))

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

    function integrand(k, omg, Sw) result(intg)
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg, Sw
      complex(kind = MK) :: intg(nIntg)
      real(kind = MK) :: kr, j0, dj0

      call grtcSetEigen(k, omg)
      call grtcCoefficient()
      call grtcKernelCoeffi()
      call grtcKernelMoment(k, omg, Sw)

      kr = k * r
      j0 = bessel_j0(kr)
      dj0 = - bessel_j1(kr)

      !ref.: eqs. (4-6a, b and c)
      intg(1) = UAC * k * dj0
      intg(2) = WAC * k * j0
      intg(3) = TAC * k * j0
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
