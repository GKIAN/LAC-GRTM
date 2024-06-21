!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Wed 23 Dec 2020 09:27:16 AM CST
!-------------------------------------------------------------------------------

module grtcMod

  use comm
  use math
  use paraMod, only: omgR => waRef, s => lSrc, j => lRec, sxyz, rxyz, nLayer, &
    & z, rho, alpha, beta, Qp, Qs
  implicit none
  private

  integer :: N
  logical :: isUp
  real(kind = MK) :: zs, zr
  complex(kind = MK), allocatable, public :: galpha(:), gbeta(:)
  complex(kind = MK), allocatable :: mu(:)
  complex(kind = MK), allocatable :: nu(:)
  !$OMP THREADPRIVATE(galpha, gbeta, mu, nu)

#ifdef SH
  complex(kind = MK), allocatable :: E22(:, :, :)
  complex(kind = MK), allocatable :: gTu(:), gTd(:), gRud(:), gRdu(:)
  complex(kind = MK) :: Eus, Eds
  complex(kind = MK) :: YSH, ZSH
  !$OMP THREADPRIVATE(E22, gTu, gTd, gRud, gRdu, Eus, Eds, YSH, ZSH)
#ifdef PeiDHS
  complex(kind = MK), allocatable :: iE22(:, :, :)
  !$OMP THREADPRIVATE(iE22)
#endif
#endif
  complex(kind = MK), public :: gSH1, qSH1, qSH2
  complex(kind = MK), public :: tSH1, bSH1, bSH2
  !$OMP THREADPRIVATE(gSH1, qSH1, qSH2, tSH1, bSH1, bSH2)

#ifdef PS
  complex(kind = MK), allocatable :: chi(:), gam(:)
  complex(kind = MK), allocatable :: E44(:, :, :)
  complex(kind = MK), allocatable :: gTu22(:, :, :), gTd22(:, :, :), &
    & gRud22(:, :, :), gRdu22(:, :, :)
  complex(kind = MK) :: Eus22(2, 2), Eds22(2, 2)
  complex(kind = MK) :: YPS(2, 2), ZPS(2, 2)
  !$OMP THREADPRIVATE(chi, gam, E44, gTu22, gTd22, gRud22, gRdu22, &
  !$OMP   & Eus22, Eds22, YPS, ZPS)
#ifdef PeiDHS
  complex(kind = MK), allocatable :: iE44(:, :, :)
  !$OMP THREADPRIVATE(iE44)
#endif
#endif
  complex(kind = MK), public :: gPS0(2), gPS1(2), qPS0(2), qPS1(2), qPS2(2)
  complex(kind = MK), public :: tPS0(2), tPS1(2), bPS0(2), bPS1(2), bPS2(2)
  !$OMP THREADPRIVATE(gPS0, gPS1, qPS0, qPS1, qPS2, tPS0, tPS1, &
  !$OMP   & bPS0, bPS1, bPS2)

  public grtcInitialize, grtcSetMedia, grtcSetEigen, grtcFinalize
#ifdef SH
  public grtcCoefficientSH, grtcKernelCoeffiSH, grtcKernelForceSH, &
    & grtcKernelMomentSH
#endif
#ifdef PS
  public grtcCoefficientPS, grtcKernelCoeffiPS, grtcKernelForcePS, &
    & grtcKernelMomentPS
#endif

  contains

    subroutine grtcInitialize()
      N = nLayer - 1
      zs = sxyz(3)
      zr = rxyz(3)
      if(zr <= zs) then
        isUp = .true.
      else
        isUp = .false.
      end if

      !$OMP PARALLEL
      Eus22 = (0.0_MK, 0.0_MK)
      Eds22 = (0.0_MK, 0.0_MK)

      allocate(gbeta (nLayer)); gbeta  = (0.0_MK, 0.0_MK)
      allocate(galpha(nLayer)); galpha = (0.0_MK, 0.0_MK)
      allocate(mu(nLayer)); mu = (0.0_MK, 0.0_MK)

      allocate(nu(nLayer)); nu = (0.0_MK, 0.0_MK)
#ifdef SH
      allocate(E22(2, 2, nLayer)); E22 = (0.0_MK, 0.0_MK)
#ifdef PeiDHS
      allocate(iE22(2, 2, nLayer)); iE22 = (0.0_MK, 0.0_MK)
#endif

      allocate(gTu (1:s - 1)); gTu  = (0.0_MK, 0.0_MK)
      allocate(gRud(0:s - 1)); gRud = (0.0_MK, 0.0_MK)
      allocate(gTd (s:N));     gTd  = (0.0_MK, 0.0_MK)
      allocate(gRdu(s:N + 1)); gRdu = (0.0_MK, 0.0_MK)
#endif

#ifdef PS
      allocate(chi(nLayer)); chi = (0.0_MK, 0.0_MK)
      allocate(gam(nLayer)); gam = (0.0_MK, 0.0_MK)
      allocate(E44(4, 4, nLayer)); E44 = (0.0_MK, 0.0_MK)
#ifdef PeiDHS
      allocate(iE44(4, 4, nLayer)); iE44 = (0.0_MK, 0.0_MK)
#endif

      allocate(gTu22 (2, 2, 1:s - 1)); gTu22  = (0.0_MK, 0.0_MK)
      allocate(gRud22(2, 2, 0:s - 1)); gRud22 = (0.0_MK, 0.0_MK)
      allocate(gTd22 (2, 2, s:N));     gTd22  = (0.0_MK, 0.0_MK)
      allocate(gRdu22(2, 2, s:N + 1)); gRdu22 = (0.0_MK, 0.0_MK)
#endif

      gSH1 = (0.0_MK, 0.0_MK)
      qSH1 = (0.0_MK, 0.0_MK)
      qSH2 = (0.0_MK, 0.0_MK)

      tSH1 = (0.0_MK, 0.0_MK)
      bSH1 = (0.0_MK, 0.0_MK)
      bSH2 = (0.0_MK, 0.0_MK)

      gPS0 = (0.0_MK, 0.0_MK)
      gPS1 = (0.0_MK, 0.0_MK)
      qPS0 = (0.0_MK, 0.0_MK)
      qPS1 = (0.0_MK, 0.0_MK)
      qPS2 = (0.0_MK, 0.0_MK)

      tPS0 = (0.0_MK, 0.0_MK)
      tPS1 = (0.0_MK, 0.0_MK)
      bPS0 = (0.0_MK, 0.0_MK)
      bPS1 = (0.0_MK, 0.0_MK)
      bPS2 = (0.0_MK, 0.0_MK)
      !$OMP END PARALLEL
    end subroutine grtcInitialize

    subroutine grtcSetMedia(omg)
      complex(kind = MK), intent(in) :: omg

      !ref.: (5.88/5.94) in P.182/175 of (Aki and Richards, 1980/2002)
      gbeta  = beta  * ( 1.0_MK + log(omg / omgR) / (pi * Qs) &
        & + (0.0_MK, 1.0_MK) / (2.0_MK * Qs) )
      galpha = alpha * ( 1.0_MK + log(omg / omgR) / (pi * Qp) &
        & + (0.0_MK, 1.0_MK) / (2.0_MK * Qp) )
      mu = rho * gbeta * gbeta
    end subroutine grtcSetMedia

    subroutine grtcSetEigen(k, omg)
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg

#ifdef PeiDHS
      complex(kind = MK) :: tau(nLayer)
#endif

      !ref.: eq. (3-16)
      nu = sqrt(k * k - (omg / gbeta) ** 2)
#ifdef SH
      !ref: eq. (3-14)
      E22(1, :, :) = (1.0_MK, 0.0_MK)
      E22(2, 1, :) = - mu * nu
      E22(2, 2, :) = - E22(2, 1, :)
#ifdef PeiDHS
      iE22(:, 1, :) = (0.5_MK, 0.0_MK)
      iE22(1, 2, :) = - 0.5_MK / (mu * nu)
      iE22(2, 2, :) = - iE22(1, 2, :)
#endif
#endif

#ifdef PS
      !ref.: eqs. (3-18), (3-18b) and (3-18c)
      chi = k * k + nu * nu
      gam = sqrt(k * k - (omg / galpha) ** 2)
      E44(1, 1, :) = galpha * k
      E44(2, 1, :) = galpha * gam
      E44(3, 1, :) = - 2.0_MK * galpha * mu * k * gam
      E44(4, 1, :) = - galpha * mu * chi
      E44(1, 2, :) = gbeta * nu
      E44(2, 2, :) = gbeta * k
      E44(3, 2, :) = - gbeta * mu * chi
      E44(4, 2, :) = - 2.0_MK * gbeta * mu * k * nu
      E44(1, 3:4, :) =   E44(1, 1:2, :)
      E44(2, 3:4, :) = - E44(2, 1:2, :)
      E44(3, 3:4, :) = - E44(3, 1:2, :)
      E44(4, 3:4, :) =   E44(4, 1:2, :)
      E44 = E44 / omg
#ifdef PeiDHS
      tau =  gbeta / (2.0_MK * galpha * mu * gam * nu * omg)
      iE44(1, 1, :) =   tau * 2.0_MK * gbeta * mu * k * gam * nu
      iE44(1, 2, :) = - tau * gbeta * mu * nu * chi
      iE44(1, 3, :) = - tau * gbeta * k * nu
      iE44(1, 4, :) =   tau * gbeta * gam * nu
      iE44(2, 1, :) = - tau * galpha * mu * gam * chi
      iE44(2, 2, :) =   tau * 2.0_MK * galpha * mu * k * gam * nu
      iE44(2, 3, :) =   tau * galpha * gam * nu
      iE44(2, 4, :) = - tau * galpha * k * gam
      iE44(3:4, 1, :) =   iE44(1:2, 1, :)
      iE44(3:4, 2, :) = - iE44(1:2, 2, :)
      iE44(3:4, 3, :) = - iE44(1:2, 3, :)
      iE44(3:4, 4, :) =   iE44(1:2, 4, :)
#endif
#endif
    end subroutine grtcSetEigen

#ifdef SH
    subroutine grtcCoefficientSH()
      complex(kind = MK) :: Lu1z0, La(2, 2) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: X(2, 2)
      !$OMP THREADPRIVATE(La)

#ifndef PeiDHS
      complex(kind = MK) :: E1(2, 2), E2(2, 2)
#endif

      integer i

      !=> for the free surface:
      !ref.: eqs. (3-20b) and (3-30)
      if(s /= 1) then
        Lu1z0 = exp( - nu(1) * (z(1) - z(0)) )
      else
        Lu1z0 = exp( - nu(1) * (zs - z(0)) )
      end if
      gRud(0) = - E22(2, 2, 1) / E22(2, 1, 1) * Lu1z0
      !=> for j = 1, 2, ..., s - 1:
      do i = 1, s - 1
        !ref.: eqs. (3-20a), (3-20b) and (3-48)
        La(1, 1) = exp( - nu(i) * (z(i) - z(i - 1)) )
        if(i /= s - 1) then
          La(2, 2) = exp( - nu(i + 1) * (z(i + 1) - z(i)) )
        else
          La(2, 2) = exp( - nu(i + 1) * (zs - z(i)) )
        end if
#ifndef PeiDHS
        !ref.: eqs. (3-32) and (3-28a)
        E1(:, 1) =   E22(:, 1, i + 1)
        E1(:, 2) = - E22(:, 2, i)
        E2(:, 1) =   E22(:, 1, i)
        E2(:, 2) = - E22(:, 2, i + 1)
        X = matmul( matmul( MatInv22(E1), E2 ), La )
        gTu(i) = X(2, 2) / (1.0_MK - X(2, 1) * gRud(i - 1))
        gRud(i) = X(1, 2) + X(1, 1) * gRud(i - 1) * gTu(i)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE22(:, :, i + 1), E22(:, :, i))
        gTu(i) = La(2, 2) / ( X(2, 1) * La(1, 1) * gRud(i - 1) + X(2, 2) )
        gRud(i) = ( X(1, 1) * La(1, 1) * gRud(i - 1) + X(1, 2) ) * gTu(i)
#endif
      end do

      !=> for the last layer (N + 1):
      !ref.: eqs. (3-31) and (3-33)
      gRdu(N + 1) = (0.0_MK, 0.0_MK)
      !=> for j = s, s + 1, ..., N:
      do i = N, s, -1
        !ref.: eqs. (3-20a), (3-20b) and (3-48)
        if(i /= s) then
          La(1, 1) = exp( - nu(i) * (z(i) - z(i - 1)) )
        else
          La(1, 1) = exp( - nu(i) * (z(i) - zs) )
        end if
        La(2, 2) = exp( - nu(i + 1) * (z(i + 1) - z(i)) )
#ifndef PeiDHS
        !ref.: eqs. (3-33) and (3-28a)
        E1(:, 1) =   E22(:, 1, i + 1)
        E1(:, 2) = - E22(:, 2, i)
        E2(:, 1) =   E22(:, 1, i)
        E2(:, 2) = - E22(:, 2, i + 1)
        X = matmul( matmul( MatInv22(E1), E2 ), La )
        gTd(i) = X(1, 1) / (1.0_MK - X(1, 2) * gRdu(i + 1))
        gRdu(i) = X(2, 1) + X(2, 2) * gRdu(i + 1) * gTd(i)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE22(:, :, i), E22(:, :, i + 1))
        gTd(i) = La(1, 1) / ( X(1, 1) + X(1, 2) * La(2, 2) * gRdu(i + 1) )
        gRdu(i) = ( X(2, 1) + X(2, 2) * La(2, 2) * gRdu(i + 1) ) * gTd(i)
#endif
      end do
    end subroutine grtcCoefficientSH

    subroutine grtcKernelCoeffiSH()
      complex(kind = MK) :: Ldj, Luj
      complex(kind = MK) :: Mj, Ms
      integer :: i

      !ref.: eqs. (3-47b) and (3-47a)
      Eus = exp( - nu(s) * (z(s) - zs) )
      Eds = exp( - nu(s) * (zs - z(s - 1)) )
      !ref.: eqs. (3-20b) & (3-20a) or (3-46c) & (3-46d)
      if(j == s .and. isUp) then
        Luj = exp( - nu(s) * (zs - zr) )
      else
        Luj = exp( - nu(j) * (z(j) - zr) )
      end if
      if(j == s .and. (.not. isUp)) then
        Ldj = exp( - nu(s) * (zr - zs) )
      else
        Ldj = exp( - nu(j) * (zr - z(j - 1)) )
      end if

      !ref.: eqs. (3-52a) - (3-52d)
      if(isUp) then
        YSH = E22(1, 1, j) * Ldj * gRud(j - 1) + E22(1, 2, j) * Luj
        ZSH = E22(2, 1, j) * Ldj * gRud(j - 1) + E22(2, 2, j) * Luj
        Ms = Eus * gRdu(s) * Eds * gRud(s - 1)
      else
        YSH = E22(1, 1, j) * Ldj + E22(1, 2, j) * Luj * gRdu(j)
        ZSH = E22(2, 1, j) * Ldj + E22(2, 2, j) * Luj * gRdu(j)
        Ms = Eds * gRud(s - 1) * Eus * gRdu(s)
      end if
      Mj = 1.0_MK / (1.0_MK - Ms)
      if(isUp) then
        do i = s - 1, j, -1
          Mj = gTu(i) * Mj
        end do
      else
        do i = s, j - 1, +1
          Mj = gTd(i) * Mj
        end do
      end if
      YSH = YSH * Mj
      ZSH = ZSH * Mj
    end subroutine grtcKernelCoeffiSH

    subroutine grtcKernelForceSH()
      ! Kernel of SH-component for single-force source
      complex(kind = MK) :: fSH, PSH
      !ref.: eq. (3-62a)
      fSH = 1.0_MK / (2.0_MK * mu(s) * nu(s))
      !ref.: eqs. (3-66a) or (3-66d)
      if(isUp) then
        PSH = Eus * gRdu(s) + 1.0_MK
      else
        PSH = 1.0_MK + Eds * gRud(s - 1)
      end if
      gSH1 = YSH * PSH * fSH
      tSH1 = ZSH * PSH * fSH
    end subroutine grtcKernelForceSH

    subroutine grtcKernelMomentSH(k)
      ! Kernel of SH-component for moment-tensor source
      real(kind = MK), intent(in) :: k
      complex(kind = MK) :: sSH, P1, P2
      !ref.: eq. (3-64a)
      sSH = 1.0_MK / (2.0_MK * mu(s) * nu(s))
      !ref.: eqs. (3-68b) and (3-68a)
      if(isUp) then
        P1 = Eus * gRdu(s)
        P2 = 1.0_MK + P1
        P1 = 1.0_MK - P1
      else
        P1 = Eds * gRud(s - 1)
        P2 = P1 + 1.0_MK
        P1 = P1 - 1.0_MK
      end if
      qSH1 = YSH * P1 * nu(s) * sSH
      bSH1 = ZSH * P1 * nu(s) * sSH
      qSH2 = YSH * P2 * k * sSH
      bSH2 = ZSH * P2 * k * sSH
    end subroutine grtcKernelMomentSH
#endif

#ifdef PS
    subroutine grtcCoefficientPS()
      complex(kind = MK) :: Lu1z0(2, 2) = (0.0_MK, 0.0_MK), &
        & La(4, 4) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: M(2, 2), X(4, 4)
      !$OMP THREADPRIVATE(Lu1z0, La)

#ifndef PeiDHS
      complex(kind = MK) :: E1(4, 4), E2(4, 4)
#endif

      integer i

      !=> for the free surface:
      !ref.: eqs. (3-24c) and (3-40a)
      if(s /= 1) then
        Lu1z0(1, 1) = exp( - gam(1) * (z(1) - z(0)) )
        Lu1z0(2, 2) = exp( -  nu(1) * (z(1) - z(0)) )
      else
        Lu1z0(1, 1) = exp( - gam(1) * (zs - z(0)) )
        Lu1z0(2, 2) = exp( -  nu(1) * (zs - z(0)) )
      end if
      M = MatInv22(E44(3:4, 1:2, 1))
      gRud22(:, :, 0) = - matmul( matmul(M, E44(3:4, 3:4, 1)), Lu1z0)
      !=> for j = 1, 2, ..., s - 1:
      do i = 1, s - 1
        !ref.: eqs. (3-24b), (3-24c) and (3-48)
        La(1, 1) = exp( - gam(i) * (z(i) - z(i - 1)) )
        La(2, 2) = exp( -  nu(i) * (z(i) - z(i - 1)) )
        if(i /= s - 1) then
          La(3, 3) = exp( - gam(i + 1) * (z(i + 1) - z(i)) )
          La(4, 4) = exp( -  nu(i + 1) * (z(i + 1) - z(i)) )
        else
          La(3, 3) = exp( - gam(i + 1) * (zs - z(i)) )
          La(4, 4) = exp( -  nu(i + 1) * (zs - z(i)) )
        end if
#ifndef PeiDHS
        !ref.: eqs. (3-40a) and (3-35)
        E1(:, 1:2) =   E44(:, 1:2, i + 1)
        E1(:, 3:4) = - E44(:, 3:4, i)
        E2(:, 1:2) =   E44(:, 1:2, i)
        E2(:, 3:4) = - E44(:, 3:4, i + 1)
        X = matmul( matmul( MatInv44(E1), E2 ), La )
        M = MatInv22( I22 - matmul(X(3:4, 1:2), gRud22(:, :, i - 1)) )
        gTu22(:, :, i) = matmul(M, X(3:4, 3:4))
        M = matmul(gRud22(:, :, i - 1), gTu22(:, :, i))
        gRud22(:, :, i) = X(1:2, 3:4) + matmul(X(1:2, 1:2), M)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE44(:, :, i + 1), E44(:, :, i))
        M = matmul(X(3:4, 1:2), La(1:2, 1:2))
        M = MatInv22( matmul(M, gRud22(:, :, i - 1)) + X(3:4, 3:4) )
        gTu22(:, :, i) = matmul(M, La(3:4, 3:4))
        M = matmul(X(1:2, 1:2), La(1:2, 1:2))
        M = matmul(M, gRud22(:, :, i - 1)) + X(1:2, 3:4)
        gRud22(:, :, i) = matmul(M, gTu22(:, :, i))
#endif
      end do

      !=> for the last layer (N + 1):
      !ref.: eq. (3-40b)
      gRdu22(:, :, N + 1) = (0.0_MK, 0.0_MK)
      !=> for j = s, s + 1, ..., N:
      do i = N, s, -1
        !ref.: eqs. (3-24b), (3-24c) and (3-48)
        if(i /= s) then
          La(1, 1) = exp( - gam(i) * (z(i) - z(i - 1)) )
          La(2, 2) = exp( -  nu(i) * (z(i) - z(i - 1)) )
        else
          La(1, 1) = exp( - gam(i) * (z(i) - zs) )
          La(2, 2) = exp( -  nu(i) * (z(i) - zs) )
        end if
        La(3, 3) = exp( - gam(i + 1) * (z(i + 1) - z(i)) )
        La(4, 4) = exp( -  nu(i + 1) * (z(i + 1) - z(i)) )
#ifndef PeiDHS
        !ref.: eqs. (3-40b) and (3-35)
        E1(:, 1:2) =   E44(:, 1:2, i + 1)
        E1(:, 3:4) = - E44(:, 3:4, i)
        E2(:, 1:2) =   E44(:, 1:2, i)
        E2(:, 3:4) = - E44(:, 3:4, i + 1)
        X = matmul( matmul( MatInv44(E1), E2 ), La )
        M = MatInv22( I22 - matmul(X(1:2, 3:4), gRdu22(:, :, i + 1)) )
        gTd22(:, :, i) = matmul(M, X(1:2, 1:2))
        M = matmul(gRdu22(:, :, i + 1), gTd22(:, :, i))
        gRdu22(:, :, i) = X(3:4, 1:2) + matmul(X(3:4, 3:4), M)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE44(:, :, i), E44(:, :, i + 1))
        M = matmul(X(1:2, 3:4), La(3:4, 3:4))
        M = MatInv22( X(1:2, 1:2) + matmul(M, gRdu22(:, :, i + 1)) )
        gTd22(:, :, i) = matmul(M, La(1:2, 1:2))
        M = matmul(X(3:4, 3:4), La(3:4, 3:4))
        M = X(3:4, 1:2) + matmul(M, gRdu22(:, :, i + 1))
        gRdu22(:, :, i) = matmul(M, gTd22(:, :, i))
#endif
      end do
    end subroutine grtcCoefficientPS

    subroutine grtcKernelCoeffiPS()
      complex(kind = MK) :: Ldj(2, 2) = (0.0_MK, 0.0_MK), &
        & Luj(2, 2) = (0.0_MK, 0.0_MK)
      complex(kind = MK) :: Mj(2, 2), Ms(2, 2)
      !$OMP THREADPRIVATE(Ldj, Luj)

      integer :: i

      !ref.: eqs. (3-56e) and (3-56f)
      Eus22(1, 1) = exp( - gam(s) * (z(s) - zs) )
      Eus22(2, 2) = exp( -  nu(s) * (z(s) - zs) )
      Eds22(1, 1) = exp( - gam(s) * (zs - z(s - 1)) )
      Eds22(2, 2) = exp( -  nu(s) * (zs - z(s - 1)) )
      !ref.: eqs. (3-24c) & (3-24b) or (3-56c) & (3-56d)
      if(j == s .and. isUp) then
        Luj(1, 1) = exp( - gam(s) * (zs - zr) )
        Luj(2, 2) = exp( -  nu(s) * (zs - zr) )
      else
        Luj(1, 1) = exp( - gam(j) * (z(j) - zr) )
        Luj(2, 2) = exp( -  nu(j) * (z(j) - zr) )
      end if
      if(j == s .and. (.not. isUp)) then
        Ldj(1, 1) = exp( - gam(s) * (zr - zs) )
        Ldj(2, 2) = exp( -  nu(s) * (zr - zs) )
      else
        Ldj(1, 1) = exp( - gam(j) * (zr - z(j - 1)) )
        Ldj(2, 2) = exp( -  nu(j) * (zr - z(j - 1)) )
      end if

      !ref.: eqs. (3-58a) - (3-58d)
      if(isUp) then
        YPS = matmul( matmul(E44(1:2, 1:2, j), Ldj), gRud22(:, :, j - 1) ) &
          & + matmul(E44(1:2, 3:4, j), Luj)
        ZPS = matmul( matmul(E44(3:4, 1:2, j), Ldj), gRud22(:, :, j - 1) ) &
          & + matmul(E44(3:4, 3:4, j), Luj)
        Ms = matmul( matmul( matmul( Eus22, gRdu22(:, :, s) ), &
          & Eds22 ), gRud22(:, :, s - 1) )
      else
        YPS = matmul(E44(1:2, 1:2, j), Ldj) &
          & + matmul( matmul(E44(1:2, 3:4, j), Luj), gRdu22(:, :, j) )
        ZPS = matmul(E44(3:4, 1:2, j), Ldj) &
          & + matmul( matmul(E44(3:4, 3:4, j), Luj), gRdu22(:, :, j) )
        Ms = matmul( matmul( matmul( Eds22, gRud22(:, :, s - 1) ), &
          & Eus22 ), gRdu22(:, :, s) )
      end if
      Mj = MatInv22(I22 - Ms)
      if(isUp) then
        do i = s - 1, j, -1
          Mj = matmul(gTu22(:, :, i), Mj)
        end do
      else
        do i = s, j - 1, +1
          Mj = matmul(gTd22(:, :, i), Mj)
        end do
      end if
      YPS = matmul( YPS, Mj )
      ZPS = matmul( ZPS, Mj )
    end subroutine grtcKernelCoeffiPS

    subroutine grtcKernelForcePS(k, omg)
      ! Kernel of PSV-component for single-force source
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg
      complex(kind = MK) :: f0(2), f1(2)
      complex(kind = MK) :: P0(2, 2), P1(2, 2)
      !ref.: eqs. (3-66c) & (3-66b) or (3-66f) & (3-66e), and (3-62c) & (3-62b)
      if(isUp) then
        P0 = matmul(Eus22, gRdu22(:, :, s))
        P1 = P0 + I22
        P0 = P0 - I22
      else
        P0 = matmul(Eds22, gRud22(:, :, s - 1))
        P1 = I22 + P0
        P0 = I22 - P0
      end if

      ! $ \mu = \rho * \beta^2 $
      f0(1) = 0.5_MK / (rho(s) * galpha(s) * omg)
      f0(2) = - k / (2.0_MK * rho(s) * gbeta(s) * omg * nu(s))
      gPS0 = matmul( matmul(YPS, P0), f0 )
      tPS0 = matmul( matmul(ZPS, P0), f0 )

      ! $ \mu = \rho * \beta^2 $
      f1(1) = k / (2.0_MK * rho(s) * galpha(s) * omg * gam(s))
      f1(2) = - 0.5_MK / (rho(s) * gbeta(s) * omg)
      gPS1 = matmul( matmul(YPS, P1), f1 )
      tPS1 = matmul( matmul(ZPS, P1), f1 )
    end subroutine grtcKernelForcePS

    subroutine grtcKernelMomentPS(k, omg)
      ! Kernel of PSV-component for monent-tensor source
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg
      complex(kind = MK) :: s0(2), s1(2), s2(2)
      complex(kind = MK) :: P0(2, 2), P1(2, 2)
      !ref.: eqs. (3-68e) - (3-68c), or (3-68j) - (3-68h), and (3-64d) - (3-64b)
      if(isUp) then
        P0 = matmul(Eus22, gRdu22(:, :, s))
        P1 = P0 - I22
        P0 = P0 + I22
      else
        P0 = matmul(Eds22, gRud22(:, :, s - 1))
        P1 = I22 - P0
        P0 = I22 + P0
      end if

      ! $ \mu = \rho * \beta^2 $
      s0(1) = gam(s) / (2.0_MK * rho(s) * galpha(s) * omg)
      s0(2) = - k / (2.0_MK * rho(s) * gbeta(s) * omg)
      qPS0 = matmul( matmul(YPS, P0), s0 )
      bPS0 = matmul( matmul(ZPS, P0), s0 )

      ! $ \mu = \rho * \beta^2 $ and $ \chi = k^2 + \nu^2 $
      s1(1) = - k / (rho(s) * galpha(s) * omg)
      s1(2) = chi(s) / (2.0_MK * rho(s) * gbeta(s) * omg * nu(s))
      qPS1 = matmul( matmul(YPS, P1), s1 )
      bPS1 = matmul( matmul(ZPS, P1), s1 )

      ! $ \mu = \rho * \beta^2 $
      s2(1) = k * k / (2.0_MK * rho(s) * galpha(s) * omg * gam(s))
      s2(2) = - k / (2.0_MK * rho(s) * gbeta(s) * omg)
      qPS2 = matmul( matmul(YPS, P0), s2 )
      bPS2 = matmul( matmul(ZPS, P0), s2 )
    end subroutine grtcKernelMomentPS
#endif

    subroutine grtcFinalize()
      !$OMP PARALLEL
      if(allocated(gbeta )) deallocate(gbeta )
      if(allocated(galpha)) deallocate(galpha)
      if(allocated(mu)) deallocate(mu)

      if(allocated(nu)) deallocate(nu)
#ifdef SH
      if(allocated(E22)) deallocate(E22)
#ifdef PeiDHS
      if(allocated(iE22)) deallocate(iE22)
#endif

      if(allocated(gTu )) deallocate(gTu )
      if(allocated(gTd )) deallocate(gTd )
      if(allocated(gRud)) deallocate(gRud)
      if(allocated(gRdu)) deallocate(gRdu)
#endif

#ifdef PS
      if(allocated(chi)) deallocate(chi)
      if(allocated(gam)) deallocate(gam)
      if(allocated(E44)) deallocate(E44)
#ifdef PeiDHS
      if(allocated(iE44)) deallocate(iE44)
#endif

      if(allocated(gTu22 )) deallocate(gTu22 )
      if(allocated(gTd22 )) deallocate(gTd22 )
      if(allocated(gRud22)) deallocate(gRud22)
      if(allocated(gRdu22)) deallocate(gRdu22)
#endif
      !$OMP END PARALLEL
    end subroutine grtcFinalize

end module grtcMod

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
