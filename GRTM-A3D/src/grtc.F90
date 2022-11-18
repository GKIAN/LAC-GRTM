!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Wed 23 Dec 2020 09:27:16 AM CST
!-------------------------------------------------------------------------------

module grtcMod

  use comm
  use math
  use paraMod, only: s => lSrc, j => lRec, sxyz, rxyz, nLayer, &
    & z, rho, alpha, Qp, Qs
  implicit none
  private

  integer :: N
  logical :: isUp
  real(kind = MK) :: zs, zr
  complex(kind = MK), allocatable, public :: galpha(:)
  complex(kind = MK), allocatable :: nu(:)
  !$OMP THREADPRIVATE(galpha, nu)

  complex(kind = MK), allocatable :: E(:, :, :)
  complex(kind = MK), allocatable :: gTu(:), gTd(:), gRud(:), gRdu(:)
  complex(kind = MK) :: Eus, Eds
  complex(kind = MK) :: YAC, ZAC
  !$OMP THREADPRIVATE(E, gTu, gTd, gRud, gRdu, Eus, Eds, YAC, ZAC)
#ifdef PeiDHS
  complex(kind = MK), allocatable :: iE(:, :, :)
  !$OMP THREADPRIVATE(iE)
#endif
  complex(kind = MK), public :: WAC, TAC, UAC
  !$OMP THREADPRIVATE(WAC, TAC, UAC)

  public grtcInitialize, grtcSetMedia, grtcSetEigen, grtcFinalize
  public grtcCoefficient, grtcKernelCoeffi, grtcKernelMoment

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
      allocate(galpha(nLayer)); galpha = (0.0_MK, 0.0_MK)
      allocate(nu(nLayer)); nu = (0.0_MK, 0.0_MK)
      allocate(E(2, 2, nLayer)); E = (0.0_MK, 0.0_MK)
#ifdef PeiDHS
      allocate(iE(2, 2, nLayer)); iE = (0.0_MK, 0.0_MK)
#endif

      allocate(gTu (1:s - 1)); gTu  = (0.0_MK, 0.0_MK)
      allocate(gRud(0:s - 1)); gRud = (0.0_MK, 0.0_MK)
      allocate(gTd (s:N));     gTd  = (0.0_MK, 0.0_MK)
      allocate(gRdu(s:N + 1)); gRdu = (0.0_MK, 0.0_MK)
      !$OMP END PARALLEL
    end subroutine grtcInitialize

    subroutine grtcSetMedia(omg, omgMax)
      complex(kind = MK), intent(in) :: omg
      real(kind = MK), intent(in) :: omgMax
      !ref.: (5.88) in P.182 of (Aki and Richards, 1980)
      galpha = alpha * ( 1.0_MK + log(omg / omgMax) / (pi * Qp) &
        & + (0.0_MK, 1.0_MK) / (2.0_MK * Qp) )
    end subroutine grtcSetMedia

    subroutine grtcSetEigen(k, omg)
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg

      nu = sqrt( k * k - (omg / galpha) ** 2 )

      E(1, :, :) = (1.0_MK, 0.0_MK)
      E(2, 1, :) = rho * omg * omg / nu
      E(2, 2, :) = - E(2, 1, :)
#ifdef PeiDHS
      iE(:, 1, :) = (0.5_MK, 0.0_MK)
      iE(1, 2, :) = 0.5_MK * nu / (rho * omg * omg)
      iE(2, 2, :) = - iE(1, 2, :)
#endif
    end subroutine grtcSetEigen

    subroutine grtcCoefficient()
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
      gRud(0) = - E(2, 2, 1) / E(2, 1, 1) * Lu1z0
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
        E1(:, 1) =   E(:, 1, i + 1)
        E1(:, 2) = - E(:, 2, i)
        E2(:, 1) =   E(:, 1, i)
        E2(:, 2) = - E(:, 2, i + 1)
        X = matmul( matmul( MatInv22(E1), E2 ), La )
        gTu(i) = X(2, 2) / (1.0_MK - X(2, 1) * gRud(i - 1))
        gRud(i) = X(1, 2) + X(1, 1) * gRud(i - 1) * gTu(i)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE(:, :, i + 1), E(:, :, i))
        gTu(i) = La(2, 2) / ( X(2, 1) * La(1, 1) * gRud(i - 1) + X(2, 2) )
        gRud(i) = ( X(1, 1) * La(1, 1) * gRud(i - 1) + X(1, 2) ) * gTu(i)
#endif
      end do

      !=> for the last layer (N + 1):
      !ref.: eqs. (3-31) and (3-33)
      gRdu(N + 1) = (0.0_MK, 0.0_MK)
      !=> for j = s, s + 1, ..., N:
      do i = N, s, - 1
        !ref.: eqs. (3-20a), (3-20b) and (3-48)
        if(i /= s) then
          La(1, 1) = exp( - nu(i) * (z(i) - z(i - 1)) )
        else
          La(1, 1) = exp( - nu(i) * (z(i) - zs) )
        end if
        La(2, 2) = exp( - nu(i + 1) * (z(i + 1) - z(i)) )
#ifndef PeiDHS
        !ref.: eqs. (3-33) and (3-28a)
        E1(:, 1) =   E(:, 1, i + 1)
        E1(:, 2) = - E(:, 2, i)
        E2(:, 1) =   E(:, 1, i)
        E2(:, 2) = - E(:, 2, i + 1)
        X = matmul( matmul( MatInv22(E1), E2 ), La )
        gTd(i) = X(1, 1) / (1.0_MK - X(1, 2) * gRdu(i + 1))
        gRdu(i) = X(2, 1) + X(2, 2) * gRdu(i + 1) * gTd(i)
#else
        !ref.: eq. (11) in (Pei, 2008, doi: 10.1785/0120070057)
        X = matmul(iE(:, :, i), E(:, :, i + 1))
        gTd(i) = La(1, 1) / ( X(1, 1) + X(1, 2) * La(2, 2) * gRdu(i + 1) )
        gRdu(i) = ( X(2, 1) + X(2, 2) * La(2, 2) * gRdu(i + 1) ) * gTd(i)
#endif
      end do
    end subroutine grtcCoefficient

    subroutine grtcKernelCoeffi()
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
        YAC = E(1, 1, j) * Ldj * gRud(j - 1) + E(1, 2, j) * Luj
        ZAC = E(2, 1, j) * Ldj * gRud(j - 1) + E(2, 2, j) * Luj
        Ms = Eus * gRdu(s) * Eds * gRud(s - 1)
      else
        YAC = E(1, 1, j) * Ldj + E(1, 2, j) * Luj * gRdu(j)
        ZAC = E(2, 1, j) * Ldj + E(2, 2, j) * Luj * gRdu(j)
        Ms = Eds * gRud(s - 1) * Eus * gRdu(s)
      end if
      Mj = 1.0_MK / (1.0_MK - Ms)
      if(isUp) then
        do i = s - 1, j, - 1
          Mj = gTu(i) * Mj
        end do
      else
        do i = s, j - 1, + 1
          Mj = gTd(i) * Mj
        end do
      end if
      YAC = YAC * Mj
      ZAC = ZAC * Mj
    end subroutine grtcKernelCoeffi

    subroutine grtcKernelMoment(k, omg, Mw)
      ! Kernel for moment-tensor source
      real(kind = MK), intent(in) :: k
      complex(kind = MK), intent(in) :: omg, Mw
      complex(kind = MK) :: gSu, gSd, lam, Ms, krw

      lam = rho(s) * galpha(s) * galpha(s)
      gSu = Mw / (2.0_MK * lam)
      gSd = - gSu

      !ref.: eqs. (3-51) or (3-51a)
      if(isUp) then
        Ms = gSu + Eus * gRdu(s) * gSd
      else
        Ms = gSd + Eds * gRud(s - 1) * gSu
      end if
      WAC = YAC * Ms
      TAC = ZAC * Ms

      krw = k / ( rho(j) * omg * omg )
      UAC = krw * TAC
      if(zr == zs) UAC = UAC + krw * Mw
    end subroutine grtcKernelMoment

    subroutine grtcFinalize()
      !$OMP PARALLEL
      if(allocated(galpha)) deallocate(galpha)
      if(allocated(nu)) deallocate(nu)
      if(allocated(E)) deallocate(E)
#ifdef PeiDHS
      if(allocated(iE)) deallocate(iE)
#endif

      if(allocated(gTu )) deallocate(gTu )
      if(allocated(gTd )) deallocate(gTd )
      if(allocated(gRud)) deallocate(gRud)
      if(allocated(gRdu)) deallocate(gRdu)
      !$OMP END PARALLEL
    end subroutine grtcFinalize

end module grtcMod

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
