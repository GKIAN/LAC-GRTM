!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Sun 27 Dec 2020 03:05:53 PM CST
!-------------------------------------------------------------------------------

module math
  
#ifdef IEEE
  use, intrinsic :: ieee_arithmetic
#endif
  use comm
  implicit none
  public

  real(kind = MK) :: inf

  real(kind = MK), parameter :: pi = 3.141592653589793238_MK
  complex(kind = MK), parameter :: I22(2, 2) = reshape([ &
    & (1.0_MK, 0.0_MK),  (0.0_MK, 0.0_MK), &
    & (0.0_MK, 0.0_MK),  (1.0_MK, 0.0_MK)  &
    & ], [2, 2])

  interface MatDet22
    module procedure rMatDet22
    module procedure cMatDet22
  end interface MatDet22

  interface MatDet33
    module procedure rMatDet33
    module procedure cMatDet33
  end interface MatDet33

  interface MatDet44
    module procedure rMatDet44
    module procedure cMatDet44
  end interface MatDet44

  interface MatDet
    module procedure rMatDet
    module procedure cMatDet
  end interface MatDet

  interface MatInv22
    module procedure rMatInv22
    module procedure cMatInv22
  end interface MatInv22

  interface MatInv33
    module procedure rMatInv33
    module procedure cMatInv33
  end interface MatInv33

  interface MatInv44
    module procedure rMatInv44
    module procedure cMatInv44
  end interface MatInv44

  interface MatInv
    module procedure rMatInv
    module procedure cMatInv
  end interface MatInv

  contains

    subroutine mathInitialize()
      real(kind = MK) :: zero
      zero = 0.0_MK
#ifndef IEEE
      inf = 1.0_MK / zero
#else
      if(ieee_support_inf(inf)) then
        inf = ieee_value(0.0_MK, ieee_positive_inf)
      else
        inf = 1.0_MK / zero
      end if
#endif
    end subroutine mathInitialize

    complex(kind = MK) function mathWavelet(omg, wType, wTime, wFreq) result(s)
      complex(kind = MK), intent(in) :: omg
      real(kind = MK), intent(in) :: wTime, wFreq
      character(len = *), intent(in) :: wType
      real(kind = MK) :: omgc, omga, omgt
      omgc = 2.0_MK * pi * wFreq
      select case(trim(adjustl(wType)))
        case('Ricker')
          s = (omg / omgc) ** 2
          s = s * exp( - s) / wFreq * 2.0_MK / sqrt(pi)
          s = s * exp( - omg * wTime * (0.0_MK, 1.0_MK) )
        case('Green')
          s = (1.0_MK, 0.0_MK)
        case('green')
          s = omg * (0.0_MK, 1.0_MK)
        case('GReen')
          omga = abs(omg)
          omgt = 2.0_MK * pi / wTime
          if(omga < omgc + omgt) then
            if(omga > omgc) then
              s = 0.5_MK * ( 1.0_MK + cos( (omga - omgc) / omgt * pi ) )
            else
              s = 1.0_MK
            end if
          else
            s = 0.0_MK
          end if
        case default
          call commErrorExcept(FAIL2CHECK, 'Unrecognized wavelet type <' &
            & // trim(adjustl(wType)) // '>.')
      end select
    end function mathWavelet

    !=================== Determinant of Matrix =================================
    real(kind = MK) function rMatDet22(A) result(det)
      real(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
    end function rMatDet22

    real(kind = MK) function rMatDet33(A) result(det)
      real(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
        & + A(1, 2) * (A(2, 3) * A(3, 1) - A(2, 1) * A(3, 3)) &
        & + A(1, 3) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
    end function rMatDet33

    real(kind = MK) function rMatDet44(A) result(det)
      real(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * MatDet33(A(2:4, 2:4)) &
        & - A(1, 2) * MatDet33(reshape([ A(2:4, 1), A(2:4, 3:4) ], [3, 3])) &
        & + A(1, 3) * MatDet33(reshape([ A(2:4, 1:2), A(2:4, 4) ], [3, 3])) &
        & - A(1, 4) * MatDet33(A(2:4, 1:3))
    end function rMatDet44

    real(kind = MK) recursive function rMatDet(A, n) result(det)
      integer, intent(in) :: n
      real(kind = MK), intent(in) :: A(:, :)
      integer i
      if(n > 4) then
        det = 0.0_MK
        do i = 1, n
          det = det + (-1) ** (i + 1) * A(1, i) &
            & * rMatDet(reshape([ A(2:n, 1:i - 1), A(2:n, i + 1:n) ], &
            & [n - 1, n - 1]), n - 1)
        end do
      else
        if(n == 4) then
          det = rMatDet44(A)
        else if(n == 3) then
          det = rMatDet33(A)
        else if(n == 2) then
          det = rMatDet22(A)
        else if(n == 1) then
          det = A(1, 1)
        else
          det = 0.0_MK
        end if
      end if
    end function rMatDet

    complex(kind = MK) function cMatDet22(A) result(det)
      complex(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
    end function cMatDet22

    complex(kind = MK) function cMatDet33(A) result(det)
      complex(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
        & + A(1, 2) * (A(2, 3) * A(3, 1) - A(2, 1) * A(3, 3)) &
        & + A(1, 3) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
    end function cMatDet33

    complex(kind = MK) function cMatDet44(A) result(det)
      complex(kind = MK), intent(in) :: A(:, :)
      det = A(1, 1) * MatDet33(A(2:4, 2:4)) &
        & - A(1, 2) * MatDet33(reshape([ A(2:4, 1), A(2:4, 3:4) ], [3, 3])) &
        & + A(1, 3) * MatDet33(reshape([ A(2:4, 1:2), A(2:4, 4) ], [3, 3])) &
        & - A(1, 4) * MatDet33(A(2:4, 1:3))
    end function cMatDet44

    complex(kind = MK) recursive function cMatDet(A, n) result(det)
      integer, intent(in) :: n
      complex(kind = MK), intent(in) :: A(:, :)
      integer i
      if(n > 4) then
        det = (0.0_MK, 0.0_MK)
        do i = 1, n
          det = det + (-1) ** (i + 1) * A(1, i) &
            & * cMatDet(reshape([ A(2:n, 1:i - 1), A(2:n, i + 1:n) ], &
            & [n - 1, n - 1]), n - 1)
        end do
      else
        if(n == 4) then
          det = cMatDet44(A)
        else if(n == 3) then
          det = cMatDet33(A)
        else if(n == 2) then
          det = cMatDet22(A)
        else if(n == 1) then
          det = A(1, 1)
        else
          det = (0.0_MK, 0.0_MK)
        end if
      end if
    end function cMatDet
    !------------------- Determinant of Matrix ---------------------------------

    !=================== Inverse of Matrix =====================================
    function rMatInv22(A) result(iA)
      real(kind = MK), intent(in) :: A(:, :)
      real(kind = MK) :: iA(2, 2)
      iA(1, 1) =   A(2, 2)
      iA(1, 2) = - A(1, 2)
      iA(2, 1) = - A(2, 1)
      iA(2, 2) =   A(1, 1)
      iA = iA / rMatDet22(A)
    end function rMatInv22

    function rMatInv33(A) result(iA)
      real(kind = MK), intent(in) :: A(:, :)
      real(kind = MK) :: iA(3, 3)
      iA(1, 1) =     A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)
      iA(1, 2) = - ( A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2) )
      iA(1, 3) =     A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
      iA(2, 1) = - ( A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1) )
      iA(2, 2) =     A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1)
      iA(2, 3) = - ( A(1, 1) * A(2, 3) - A(1, 3) * A(2, 1) )
      iA(3, 1) =     A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1)
      iA(3, 2) = - ( A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1) )
      iA(3, 3) =     A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
      iA = iA / rMatDet33(A)
    end function rMatInv33

    function rMatInv44(A) result(iA)
      real(kind = MK), intent(in) :: A(:, :)
      real(kind = MK) :: iA(4, 4), M(3, 3)
      integer i, j
      do i = 1, 4
        do j = 1, 4
          M(1:i - 1, 1:j - 1) = A(1:i - 1, 1:j - 1)
          M(1:i - 1, j:3) = A(1:i - 1, j + 1:4)
          M(i:3, 1:j - 1) = A(i + 1:4, 1:j - 1)
          M(i:3, j:3) = A(i + 1:4, j + 1:4)
          iA(j, i) = (-1) ** (i + j) * rMatDet33(M)
        end do
      end do
      iA = iA / rMatDet44(A)
    end function rMatInv44

    function rMatInv(A, n) result(iA)
      integer, intent(in) :: n
      real(kind = MK), intent(in) :: A(:, :)
      real(kind = MK) :: iA(n, n), M(n - 1, n - 1)
      integer i, j
      do i = 1, n
        do j = 1, n
          M(1:i - 1, 1:j - 1) = A(1:i - 1, 1:j - 1)
          M(1:i - 1, j:n - 1) = A(1:i - 1, j + 1:n)
          M(i:n - 1, 1:j - 1) = A(i + 1:n, 1:j - 1)
          M(i:n - 1, j:n - 1) = A(i + 1:n, j + 1:n)
          iA(j, i) = (-1) ** (i + j) * rMatDet(M, n - 1)
        end do
      end do
      iA = iA / rMatDet(A, n)
    end function rMatInv

    function cMatInv22(A) result(iA)
      complex(kind = MK), intent(in) :: A(:, :)
      complex(kind = MK) :: iA(2, 2)
      iA(1, 1) =   A(2, 2)
      iA(1, 2) = - A(1, 2)
      iA(2, 1) = - A(2, 1)
      iA(2, 2) =   A(1, 1)
      iA = iA / cMatDet22(A)
    end function cMatInv22

    function cMatInv33(A) result(iA)
      complex(kind = MK), intent(in) :: A(:, :)
      complex(kind = MK) :: iA(3, 3)
      iA(1, 1) =     A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)
      iA(1, 2) = - ( A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2) )
      iA(1, 3) =     A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
      iA(2, 1) = - ( A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1) )
      iA(2, 2) =     A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1)
      iA(2, 3) = - ( A(1, 1) * A(2, 3) - A(1, 3) * A(2, 1) )
      iA(3, 1) =     A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1)
      iA(3, 2) = - ( A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1) )
      iA(3, 3) =     A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
      iA = iA / cMatDet33(A)
    end function cMatInv33

    function cMatInv44(A) result(iA)
      complex(kind = MK), intent(in) :: A(:, :)
      complex(kind = MK) :: iA(4, 4), M(3, 3)
      integer i, j
      do i = 1, 4
        do j = 1, 4
          M(1:i - 1, 1:j - 1) = A(1:i - 1, 1:j - 1)
          M(1:i - 1, j:3) = A(1:i - 1, j + 1:4)
          M(i:3, 1:j - 1) = A(i + 1:4, 1:j - 1)
          M(i:3, j:3) = A(i + 1:4, j + 1:4)
          iA(j, i) = (-1) ** (i + j) * cMatDet33(M)
        end do
      end do
      iA = iA / cMatDet44(A)
    end function cMatInv44

    function cMatInv(A, n) result(iA)
      integer, intent(in) :: n
      complex(kind = MK), intent(in) :: A(:, :)
      complex(kind = MK) :: iA(n, n), M(n - 1, n - 1)
      integer i, j
      do i = 1, n
        do j = 1, n
          M(1:i - 1, 1:j - 1) = A(1:i - 1, 1:j - 1)
          M(1:i - 1, j:n - 1) = A(1:i - 1, j + 1:n)
          M(i:n - 1, 1:j - 1) = A(i + 1:n, 1:j - 1)
          M(i:n - 1, j:n - 1) = A(i + 1:n, j + 1:n)
          iA(j, i) = (-1) ** (i + j) * cMatDet(M, n - 1)
        end do
      end do
      iA = iA / cMatDet(A, n)
    end function cMatInv
    !------------------- Inverse of Matrix -------------------------------------

    recursive function fft(N, x, fob) result(y)
      integer, intent(in) :: N, fob
      complex(kind = MK), intent(in) :: x(0:N - 1)
      complex(kind = MK) :: y(0:N - 1)
      complex(kind = MK) :: g(N/2), h(N/2), W(N/2)
      integer :: l
      if (N > 2) then
        W = [ ( exp( - (0.0_MK, 2.0_MK) * pi * l * fob / N ), l = 0, N/2 - 1 ) ]
        g =   x(0:N/2 - 1) + x(N/2:N - 1)
        h = ( x(0:N/2 - 1) - x(N/2:N - 1) ) * W
        y(0:N - 1:2) = fft(N/2, g, fob)
        y(1:N - 1:2) = fft(N/2, h, fob)
      else
        y(0) = x(0) + x(1)
        y(1) = x(0) - x(1)
      end if
    end function fft

end module math

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
