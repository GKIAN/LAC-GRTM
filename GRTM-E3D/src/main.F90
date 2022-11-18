!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 22 Dec 2020 10:12:11 AM CST
!-------------------------------------------------------------------------------

program main

  !$ use omp_lib
  use comm
  use math
  use paraMod
  use grtcMod
  use dwimMod
  implicit none

  character(len = LLS) :: UxyzFile, TxyzFile
  !$ real(kind = MK) :: tbeg, tend

  call paraGetArguments()
  !$ tbeg = omp_get_wtime()

  call mathInitialize()
  call paraInitialize()
  call grtcInitialize()
  call dwimInitialize()

  UxyzFile = trim(adjustl(outPref)) // '.gu'
  TxyzFile = trim(adjustl(outPref)) // '.gt'

  call dwimRun()

  !$ tend = omp_get_wtime()
  !$ write(*, '(A, F10.6, A)') 'OpenMP elapsed time is ', tend - tbeg, ' s.'

  if(toxyz) then
    call dwimTransform()
    call output(UxyzFile, ux, uy, uz, ['Ux ', 'Uy ', 'Uz '])
    call output(TxyzFile, tx, ty, tz, ['Txz', 'Tyz', 'Tzz'])
  else
    call output(UxyzFile, ur, ut, uz, ['Ur ', 'Ut ', 'Uz '])
    call output(TxyzFile, tr, tt, tz, ['Trz', 'Ttz', 'Tzz'])
  end if

  call dwimFinalize()
  call grtcFinalize()
  call paraFinalize()

  contains

    subroutine output(fileName, ua, ub, uc, labels)
      character(len = *), intent(in) :: fileName
      real(kind = MK), intent(in) :: ua(:), ub(:), uc(:)
      character(len = *), intent(in) :: labels(:)
      integer :: fileID, i

      open(newunit = fileID, file = fileName)
        write(fileID, '(4(13X, A, 11X))') 't  ', labels
        do i = 1, ntRec
          write(fileID, '(4(2X, ES25.17E3))') t(i), ua(i), ub(i), uc(i)
        end do
      close(fileID)
    end subroutine output

end program

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
