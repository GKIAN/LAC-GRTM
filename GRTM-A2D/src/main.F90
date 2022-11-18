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
  integer :: fileID, i
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

  open(newunit = fileID, file = UxyzFile)
    write(fileID, '(12X, A, 2(21X, A))') 't (s)', 'Ux (m)', 'Uz (m)'
    do i = 1, ntRec
      write(fileID, '(3(2X, ES25.17E3))') t(i), ux(i), uz(i)
    end do
  close(fileID)

  open(newunit = fileID, file = TxyzFile)
    write(fileID, '(12X, A, 21X, A)' ) 't (s)', 'Tzz (Pa)'
    do i = 1, ntRec
      write(fileID, '(2(2X, ES25.17E3))') t(i), tzz(i)
    end do
  close(fileID)

  call dwimFinalize()
  call grtcFinalize()
  call paraFinalize()

end program

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
