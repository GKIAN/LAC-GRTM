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
  character(len = 3) :: StrLab(6)
  !$ real(kind = MK) :: tbeg, tend

  call paraGetArguments()
  !$ tbeg = omp_get_wtime()

  call mathInitialize()
  call paraInitialize()
  call grtcInitialize()
  call dwimInitialize()

  write(*, '(A)') 'Output waveform file prefix is: ' // trim(adjustl(outPref))
  UxyzFile = trim(adjustl(outPref)) // '.gu'
  TxyzFile = trim(adjustl(outPref)) // '.gt'

  !$ if(ompNthd <= 0) ompNthd = omp_get_num_procs()
  !$ call omp_set_num_threads(ompNthd)

  call dwimRun()

  !$ tend = omp_get_wtime()
  !$ write(*, '(A, F10.6, A)') 'OpenMP elapsed time is ', tend - tbeg, ' s.'

  if(toxyz) then
#ifndef STRAIN
      StrLab = ['Txx', 'Tyy', 'Tzz', 'Txy', 'Txz', 'Tyz']
#else
      StrLab = ['Exx', 'Eyy', 'Ezz', 'Exy', 'Exz', 'Eyz']
#endif
    call dwimTransform()
    call outputDispla(UxyzFile, ux, uy, uz, ['Ux ', 'Uy ', 'Uz '])
    call outputStress(TxyzFile, txx, tyy, tzz, txy, txz, tyz, StrLab)
  else
#ifndef STRAIN
    StrLab = ['Trr', 'Ttt', 'Tzz', 'Trt', 'Trz', 'Ttz']
#else
    StrLab = ['Err', 'Ett', 'Ezz', 'Ert', 'Erz', 'Etz']
#endif
    call outputDispla(UxyzFile, ur, ut, uz, ['Ur ', 'Ut ', 'Uz '])
    call outputStress(TxyzFile, trr, ttt, tzz, trt, trz, ttz, StrLab)
  end if

  call dwimFinalize()
  call grtcFinalize()
  call paraFinalize()

  contains

    subroutine outputDispla(fileName, ua, ub, uc, labels)
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
    end subroutine outputDispla

    subroutine outputStress(fileName, ua, ub, uc, ud, ue, uf, labels)
      character(len = *), intent(in) :: fileName
      real(kind = MK), intent(in) :: ua(:), ub(:), uc(:), ud(:), ue(:), uf(:)
      character(len = *), intent(in) :: labels(:)
      integer :: fileID, i

      open(newunit = fileID, file = fileName)
        write(fileID, '(7(13X, A, 11X))') 't  ', labels
        do i = 1, ntRec
          write(fileID, '(7(2X, ES25.17E3))') t(i), ua(i), ub(i), uc(i), &
            & ud(i), ue(i), uf(i)
        end do
      close(fileID)
    end subroutine outputStress

end program

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
