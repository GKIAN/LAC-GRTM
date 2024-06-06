!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Tue 22 Dec 2020 10:37:01 AM CST
!-------------------------------------------------------------------------------

module paraMod

  use comm
  use math
  implicit none
  public

  integer :: ompNthd
  character(len = LLS) :: inputFile = 'input.conf', modelFile, outPref

  real(kind = MK) :: tRec, dtRec, faLim
  real(kind = MK) :: kLim, kfCri, dkLim

  real(kind = MK) :: xs, zs, xr, zr

  real(kind = MK) :: svInty
  real(kind = MK) :: swTime, swFreq, srTime
  character(len = LSS) :: swType

  integer :: nLayer = 0
  integer :: lRec = 0, lSrc = 0
  real(kind = MK), allocatable :: z(:), rho(:), alpha(:), Qs(:), Qp(:)

  contains

    subroutine paraGetArguments()
      integer :: nArgs
      character(len = LLS) :: argStr = 'NULL'
      nArgs = command_argument_count()
      if(nArgs >= 1) then
        call get_command_argument(1, argStr)
        if(argStr == '-h' .or. argStr == '--help') then
          call paraPrintHelp()
          stop
        else
          inputFile = argStr
        end if
      end if
      if(nArgs > 1) call commErrorExcept(NOFATALERR, &
        & 'Some invalid arguments found.')

#ifdef DEBUG
      write(*, '(A)') 'paraGetArguments: inputFile = ' // trim(inputFile)
#endif
    end subroutine paraGetArguments

    subroutine paraPrintHelp()
      character(len = LSS) :: exeName
      call get_command_argument(0, exeName)
      write(*, *)
      write(*, '(A)') 'Usage:'
      write(*, '(2X, A)') trim(adjustl(exeName)) // ' [inputFile]'
      write(*, '(2X, A)') '             calculate 2-D acoustic waveform ' &
        & // 'for an explosive source in a multi-layer half-space.'
      write(*, *)
      write(*, '(2X, A)') 'inputFile    the input configuring file. ' &
        & // '[default: ' // trim(adjustl(inputFile)) // ']'
      write(*, *)
    end subroutine paraPrintHelp

    subroutine paraInitialize()
      integer :: fileID, ioStatus, nLine
      integer :: i, iTemp
      real(kind = MK) :: rTemp
      character(len = LSS) :: sTemp

      ! read input parameters from file
      open(newunit = fileID, file = inputFile, status = 'old')
        call commGetConf(fileID, toGoon, 'OpenMP_num_threads', 1, ompNthd)
        call commGetConf(fileID, toRwnd, 'model_file', 1, modelFile)
        call commGetConf(fileID, toRwnd, 'output_prefix', 1, outPref)
        call commGetConf(fileID, toRwnd, 'record_time_length', 1, tRec)
        call commGetConf(fileID, toRwnd, 'record_time_step', 1, dtRec)
        call commGetConf(fileID, toRwnd, 'frequency_limit_amplitude_ratio', &
          & 1, faLim)
        call commGetConf(fileID, toRwnd, 'integrate_limit_k-step', 1, dkLim)
        call commGetConf(fileID, toRwnd, 'integrate_limit_k-value', 1, kLim)
        call commGetConf(fileID, toRwnd, 'integrate_critical_k-factor', &
          & 1, kfCri)
        call commGetConf(fileID, toRwnd, 'coordinate_source', 1, xs)
        call commGetConf(fileID, toBack, 'coordinate_source', 2, zs)
        call commGetConf(fileID, toRwnd, 'coordinate_receiver', 1, xr)
        call commGetConf(fileID, toBack, 'coordinate_receiver', 2, zr)
        call commGetConf(fileID, toRwnd, 'source_vibrate_intensity', 1, svInty)
        call commGetConf(fileID, toRwnd, 'source_wavelet_type', 1, swType)
        call commGetConf(fileID, toRwnd, 'source_wavelet_time', 1, swTime)
        call commGetConf(fileID, toRwnd, 'source_wavelet_frequency', 1, swFreq)
        call commGetConf(fileID, toRwnd, 'source_rise_time', 1, srTime)
      close(fileID)
      if(ompNthd < 0) then
        call get_environment_variable('OMP_NUM_THREADS', sTemp)
        read(sTemp, *, iostat = ioStatus) ompNthd
      end if
      do while(.true.)
        i = index(outPref, '@')
        if(i == 0) exit
        select case(outPref(i + 1:i + 1))
          case('M')
            outPref = outPref(:i - 1) // trim(commStringStrip(modelFile, &
              & .true., .true.)) // outPref(i + 2:)
          case('C')
            outPref = outPref(:i - 1) // trim(commStringStrip(inputFile, &
              & .true., .true.)) // outPref(i + 2:)
          case default
            call commErrorExcept(FAIL2CHECK, &
              & 'Unresolved output prefix format <' // outPref(i:i + 1) // '>')
        end select
      end do

      ! read model parameters from file
      open(newunit = fileID, file = modelFile, status = 'old')
        !=> get number of lines of file
        nLine = - 1
        read(fileID, *, iostat = ioStatus)
        do while(ioStatus == 0)
          nLine = nLine + 1
          read(fileID, *, iostat = ioStatus) iTemp, rTemp
        end do
        !=> determine whether receiver or source is in the half space
        if( (zr > rTemp) .or. (zs > rTemp) ) then
          nLayer = nLine + 1
        else
          nLayer = nLine
        end if
        !=> allocate model parameters
        allocate(z(0:nLayer))
        allocate(rho  (nLayer))
        allocate(alpha(nLayer))
        allocate(Qs(nLayer))
        allocate(Qp(nLayer))
        !=> read model parameters
        rewind(fileID)
        read(fileID, *)
        do i = 1, nLine
          read(fileID, *) iTemp, z(i - 1), rho(i), alpha(i), Qs(i), Qp(i)
        end do
      close(fileID)
      ! insert a fictitious interface, to ensure source not in half space
      if(nLayer > nLine) then
        call commErrorExcept(NOFATALERR, 'Insert a fictitious interface.')
        z(nLine) = max(zr, zs) * 2.0_MK - z(nLine - 1)
        rho  (nLayer) = rho  (nLine)
        alpha(nLayer) = alpha(nLine)
        Qs(nLayer) = Qs(nLine)
        Qp(nLayer) = Qp(nLine)
      end if
      z(nLayer) = inf
      ! check some keypoints
      if( any( z(1:nLayer) < z(0:nLayer - 1) ) ) &
        & call commErrorExcept(FAIL2CHECK, 'The lower interface MUST be ' &
        & // 'deeper than the upper interface in any layer.')
      if( (zr < z(0)) .or. (zs < z(0)) ) &
        & call commErrorExcept(FAIL2CLOSE, 'Receiver and source MUST be ' &
        & // 'on/below the free surface.')

      ! other parameters
      if(kfCri < 0.0_MK) kfCri = 1.0_MK
      lSrc = - 1
      lRec = - 1
      do i = 1, nLayer
        if(zs >= z(i - 1) .and. zs < z(i)) lSrc = i
        if(zr >= z(i - 1) .and. zr < z(i)) lRec = i
      end do

#ifdef DEBUG
      write(*, '(A)') 'paraInitialize: modelFile = ' // trim(modelFile)
      write(*, '(A, 3(1X, G0))') 'paraInitialize: s =', xs, zs, lSrc
      write(*, '(A, 3(1X, G0))') 'paraInitialize: r =', xr, zr, lRec
#endif
    end subroutine paraInitialize

    subroutine paraFinalize()
      if(allocated(z)) deallocate(z)
      if(allocated(rho  )) deallocate(rho  )
      if(allocated(alpha)) deallocate(alpha)
      if(allocated(Qs)) deallocate(Qs)
      if(allocated(Qp)) deallocate(Qp)
    end subroutine paraFinalize

end module paraMod

! vim:ft=fortran tw=80 ts=2 sw=2 et ai 
