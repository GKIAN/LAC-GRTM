!===============================================================================
! Author: Tche L., USTC, seistche@gmail.com
! Created at: Sat 26 Jun 2021 09:27:05 PM CST
!-------------------------------------------------------------------------------

module comm
  use iso_fortran_env, only: output_unit, error_unit
  implicit none
  public

#ifndef DOUBLETYPE
  integer, parameter :: MK = 4
#else
  integer, parameter :: MK = 8
#endif

  integer, parameter :: LXSS = 4, LSS = 16, LMS = 64, LLS = 256, LXLS = 1024

  integer, parameter :: NOFATALERR = 0
  integer, parameter :: FAIL2CHECK = 1
  integer, parameter :: FAIL2FILE = 10, FAIL2OPEN = 11, FAIL2CLOSE = 12, &
    & FAIL2READ = 13, FAIL2WRITE = 14

  character(len = 1) :: charDirConntr = '/'
  character(len = 1), private :: charCommentor = '#'
  character(len = 5), private :: charSeparator = '=|,<>'

  integer, parameter :: toRwnd = - 1, toBack = 0, toGoon = 1

  integer, parameter :: LenPB = 50
  character(len = LenPB), private :: charPBar = ''

  interface commGetConf
    module procedure commGetConf_string
    module procedure commGetConf_integer
    module procedure commGetConf_real
    module procedure commGetConf_logical
  end interface commGetConf

  contains

    subroutine commErrorExcept(errorCode, suppInfo)
      integer, intent(in) :: errorCode
      character(len = *), intent(in) :: suppInfo
      if(errorCode == NOFATALERR) then
#ifndef COLORPRINT
        write(error_unit, '(A)') repeat('-', len(suppInfo) + 13)
        write(error_unit, '(A)') '>>> Warning: ' // suppInfo
        write(error_unit, '(A)') repeat('-', len(suppInfo) + 13)
#else
        write(error_unit, '(A)') char(27) // '[00;93;100mWarning:' &
          & // char(27) // '[0m ' // suppInfo
#endif
      else
#ifndef COLORPRINT
        write(error_unit, '(A)') repeat('=', len(suppInfo) + 14)
        write(error_unit, '(A, A)') '>>>>>> ERROR: ', suppInfo
        write(error_unit, '(A)') repeat('=', len(suppInfo) + 14)
#else
        write(error_unit, '(A)') char(27) // '[01;91;100mERROR:' &
          & // char(27) // '[04;39;49m ' // suppInfo // char(27) // '[0m'
#endif
        stop errorCode
      end if
    end subroutine commErrorExcept

    character(len = 12) function commFormatTime(inSecond) result(time)
      real(kind = MK), intent(in) :: inSecond
      integer :: numHour, numMinute, numSecond, leftPoint
      if( inSecond < 3600 * 1000.0_MK ) then
        numSecond = int(inSecond)
        leftPoint = int( (inSecond - numSecond) * 100 )
        numHour = numSecond / 3600
        numMinute = ( numSecond - numHour * 3600 ) / 60
        numSecond = numSecond - numHour * 3600 - numMinute * 60
        write(time, '(I3.3, 3(A, I2.2))') numHour, ':', numMinute, ':', &
          & numSecond, '.', leftPoint
      else
        time = ' >= 1000 hr '
      end if
    end function commFormatTime

    character(len = LLS) function commStringStrip(inStr, offSuff, offPath) &
      & result(outStr)
      character(len = *), intent(in) :: inStr
      logical, intent(in), optional :: offSuff, offPath
      integer :: iPath, iSuff
      if(len_trim(adjustl(inStr)) > LLS) call commErrorExcept(NOFATALERR, &
        & 'Too long string <' // trim(adjustl(inStr)) // '> for stripping.')
      outStr = adjustl(inStr)
      iPath = index(outStr, charDirConntr, .true.)
      if(present(offSuff)) then
        if(offSuff) then
          iSuff = index(outStr, '.', .true.)
          if(iSuff > iPath + 1) then
            outStr = outStr(:iSuff - 1)
          else if(iSuff == iPath + 1) then
            ! if a filename starts by '.' and has no any suffix, remove this '.'
            outStr = outStr(:iPath) // outStr(iSuff + 1:)
          end if
        end if
      end if
      if(present(offPath)) then
        if(offPath) then
          outStr = outStr(iPath + 1:)
        end if
      end if
    end function commStringStrip

    subroutine commProgressBar(iItr, nItr, suppInfo)
      integer, intent(in) :: iItr, nItr
      character(len = *), intent(in), optional :: suppInfo
      integer :: iPBar, i
      real(kind = MK) :: pBar
      character(len = 13) :: strPref
      pBar = iItr * 1.0_MK / nItr
      iPBar = int(pBar * LenPB) + 1
      do i = 1, iPBar - 1
        charPBar(i:i) = '='
      end do
      if(iPBar < LenPB) charPBar(iPBar:iPBar) = '>'
      write(strPref, '(A6, F6.2, A1)') 'DONE: ', pBar * 100.0_MK, '%'
      if(present(suppInfo)) then
        write(unit = output_unit, fmt = '(A)', advance = 'no') &
          & char(13) // strPref // ' [' // charPBar // '] in <' &
          & // trim(adjustl(suppInfo)) // '>'
      else
        write(unit = output_unit, fmt = '(A)', advance = 'no') &
          & char(13) // strPref // ' [' // charPBar // ']'
      end if
      if(iPBar > LenPB) then
        write(unit = output_unit, fmt = *)
      else
        flush(unit = output_unit)
      end if
    end subroutine commProgressBar

    !------------------- commGetConf -------------------------------------------
    subroutine commGetConf_string(fileID, firstDo, keyword, pivot, theValue)
      integer, intent(in) :: fileID, firstDo, pivot
      character(len = *), intent(in) :: keyword
      character(len = *), intent(out) :: theValue
      character(len = LLS) :: line
      character(len = LMS) :: word
      integer :: stat, indx, i

      stat = 0
      if(firstDo < 0) then
        rewind(fileID)
      else if(firstDo == 0) then
        backspace(fileID)
      end if

      do while(.true.)
        read(fileID, '(A256)', iostat = stat) line
        if(stat /= 0) exit

        line = adjustl(line)
        indx = index(line, charCommentor)
        if(indx > 0) line(indx:) = ' '
        do i = 1, len_trim(line)
          if(index(charSeparator, line(i:i)) > 0) line(i:i) = ' '
        end do
        if(len_trim(line) == 0) cycle

        read(line, *) word
        if(word == keyword) then
          ! the two numbers of heading blanks must be identical
          read(line, *, iostat = stat) (theValue, i = 0, pivot)
          exit
        end if
      end do

      if(stat /= 0) then
        write(word, *) pivot
        call commErrorExcept(FAIL2READ, 'There is no the <' &
          & // trim(adjustl(word)) // '>-th pivot in the <' &
          & // keyword // '> key-line to read.')
      end if
    end subroutine commGetConf_string

    subroutine commGetConf_integer(fileID, firstDo, keyword, pivot, theValue)
      integer, intent(in) :: fileID, firstDo, pivot
      character(len = *), intent(in) :: keyword
      integer, intent(out) :: theValue
      character(len = LMS) :: strValue
      call commGetConf_string(fileID, firstDo, keyword, pivot, strValue)
      read(strValue, *) theValue
    end subroutine commGetConf_integer

    subroutine commGetConf_real(fileID, firstDo, keyword, pivot, theValue)
      integer, intent(in) :: fileID, firstDo, pivot
      character(len = *), intent(in) :: keyword
      real(kind = MK), intent(out) :: theValue
      character(len = LMS) :: strValue
      call commGetConf_string(fileID, firstDo, keyword, pivot, strValue)
      read(strValue, *) theValue
    end subroutine commGetConf_real

    subroutine commGetConf_logical(fileID, firstDo, keyword, pivot, theValue)
      integer, intent(in) :: fileID, firstDo, pivot
      character(len = *), intent(in) :: keyword
      logical, intent(out) :: theValue
      character(len = LMS) :: strValue
      call commGetConf_string(fileID, firstDo, keyword, pivot, strValue)
      read(strValue, *) theValue
    end subroutine commGetConf_logical
    !---------------------------------------------------------------------------

    subroutine commSetConfPoint(fileID, anchor)
      integer, intent(in) :: fileID
      character(len = *), intent(in) :: anchor
      character(len = LMS) :: strValue
      call commGetConf_string(fileID, toRwnd, anchor, 0, strValue)
    end subroutine commSetConfPoint

end module comm

! vim:ft=fortran tw=80 ts=4 sw=2 et ai
