      subroutine amopen(lun,fname,fstat,fform,facc)
!
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
!                Regents of the University of California               **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************
!
!  When this is converted to Fortran 90
!  these codes can be replaced with module types, eg,
!       Character(*), public, parameter :: unknown = 'unknown'
!  facc is not used.

      implicit none
!
!     INPUT:
!
      integer lun
!        ... logical unit number
      character*(*) fname
!        ... file name (not used in VAX/VMS implementation)
      character*1 fstat
!        ... status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
      character*1 fform
!        ... format code: 'F', 'U' = formatted, unformatted
      character*1 facc
!        ... access code: 'R', 'W', 'A' = read, read/write, append
!
!     INTERNAL:

      character*7 stat
!        ... status keyword
      character*11 kform
!        ... form keyword
      integer ios
!        ... i/o status variable
      integer errlun
!        ... where to log error messages (STDOUT or STDERR)

!
      integer FILE_UNIT_FIRST; parameter(FILE_UNIT_FIRST=20)
      integer FILE_UNIT_LAST; parameter(FILE_UNIT_LAST=50)
      logical in_use
!

      if (fstat == 'N') then
           stat = 'NEW'
      elseif (fstat == 'O') then
           stat = 'OLD'
      elseif (fstat == 'R') then
           stat = 'REPLACE'
      elseif (fstat == 'U') then
           stat = 'UNKNOWN'
      else
           write(stdout,'(/,2x,a,i4)')
     &           'amopen: bogus fstat, unit ', lun
           call mexit(6, 1)
      endif
!
      if (fform == 'U') then
           kform = 'UNFORMATTED'
      elseif (fform == 'F') then
           kform = 'FORMATTED'
      else
           write(stdout,'(/,2x,a,i4)')
     &           'amopen: bogus fform, unit', lun
           call mexit(6, 1)
      endif
!
! If lun is -1, find an available unit number:
      if ( lun < 0 ) then
        lun = FILE_UNIT_FIRST
        in_use = .true.
        do while (lun <= FILE_UNIT_LAST .and. in_use)
          inquire(unit=lun,opened=in_use)
        enddo
! If we have checked all of them, and none available, return error
        if (in_use) then
          write(stdout,'(A)')                                                
     &'ERROR in AMOPEN: no File LUNs available for automatic allocation'
          call mexit(6,1)       
        endif
      endif

      open(unit=lun,file=fname,status=stat,form=kform,iostat=ios)

      if (ios /= 0) then
        if (lun == 6) then
          errlun=0 !STDERR
        else
          errlun=6 !STDOUT
        endif
        write(stderr,'(/,2x,a,i4,a,a)') 'Unit ',lun,' Error on OPEN: ',fname
        !close(unit=6) !FIXME: why close STDOUT??
        call mexit(6, 1)
      endif
      rewind(lun)
      return
      end


