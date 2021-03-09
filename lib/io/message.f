!:: OPT_FLAGS=-O1
#include "io.h"

!!
!! NOTE: These are implemented as macros, so this file must
!!       be included with a cpp #include, and must be cpp'ed at
!!       compile time; automatic for source with the ".F" extension.
!!       Or, for g77 use the flag "-xf77-cpp-input"
!! ALSO NOTE: Macros can extend line length, so fixed lengths should
!!       be disabled. g77: "-ffixed-line-length-none"
!!               IRIX f77: "-extend_source"
!!
!! MSG(type,format) ... ! This macro updates the log-buffer pointer,
!!                      ! with syntax like a standard write staement, but
!!                      ! the first argument is a log class. Negative
!!                      ! values are errors, positive are standard information.
!!
!!     Values for type:
!!       INFO_CORE == 0      Always printed.
!!       INFO_STANDARD == 1  Information for each big cycle (rebuild non-bond list, etc.)
!!       INFO_VERBOSE == 2   Information for each small cycle.
!!       INFO_DETAILS == 3   Too much info for normal runs, but sometimes helpful.
!!       INFO_DEBUG == 4     Loads of stuff, even dumping some arrays, etc.
!!      
!!       WARN_CAUTION == -1  Something looks unusual.
!!       WARN_SEVERE == -2   I think it's wrong, but I may be able to continue.
!!       WARN_FATAL == -3    Cannot continue. Bye-bye.
!!
!! INFO(format) ...     ! Shortcut for MSG(0,format)
!!
!! call Error()         ! Call error handler. Dumps trace and timing information
!!                      ! then exits.
!! call ErrorMsg(text)  ! Same as Error(), but with text. 
!----------------------------------------------------
      subroutine Message_Init(log_verbosity,warn_max)
      implicit none
      integer log_verbosity, warn_max
      include "message.inc"
      integer i

      p_log_data=Buffer_Alloc(
     &    0,LOG_BUFF_SIZE,LOG_LINE_LENGTH,'$LogBuff')
      p_log_verbosity=DB_Add(p_log_data,DB_INT,'Verbosity',1)
      IH(p_log_verbosity)=log_verbosity
      p_log_warnmax=DB_Add(p_log_data,DB_INT,'WarnMax',1)
      IH(p_log_warnmax)=warn_max

      end 

!----------------------------------------------------
      function Message_Newline(level)
      integer level
      include 'message.inc'

      if (level.gt.IH(p_log_verbosity)) then
        Message_Newline=.false.

      else
        Message_Newline=.true.
        log_ptr=log_ptr+
     $      Len_Trim(SH(log_ptr:log_ptr+LOG_LINE_LENGTH))
        SH(log_ptr:log_ptr)='\n'
        log_ptr=log_ptr+1

        if (log_ptr+LOG_LINE_LENGTH .ge.
     &      p_log_buf+LOG_BUFF_SIZE) then
          call Message_Flush()
        endif
        log_eol=log_ptr+LOG_LINE_LENGTH-1
      endif
!     write(0,*)'WARN line=SH(',log_ptr-p_log_buf,
!    &    log_eol-p_log_buf,')'
      end

!----------------------------------------------------
      subroutine Message_Clear()
      include "message.inc"
      integer i, s1,s2

      !do i=0,(LOG_BUFF_SIZE/10)-1
      !  s1=p_log_buf+i*10
      !  s2=s1+LOG_BUFF_SIZE
      !  s2=max(s2,p_log_buf+LOG_BUFF_SIZE-1)
      !  write(SH(s1:s2),'("[--",I4,"--]")'),i
      !enddo

      SH(p_log_buf:p_log_buf+LOG_BUFF_SIZE-1)=' '

      log_ptr=p_log_buf
      log_eol=p_log_buf+LOG_LINE_LENGTH-1
      end

!----------------------------------------------------
      subroutine Message_Flush()
      include "message.inc"
      integer iunit; parameter(iunit=STDERR)

      !write(iunit,'(A)') '========== START WARN BUFFER ============'
      write(iunit,'(A)') SH(p_log_buf:p_log_buf+
     &    Len_Trim(SH(p_log_buf:p_log_buf+LOG_BUFF_SIZE-1))-1)
      !write(iunit,'(A)') '========== END OF WARN BUFFER ==========='
      !FIXME: call Warn_Clear()
      log_ptr=p_log_buf
      end

!----------------------------------------------------
      subroutine ErrorMsg(text)
      implicit none
      character*(*) text
      include "message.inc"
      MSG(WARN_FATAL,'(A)') text
      call Error()
      end
!----------------------------------------------------
      subroutine Error()
      implicit none
      include "message.inc"

      ERRMSG('(36("="),A,36("="),A)') ' ERROR ','\n'
      call Timer_Dump()
      ERRMSG('(A)') '\n'
      call Trace_Dump()
      ERRMSG('(A)') '\n'
      call Mem_DumpStats()
      ERRMSG('(A)') '\n'

      call Message_flush()

      STOP

      end

!--------------------------------------------------------------------
      ! write a string to a file unit
      ! Maybe we could use negative 'units' to indicate socket streams.
      subroutine Write_Text(iunit,string)
      integer iunit
      character*(*) string
      write (iunit,'(A)') string
      end
!--------------------------------------------------------------------
