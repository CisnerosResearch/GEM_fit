! General I/O routines. Mostly for opening standard FILE UNITs, but may add
! support for buffered internal files for logging over a network socket, etc.
!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Open a File, and return a UNIT number. Similar to C fopen().
!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function file_open here]
FUNCTION File_Open(filename, mode)
  Implicit none
  Character(*) filename !Input, Name of file to open.
  Character(*) mode     !Input, Method of file access.
  !
  !   The first letter of mode (case insignificant) must be one of:
  !     'r'   open for reading  (status='OLD', action='READ')
  !     'w'   overwrite or create for writing
  !              (status='UNKNOWN', action='WRITE')
  !     'a'   append to end of existing file, or create
  !              (status='UNKNOWN', action='WRITE',ACCESS='APPEND')
  !
  !   The following flags may be added:
  !     'b'    open for binary access (form='UNFORMATTED')
  !     '+'    Also open for read access. (action='READWRITE')
  !     '!'    Allow overwrite of existing files.
  !     '>'    Rename existing files.
  !
  ! Example: 'w+b?' means open for read/write of a binary (unformatted)
  ! file, but if the file already exists, return an error.
  !
  ! RETURN VALUE:
  ! A free UNIT is found automatically, and returned. A negative UNIT
  ! is returned on error.
  ! The error value is the negative of the compiler/OS-dependent IOSTAT,
  ! or an internal error number:
  !   -1001 for an empty filename.
  !   -1002 if no free UNITs are found
  !   -1003 for a syntax error in the mode string,
  !
  ! REMARKS:
  !   File_Open does not keep track of which UNITs are in use.
  !   It just does an INQUIRE to check if a UNIT is available.
  !
  !   We could add a 'scratch file' option, i.e. an 's' as the first
  !   character of mode. But, I suspect this is not necessary.
  !-------------------------------------------------------------------
  include "file_io.inc"
  logical in_use
  character(len=10) file_access ! SEQUENTIAL, APPEND
  character(len= 7) file_status ! NEW, OLD, UNKNOWN
  ! Also REPLACE, SEARCH (F90?)
  character(len=11) file_form   ! FORMATTED, UNFORMATTED
  character(len= 9) file_action ! READ, WRITE, READWRITE
  integer ios

  character(len=FILENAME_MAX_LEN) filename_bak
  integer v
  integer backup_files
  logical noclobber

  ! Strange getenv() errors with RH9/IFC
  !call GetEnv("AMBER_FILEOPEN",filename_bak)
  !if (filename_bak(1:6).eq.'BACKUP') then
  !  read(filename_bak(6:),*) backup_files
  !else
  backup_files=0
  !endif
  if (index(mode,'>') /= 0) then
    backup_files=20
  end if

  if (Len_Trim(filename) == 0) then
    write(stderr,'(A)') &
      'ERROR in OpenFile(): zero-length filename!'
    !call error_handler; stop
    stop 
    !File_Open=-1001
    !return
  end if ! (Len_Trim(filename) == 0)

  !       Derive Fortran OPEN parameters from the mode string.
  file_access = 'SEQUENTIAL'
  file_status = 'UNKNOWN'
  file_action = 'READWRITE'
  if (index(mode,'!') /= 0) then
    noclobber=.false.
  else ! (index(mode,'!') /= 0)
    noclobber=.true.
  end if ! (index(mode,'!') /= 0)
  if (index(mode,'b')+index(mode,'B') > 0) then
    file_form = 'UNFORMATTED'
  else ! (index(mode,'b')+index(mode,'B') > 0)
    file_form = 'FORMATTED'
  end if ! (index(mode,'b')+index(mode,'B') > 0)

  if (mode(1:1) == 'R' .or. mode(1:1) == 'r') then
    file_action='READ'
    file_status='OLD'
  else ! (mode(1:1) == 'R' .or. mode(1:1) == 'r')
    if (index(mode,'+') == 0) then
      file_action='WRITE'
    end if ! (index(mode,'+') == 0)
    if (mode(1:1) == 'W' .or. mode(1:1) == 'w') then
      if (noclobber) then
        file_status='NEW'
      end if ! (noclobber)
      !           Backup files:
      if (noclobber .and. backup_files>0 .and. &
          File_Exists(filename)) then
        filename_bak=filename
        v=0
        do while (File_Exists(filename_bak) .and. v < 99)
          v=v+1
          write(filename_bak,'(2A,I2.2)') filename,'_',v
        end do !  while (File_Exists(filename_bak) .and. v < 99)
        call rename(filename,filename_bak,ios)
        if (ios /= 0) then
          write(stderr,*)'Error creating backup file.'
        end if ! (ios /= 0)
      end if ! (noclobber .and. backup_files .and.
      !
    else if (mode(1:1) == 'A' .or. mode(1:1) == 'a') then ! (mode(1:1) == 'W' .or. mode(1:1) == 'w')
      file_access='APPEND'
      if (noclobber) then
        file_status='OLD'
      end if ! (noclobber)
    else ! (mode(1:1) == 'W' .or. mode(1:1) == 'w')
      !           ERROR: bad mode syntax
      write(stderr,'(A)') &
        'ERROR in OpenFile(): bad file mode syntax: ',mode
      !call error_handler; stop
      stop
      !File_Open=-1003
      !return
    end if ! (mode(1:1) == 'W' .or. mode(1:1) == 'w')
  end if ! (mode(1:1) == 'R' .or. mode(1:1) == 'r')

  ! Find an available unit number:
  File_Open = FILE_UNIT_FIRST-1
  in_use=.true.
  do while (File_Open < FILE_UNIT_LAST .and. in_use)
    File_Open=File_Open+1
    inquire(unit=File_Open,opened=in_use)
  end do

  ! If we have checked all of them, and none available, return error
  if (in_use) then
    write(stderr,'(A)') &
        'ERROR: OpenFile() ran out of available File UNITs'
    File_Open=-1002
    !call error_handler; stop
    stop
  end if ! (in_use)

  !OPEN(UNIT=File_Open,NAME=filename,FORM=file_form, &
  OPEN(UNIT=File_Open,file=filename,FORM=file_form, &
      status=file_status, ACCESS=file_access, IOSTAT=ios)

  if (ios /= 0) then
    File_Open = -ios
    !         Just in case a some systems can return a negative IOSTAT:
    if (File_Open > 0) then
      File_Open = ios - 10000
    end if ! (File_Open > 0)
  end if ! (ios /= 0)
  return
end function file_open 
!---------------------------------------------------------
! Nothing special here, unless we add something for socket I/O
! or tracking of files opened via File_Open().

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine file_close here]
subroutine File_Close(iunit)
  implicit none
  integer iunit
  close(unit=iunit)
end subroutine file_close 

!---------------------------------------------------------------
! Close all open file units in the range allowed for File_Open()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine file_close_all here]
subroutine File_Close_All()
  implicit none
 include "file_io.inc"
  integer iunit
  logical in_use
  character(len=FILENAME_MAX_LEN) filename

  do iunit = FILE_UNIT_FIRST, FILE_UNIT_LAST
    INQUIRE(unit=iunit,opened=in_use,name=filename)
    if (in_use) then
      write(stdout,'(A,I2,2A)') &
          'Closing file: UNIT=',iunit, &
          ' FILENAME=',filename(1:len_trim(filename))
      close(unit=iunit)
    end if ! (in_use)
  end do !  iunit = FILE_UNIT_FIRST, FILE_UNIT_LAST
end subroutine file_close_all 

!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function file_exists here]
function File_Exists(filename)
  implicit none
  logical File_Exists
  character(len=*) filename
  !#include "file_io.inc"
  INQUIRE(FILE=filename,EXIST=File_Exists)
end function
!---------------------------------------------------------
