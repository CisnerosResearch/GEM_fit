
include "io.fh"
      !include 'memory.fh'
      !include 'database.fh'
      !include 'strings.fh'

      ! File_Open will allocate the next unused UNIT number in
      ! the range of FILE_UNIT_FIRST to FILE_UNIT_LAST.
      ! I start at 20 to leave low-numbered units for direct OPENs.
      ! I stop at 50, because that's originally the highest allowed,
      ! and 30 file units should be plenty.
      integer FILE_UNIT_FIRST; parameter(FILE_UNIT_FIRST=20)
      integer FILE_UNIT_LAST; parameter(FILE_UNIT_LAST=50)
