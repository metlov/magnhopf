module debug
  implicit none

  public :: CHKASSERT, ASSERT
  public :: ASSERTS_ENABLED

  private
  logical, parameter :: ASSERTS_ENABLED = .true.

contains
  subroutine CHKASSERT(COND, FNAME, LINE, MSG)
    logical COND
    character(*) MSG,FNAME
    integer LINE
    if (.not.COND) then
       write(0,'(A,A,A,I0,A,A)') &
            "ASSERTION FAILED in ",FNAME," at ", LINE, ": ", MSG
       error stop
    end if
  end subroutine CHKASSERT

  subroutine ASSERT(COND)
    logical, intent(IN) :: COND
    if (.not.COND) then
       error stop "ASSERTION FAILED"
    end if
  end subroutine ASSERT

end module debug
