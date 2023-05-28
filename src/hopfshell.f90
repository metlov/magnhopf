!
!     MAGNHOPF -- MAGNETIC HOPFIONS LIBRARY AND A SET OF TOOLS
!
!     (c) 2023 Konstantin L. Metlov <metlov@donfti.ru>
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!     THIS FILE CONTAINS PHASE BOUNDARY TRACING UTILITY PHB_TRACE
!
program HOPFSHELL
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit
  use jsu_readline
  use hopf_state
  implicit none
  class(skyrmion_state), allocatable :: currstate
  integer :: currstate_type ! 1 for hopfion (of either type), 2 for skyrmion
  integer, parameter :: lineMaxLen=1024
  character(len=lineMaxLen):: line
  character(len=*), parameter :: last_state_filename = ".last_hopfion_state"
  logical display
  character (lineMaxLen) command, rest

  call initialize_state
  print *,"Enter command ('q' to quit, 'h' to display help):"
  display=.true.
  do
     if(display) call display_state
     display=.true.
     call iso_readline(line,"HOPFSHELL> ")
     call parse_first_word(line, command, rest)
     call tolowercase(command)
     if(command.eq.'q'.or.command.eq.'quit') then
        exit
     elseif(command.eq.'h'.or.command.eq.'help') then
        display=cmd_help(rest)
     elseif(command.eq.'e'.or.command.eq.'evolve') then
        display = cmd_evolve(rest)
     elseif(command.eq.'n'.or.command.eq.'namelist_export') then
        display = cmd_namelist_export(rest)
     elseif(command.eq.'p'.or.command.eq.'predicate_set') then
        display = cmd_predicate_set(rest)
     elseif(command.eq.'t'.or.command.eq.'type_set') then
        display = cmd_type_set(rest)
     end if
     call save_state(last_state_filename)
  enddo
  ! delete last state on successful exit
  call delete_last_state
  deallocate(currstate)
contains
  subroutine tolowercase(str)
    character(len=*), intent(inout) :: str   ! input/output string
    integer :: i
    do i = 1, len(str)
       if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
          str(i:i) = char(ichar(str(i:i)) + 32) ! convert to lowercase
       endif
    enddo
  end subroutine tolowercase
  ! parses the first word in line
  subroutine parse_first_word(s, command, rest)
    implicit none
    character(len=*), intent(IN) :: s
    character(len=*), intent(OUT) :: command, rest
    integer start_idx, end_idx

    start_idx=1
    do while (start_idx.LE.LEN_TRIM(s).AND.s(start_idx:start_idx).EQ.' ')
       start_idx = start_idx + 1
       if (start_idx.GT.LEN_TRIM(s)) exit
    end do
    end_idx=start_idx
    do while ((end_idx.LE.LEN_TRIM(s)).AND.(s(end_idx:end_idx).NE.' '))
       end_idx = end_idx + 1
       if (end_idx.GT.LEN_TRIM(s)) exit
    end do
    if (start_idx.eq.LEN_TRIM(s)+1) then
       ! line is an empty string
       command=""
       rest=""
    else
       command=s(start_idx:end_idx-1)
       if (end_idx.GT.LEN_TRIM(s)) then
          rest = ""
       else
          rest=s(end_idx:LEN_TRIM(s))
       end if
    endif
  end subroutine parse_first_word

  subroutine initialize_state
    implicit none
    logical last_state_exists, load_last_state
    inquire(file=last_state_filename, exist=last_state_exists)
    load_last_state=last_state_exists
    if (last_state_exists) then
       print '(A)', "The file .last_hopfion_state exists.", "", &
            "It probably means that the last computation in the HOPFSHELL "// &
            "was aborted.", "Would you like to load this last state into "//&
            "the shell or start from scratch ?", "The state will be deleted"// &
            " in both cases.", ""
       do
          call iso_readline(line,"Please answer Y(es) or N(o) : ")
          line=adjustl(line)
          if (line.eq."y".or.line.eq."Y".or. &
               line.eq."n".or.line.eq."N") exit
       end do
       if(line.eq.'n'.or.line.eq.'N') load_last_state=.false.
    end if
    if (load_last_state) then
       call load_state(last_state_filename)
    else
       currstate_type = 1
       allocate(hopfion_state :: currstate)
       call currstate%init(0.0D0,0.0D0)
    end if
    if (last_state_exists) call delete_last_state
  end subroutine initialize_state

  subroutine load_state(state_filename)
    implicit none
    character(len=*) :: state_filename
    integer unit, ios
    open(newunit=unit, file=last_state_filename, access="stream", &
         form="unformatted", status="old", action="read", iostat=ios)
    CALL ASSERT(ios.eq.0)  ! we always must end up in the correct state
    read (unit, iostat=ios) currstate_type
    CALL ASSERT(ios.eq.0)
    if (currstate_type.EQ.1) then
       allocate(hopfion_state :: currstate)
    elseif (currstate_type.EQ.2) then
       allocate(skyrmion_state :: currstate)
    else
       CALL ASSERT(.false.)
    end if
    call currstate%load(unit,ios)
    CALL ASSERT(ios.eq.0)  ! we always must end up in the correct state
    close(unit)
    if (ios.ne.0) print *, "File close error."
  end subroutine load_state

  subroutine save_state(state_filename)
    implicit none
    character(len=*) :: state_filename
    integer unit, ios
    open(newunit=unit, file=last_state_filename, access="stream", &
         form="unformatted", status="replace", action="write", iostat=ios)
    if (ios.ne.0) then
       print '("File open error ios=",i0,".")', ios
    else
       write (unit, iostat=ios) currstate_type
       if (ios.ne.0) print *, "File write error."
       call currstate%write(unit,ios)
       if (ios.ne.0) print *, "File write error."
       close(unit, iostat=ios)
       if (ios.ne.0) print *, "File close error."
    end if
  end subroutine save_state

  subroutine delete_last_state()
    implicit none
    integer unit, ios
    open(newunit=unit, iostat=ios, file=last_state_filename, status='old')
    if (ios == 0) close(unit, status='delete')
  end subroutine delete_last_state

  subroutine display_state
    use groundstate
    logical, parameter :: fancy=.true.
    integer, parameter :: DW=79, DH=20
    double precision :: D(DW,DH), Dmax
    integer i, j, ri, fi
    double precision q, h
    logical out_of_range
    integer, parameter :: maxPredicateHelpLen=70
    character(len=maxPredicateHelpLen) :: buf, bufShort
    character(len=2) :: code
    integer stateCode
    double precision currE, stateE
    character(len=3) :: stateSymb, exclamation

    D = 0.0D0

    out_of_range=.false.
    do i=1,currstate%NFPTS
       ri=1+nint(currstate%FPTS(i,1)*(DW-1))
       fi=1+nint(currstate%FPTS(i,2)*(DH-1))
       if (ri.ge.1.and.ri.le.DW.and.fi.ge.1.and.fi.le.DH) then
          D(ri,fi) = D(ri,fi) + 1
       else
          out_of_range = .true.
       end if
    end do
    if (out_of_range) print *, "WARNING: some values vere out of range."

    Dmax=-1
    do i=DH,1,-1
       do j=1,DW
          if (D(j,i).gt.Dmax) Dmax=D(j,i)
       end do
    end do
    do i=DH,1,-1
       do j=1,DW
          if (D(j,i).gt.0.0D0) then
             if (fancy) then
                select case (1+nint(3.0D0*(D(j,i)/Dmax)))
                case (1)
                   write (*,'(a)', ADVANCE="NO") '░'
                case (2)
                   write (*,'(a)', ADVANCE="NO") '▒'
                case (3)
                   write (*,'(a)', ADVANCE="NO") '▓'
                case (4)
                   write (*,'(a)', ADVANCE="NO") '█'
                case default
                   write (*,'(a)', ADVANCE="NO") '?'
                end select
             else
                write (*,'(a)', ADVANCE="NO") '.'
             end if
          else
             write (*,'(a)', ADVANCE="NO") ' '
          end if
       end do
       write (*,*)
    end do

    ! status line
    call currstate%getxy(q,h)
    call currstate%predicateInfo(-1, maxPredicateHelpLen, buf, bufShort)
    call currstate%get_short_code(code)
    call getGroundState(q, h, currstate%mu, stateCode, stateE)
    call getStateSymbol(stateCode, stateSymb)
    currE=currstate%energy()
    exclamation=""
    if (currE.LT.stateE) then ! lower than the previously known ground state
       if (.NOT.fancy) then
          exclamation="!"
       else
          exclamation="‼"
       end if
    end if
    print '(A,":",A,": {q,h,mu}={",G0.4,",",G0.4,",",G0.4,"}",&
         &", eTot=",G0.5, " [",A,"=", G0.5,A,"]")', &
         TRIM(bufShort),code, q, h, &
         currstate%mu, currE, &
         stateSymb, stateE, TRIM(exclamation)
  end subroutine display_state

!-------------------------------------------------------
! commands implementation
  recursive function cmd_help(params) result(display)
    implicit none
    character(len=*), intent(IN) :: params
    character(len=lineMaxLen) command,rest
    logical display
    integer i
    integer, parameter :: maxPredicateHelpLen=70
    character(len=maxPredicateHelpLen) :: buf, bufShort
    call parse_first_word(params, command, rest)
    call tolowercase(command)
    if (LEN_TRIM(command).GT.0.AND.LEN_TRIM(rest).GT.0) then
       display = cmd_help("help")
    else
       if (LEN_TRIM(command).EQ.0) then
          print '(A)'," Here is the list of available commands:", &
               "   quit  -- quits the program.", &
               "   help  -- prints this help or individual command usage.", &
               "   evolve  -- tries to evolve the hopfion towards a point"// &
               " in 3D parameter space.", &
               "   namelist_export  -- exports the current state in "// &
               "Fortran namelist format", &
               "   type -- set hopfion type and reinit to h=0, q=0, mu=0 ."
       else
          if (command.eq."q".or.command.eq."quit") then
             print '(A)',"Usage: quit"
          end if
          if (command.eq."h".or.command.eq."help") then
             print '(A)',"Usage: help [command]"
          end if
          if (command.eq."e".or.command.eq."evolve") then
             print '(A)',"Usage: evolve q h mu"
          end if
          if (command.eq."n".or.command.eq."namelist_export") then
             print '(A)',"Usage: namelist_export [<filename.nml>]"
          end if
          if (command.eq."p".or.command.eq."predicate_set") then
             print '(A)',"Usage: predicate_set <predicate_idx>", &
                  "where predicate idx is one of:"
             do i=0,currstate%predicateMax()
                call currstate%predicateInfo(i,maxPredicateHelpLen,buf,bufShort)
                print '(i2," - ",A)', i, buf
             end do
          end if
          if (command.eq."t".or.command.eq."type_set") then
             print '(A)',"Usage: type_set < 1 | 2 | 3 >", &
                  "where type stands for:", &
                  "  1  -- hopfion of the Type 1,", &
                  "  2  -- hopfion of the Type 2,", &
                  "  3  -- skyrmion."
          end if
       end if
       display = .false.
    endif
  end function cmd_help

  logical function cmd_evolve(params)
    implicit none
    character(len=*), intent(IN) :: params
    double precision q, h, mu
    double precision qt, ht, mut, t
    integer ios
    read (params,*,iostat=ios) qt, ht, mut
    if (ios.ne.0) then
       cmd_evolve = cmd_help("evolve")
    else
       call currstate%getxy(q,h)
       mu=currstate%mu
       t = currstate%evolve3D(qt,ht,mut,.true.)
       if (t.eq.1.0D0) then
          print '("Successfully reached the target state.")'
          cmd_evolve = .true.
       else
          print '("Went only ", G0.3, " of the way to {q,h,mu}={", '// &
               'G0.5, ",", G0.5, ",", G0.5, "}. Old state restored.")', &
               t, q+t*(qt-q), h+t*(ht-h), mu+t*(mut-mu)
          cmd_evolve = .false.
       end if
    end if
  end function cmd_evolve
  logical function cmd_namelist_export(params)
    implicit none
    character(len=*), intent(IN) :: params
    character(len=4096) filename
    logical file_exists, do_export
    integer unit, ios

    cmd_namelist_export = .false.

    unit = stdout
    do_export= .true.
    filename = TRIM(ADJUSTL(params))
    if (LEN(filename).GT.0) then
       inquire(file=filename, exist=file_exists)
       if (file_exists) then
          print '("The file """,A,""" exists.")', TRIM(filename)
          print *, "You can answer Y(es) to overwrite, N(o) to abort, "//&
               "or A(ppend) to append the data to the existing file."
          do
             call iso_readline(line,"Please answer Y(es), N(o) or A(ppend) : ")
             line=adjustl(line)
             if (line.eq."y".or.line.eq."Y".or. &
                  line.eq."n".or.line.eq."N".or. &
                  line.eq."a".or.line.eq."A") exit
          end do
          if (line.eq."y".or.line.eq."Y") then
             open(newunit=unit, file=filename, status="replace", &
                  action="write", iostat=ios)
          elseif (line.eq.'a'.or.line.eq.'A') then
             open(newunit=unit, file=filename, status="old", &
                  position="append", action="write", iostat=ios)
          else
             do_export=.false.
             ios = 0
          endif
          if (ios.ne.0) then
             print '("File open error ios=",i0,".")', ios
          end if
       else
          open(newunit=unit, file=filename, status="new", &
               action="write", iostat=ios)
       end if
    end if
    if (do_export) then
       call currstate%write_nml(unit, ios)
       if (ios.ne.0) then
          print '("File write error ios=",i0,".")', ios
       end if
    end if
    if (unit.NE.stdout) close(unit)
  end function cmd_namelist_export

  logical function cmd_predicate_set(params)
    implicit none
    character(len=*), intent(IN) :: params
    integer pred
    integer ios
    read (params,*,iostat=ios) pred
    if (ios.ne.0) then
       cmd_predicate_set = cmd_help("predicate_set")
    else
       if (pred.LT.0.OR.pred.GT.currstate%predicateMax()) then
          print '("ERROR: the predicate index must be in range [0,",I0,"]")', &
               currstate%predicateMax()
          cmd_predicate_set = .false.
       else
          if (.not.currstate%predicateSet(pred)) then
             print *, "ERROR: the current state must satisfy the predicate."
             cmd_predicate_set = .false.
          else
             cmd_predicate_set = .true.
          end if
       end if
    end if
  end function cmd_predicate_set
  logical function cmd_type_set(params)
    implicit none
    character(len=*), intent(IN) :: params
    integer type
    integer ios
    read (params,*,iostat=ios) type
    if ((ios.ne.0).or.(type.lt.1).or.(type.gt.3)) then
       cmd_type_set = cmd_help("type_set")
    else
       select case(type)
       case (1:2)
          if (currstate_type.ne.1) then
             deallocate(currstate)
             allocate(hopfion_state :: currstate)
             currstate_type = 1
          end if
          call currstate%initHopfion(0.0D0, 0.0D0, type, .false.)
          cmd_type_set = .true.
       case (3)
          if (currstate_type.ne.2) then
             deallocate(currstate)
             allocate(skyrmion_state :: currstate)
             currstate_type = 2
          end if
          call currstate%init(0.0D0, 0.0D0)
       case default
          call ASSERT(.false.)
       end select
    end if
  end function cmd_type_set
end program HOPFSHELL
