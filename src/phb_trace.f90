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
program PHB_TRACE
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit
  use debug
  implicit none
  integer version, state_selector, display_progress, delay
  integer term_width, term_height
  double precision :: gridx1, gridx2, gridy1, gridy2
  integer :: gridnx, gridny
  logical :: CCW
  character(len=4096) :: log_file_name
  integer, parameter :: DISPLAY_NONE=0, DISPLAY_PLAIN=1, DISPLAY_ANSI=2
  integer FRAME_START_TIME
  integer ios
  integer total_steps, total_contour
  logical last_done

  namelist /TRACER/ version, log_file_name, state_selector, display_progress, &
       delay, term_width, term_height, gridx1, gridx2, gridy1, gridy2, &
       gridnx, gridny, CCW

  write (stderr,'(A)') "This program reads the problem definition from STD"// &
       "IN, outputs the resulting", "traced contour to STDOUT and prints "// &
       "various progress indicators on STDERR. If", "you are on a terminal,"// &
       " just press Ctrl-D to send EOF and see the program help."

  version = -1
  state_selector = 1
  display_progress = DISPLAY_PLAIN
  delay = 0
  term_width=80
  term_height=25
  log_file_name = ""
  gridx1 = -5.0D0
  gridx2 =  5.0D0
  gridy1 = -5.0D0
  gridy2 =  5.0D0
  gridnx =  50
  gridny =  50
  CCW = .true.

  read (unit=stdin, nml=tracer, iostat=ios)

  total_steps = 0
  total_contour = 0
  last_done = .FALSE.   ! last point was a contour point (for deciding on
                        ! whether to introduce a time delay)

  if ((ios.ne.0).or.(version.EQ.-1)) then
     write (stderr,'(A)') "", "         PHB_TRACE -- program for tracing "// &
          "the phase diagram boundaries of","                 various "// &
          "interesting (mainly micromagnetic) problems.", "", "Copyright "// &
          "(C) 2023 Konstantin L. Metlov <metlov@donfti.ru>.", "", "License"// &
          " GPLv3+: GNU GPL version 3 or later "// &
          "<https://gnu.org/licenses/gpl.html>", "This is free software: "// &
          "you are free to change and redistribute it.", "There is "// &
          "NO WARRANTY, to the extent permitted by law.", "", "", "The "// &
          "program searches for (and traces by moving along) a boundary "// &
          "on a generic", "X-Y plane given by the equation f(X,Y)=0 of "// &
          "a certain (usually very complex and", "hard to compute) "// &
          "function f(X,Y). Specifics of this particular program is that", &
          "it does not assume that a function f(X,Y) can be evaluated at "// &
          "a particular", "point (X,Y) straight away (that would be too "// &
          "easy ;-). Instead, most of the", "considered problems start "// &
          "from a single point inside the region f(X,Y)<0 and", "then "// &
          "carefully move from this point to a neighbouring one in a "// &
          "small (and", "decreasing, if necessary) steps. At first, the "// &
          "search moves along a straight", "line in a particular direction "// &
          "and then, once the boundary f(X,Y)=0 is reached,", "along the "// &
          "said boundary until returns back to the initial boundary point.", &
          "", "The search employs a finite grid to specify the expected "// &
          "feature size of", "the boundary. The considered internal points "// &
          "of the region are restricted to", "this grid, while the "// &
          "boundary points (the ones, the program outputs) are not", &
          "and are always a solution of the f(X,Y)=0 equation to a "// &
          "certain (very small)", "numerical precision. Finer the "// &
          "resolution of the grid, more boundary points", "will be "// &
          "computed and printed to STDOUT. Each output line contains "// &
          "the X and", "Y coordinates of a boundary point as a space-"// &
          "separated ASCII pair of numbers.", "", "Grid size is finite. "// &
          "Once the search reaches the grid boundary, it continues", &
          "along as if the grid boundary was the phase boundary.", "", &
          "Most problems enable the checkpointing feature, storing the "// &
          "intermediate", "computation result to a binary platform-"// &
          "specific LOG file. If the tracing", "was interrupted, a simple "// &
          "restart would continue from the last computed point.", &
          "If the tracing was completed and the complete log file is "// &
          "available, no", "computation is done and the program simply "// &
          "prints the previously found boundary.", "", "The STDIN should "// &
          "contain the problem definition as a series of Fortran", &
          "namelists. The first namelist defines the grid, the problem to "// &
          "be solved and", "the level of progress reporting. The second "// &
          "(problem-specific) namelist", "defines the initial point and "// &
          "the meaning of the coordinates X and Y. See", "the example "// &
          "parameter files for specifics."
     stop 9
  end if

  if (version.NE.1) then
     write (stderr,'(/,"ERROR: version = ", I0, " is invalid.")') version
     stop 9
  end if

  if ((display_progress.LT.0).or.(display_progress.GT.2)) then
     write (stderr,'(/,"ERROR: display_progress = ", I0, " is invalid.")') &
          display_progress
     stop 9
  end if

  if (delay.LT.0) then
     write (stderr,'(/,"ERROR: delay = ", I0, " is invalid.")') delay
     stop 9
  end if

  if (display_progress.EQ.DISPLAY_ANSI) then
     write(stderr,'(a)') char(27) // "[2J"  ! clear screen
  end if

  if (delay.GT.0) call system_clock(FRAME_START_TIME)

  select case (state_selector)
  case (1)
    call process_test_state
  case (2)
    call process_hopfion
  case (3)
    call process_skyrmion
  case default
     write (stderr,'(/,"ERROR: state_selector = ", I0, " is invalid.")') &
          state_selector
     stop 9
  end select

  if (display_progress.EQ.DISPLAY_PLAIN) write(stderr,*)

contains

  subroutine process_test_state
    use test_state
    type (testing_state) :: s
    double precision sx,sy
    namelist /INITIAL_TESTING_STATE/ sx,sy

    sx=0
    sy=-3
    read (unit=stdin, nml=INITIAL_TESTING_STATE, iostat=ios)
    if (ios.NE.0) then
       write (stderr,'(/,"ERROR: can''t load INITIAL_TESTING_STATE.")')
       stop 9
    end if
    call s%init(sx,sy)
    call do_trace(s)
  end subroutine process_test_state

  subroutine process_hopfion
    use hopf_state
    type (hopfion_state) :: s
    integer ios

    call s%load_nml(stdin, ios)
    if (ios.NE.0) then
       write (stderr,'(/,"ERROR: can''t load HOPFION_STATE.")')
       stop 9
    end if

    call do_trace(s)
  end subroutine process_hopfion

  subroutine process_skyrmion
    use skyrme_state
    type (skyrmion_state) :: s
    integer ios

    call s%load_nml(stdin, ios)
    if (ios.NE.0) then
       write (stderr,'(/,"ERROR: can''t load SKYRMION_STATE.")')
       stop 9
    end if

    call do_trace(s)
  end subroutine process_skyrmion

  subroutine REPORT_CALLBACK(X, Y, done)
    implicit none
    double precision, intent(IN):: X, Y
    logical, intent(IN):: done
    double precision aspect, IW, IH, scale
    integer tx, ty
    integer t1, dt
    character(len=8) :: date_string
    character(len=8) :: time_string

    if (done) then
       write(stdout,'(G0.10,1X,G0.10)') X, Y
       flush(stdout)
    end if
    if (done) then
       total_contour = total_contour + 1
    else
       total_steps  = total_steps + 1
    end if
    select case (display_progress)
    case (DISPLAY_NONE)
    case (DISPLAY_PLAIN)
       if (done) then
          write(stderr,'("BOUNDARY ")',ADVANCE="NO")
       else
          write(stderr,'("TO ")',ADVANCE="NO")
       end if
       write(stderr,'("{",G0.2,",",G0.2,"}; ")',ADVANCE="NO") X,Y
       if (done) write(stderr,*)
    case (DISPLAY_ANSI)
       aspect = (1.0D0*gridny)/(1.0D0*gridnx)
       IW = gridx2-gridx1
       IH = gridy2-gridy1
       scale = min(real(term_width/2)/real(IW), real(term_height-1)/real(IH))
       tx=NINT((X-gridx1)*2*scale)
       ty=NINT((gridy2-Y)*scale)
       write(stderr,'(a,"[",I0,a,I0,a)',ADVANCE="NO") &
            char(27), ty, ";", tx, "H"
       if (done) then
          write(stderr,'(a)',ADVANCE="NO") "++"
       else
          write(stderr,'(a)',ADVANCE="NO") ".."
       end if
       write(stderr,'(a,"[",I0,a,I0,a)',ADVANCE="NO") &
            char(27), term_height, ";", 0, "H"
       call DATE_AND_TIME(date_string, time_string)
       write(stderr,'(" updated at ", a,":",a,":",a, " contour/steps=", &
            &I0, "/", I0)', ADVANCE="NO") &
            time_string(1:2), time_string(3:4), time_string(5:8), &
            total_contour, total_steps
       write(stderr,'(a,"[",I0,a,I0,a)',ADVANCE="NO") &
            char(27), term_height, ";", 0, "H"
    case default
       CALL ASSERT(.FALSE.)
    end select

    ! wait if frame is computed faster than delay
    if (.not.last_done .and. (delay.GT.0)) then
       do
          call system_clock(t1)
          dt = t1 - FRAME_START_TIME
          if (dt >= delay) exit
       end do
       FRAME_START_TIME = t1
    end if
    last_done = done
  end subroutine REPORT_CALLBACK

  subroutine do_trace(s)
    use tracer
    use iso_c_binding
    implicit none
    class (state), intent(INOUT) :: s
    integer lu
    integer, parameter :: NMAX=5000
    integer N
    double precision X(NMAX), Y(NMAX)
    integer IERROR

    if (LEN_TRIM(log_file_name).GT.0) then
       open(newunit=lu, file=log_file_name, access="stream", &
            form="unformatted", status="UNKNOWN", iostat=ios)
       if (ios.NE.0) then
          write (stderr,'(/,"ERROR: can''t open the specified LOG file.")')
          stop 9
       end if
       call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, CCW, &
            NMAX, N, X, Y, IERROR, lu, REPORT_CALLBACK)
       close(lu)
    else
       call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, CCW, &
            NMAX, N, X, Y, IERROR, NULL() , REPORT_CALLBACK)
    end if
    if (IERROR.EQ.15) then
       write (stderr, '("Log file """,A,""" exists and is incompatible with&
            & the current grid.")') TRIM(log_file_name)
       stop 9
    end if
    if (IERROR.NE.0) then
       write (stderr,'(/,"ERROR: tracer returned IERROR= ", I0, ".")') IERROR
       stop 9
    end if
  end subroutine do_trace

end program PHB_TRACE
