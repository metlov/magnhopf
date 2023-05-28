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
!     SUBROUTINES IN THIS FILE ORGANIZE INCREMENTAL TRACING OF THE PHASE
!     DIAGRAM BOUNDARIES.
module tracer
  use debug
  implicit none

  ! an abstract type for the computation state, representing a point on
  ! phase diagram
  type, abstract :: state
     private
     double precision :: X, Y
   contains
     ! constructor
     procedure :: init => state_create

     ! returns the coordinates of the present state on the phase diagram
     procedure :: getxy => state_getxy

     ! evolves the state towards a specified point on the phase
     ! diagram. Returns the fraction of the distance, which can be
     ! successfully traversed. A special value of 1.0D0 on return means
     ! that the requested point on the phase diagram was successfully
     ! reached. Only in that case (when the destination was actually reached),
     ! this function must advance the state of this object, otherwise the
     ! state must be left unchanged. This condition is necessary to ensure that
     ! the state always stays at grid points, which are fully within the
     ! traced region. If this contract is broken, an assertion will fail
     ! in tracer.
     procedure :: evolve => state_evolve

     ! appends the current state of the computation to a binary
     ! file at the current position
     procedure :: write => state_write

     ! loads the previously written state from the current position in a
     ! binary unit into the current object
     procedure :: load => state_load

     ! prints the state on stdout for debugging purposes
     procedure :: print => state_print
  end type state

  public :: TRACE

contains
  ! constructor to initialize the object
  subroutine state_create(this, xs, ys)
    class(state), intent(OUT) :: this
    double precision, intent(IN) :: xs, ys
    this%X = xs
    this%Y = ys
  end subroutine state_create

  ! a getter method for state coordinates
  subroutine state_getxy(this, xs, ys)
    implicit none
    class(state), intent(IN) :: this
    double precision, intent(OUT) :: xs, ys
    xs = this%X
    ys = this%Y
  end subroutine state_getxy

  ! a setter method for state coordinates
  double precision function state_evolve(this, xt, yt)
    implicit none
    class(state), intent(INOUT) :: this
    double precision, intent(IN) :: xt, yt
    this%X = xt
    this%Y = yt
    state_evolve = 1.0D0    ! always a success here
  end function state_evolve

  ! state writer
  subroutine state_write(this, unit, ios)
    implicit none
    class(state), intent(IN) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    character(len=12) :: u_form, u_act, u_access
    logical ok_form, ok_act, ok_access

    if (ASSERTS_ENABLED) then
       inquire(unit, access=u_access, action=u_act, form=u_form, iostat=ios)
       if (ios /= 0) return
       ok_form =trim(u_form) == "UNFORMATTED"
       CALL ASSERT(ok_form)
       ok_act = trim(u_act)=="READWRITE" .OR. trim(u_act)=="WRITE"
       CALL ASSERT(ok_act)
       ok_access = trim(u_access).eq."STREAM"
       CALL ASSERT(ok_access)
    end if

    write(unit,iostat=ios) this%X, this%Y
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine state_write

  ! state loader
  subroutine state_load(this, unit, ios)
    class(state), intent(OUT) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    character(len=12) :: u_form, u_act, u_access
    logical ok_form, ok_act, ok_access

    if (ASSERTS_ENABLED) then
       inquire(unit, access=u_access, action=u_act, form=u_form, iostat=ios)
       if (ios /= 0) return
       ok_form =trim(u_form) == "UNFORMATTED"
       CALL ASSERT(ok_form)
       ok_act = trim(u_act)=="READWRITE" .OR. trim(u_act)=="READ"
       CALL ASSERT(ok_act)
       ok_access = trim(u_access).eq."STREAM"
       CALL ASSERT(ok_access)
    end if

    read(unit,iostat=ios) this%X, this%Y
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine state_load

  ! a getter method for state coordinates
  subroutine state_print(this)
    implicit none
    class(state), intent(IN) :: this
    print '("X=",G0.5,", Y=",G0.5)', this%X, this%Y
  end subroutine state_print


  ! The main phase diagram boundary tracing routine
  !
  ! The phase diagram is traced in a generic X-Y coordinates. Their specific
  ! meaning is up to the definition of a particular implementation of the
  ! state type.
  !
  ! The tracing uses grid, described by gridx1, gridy1, gridx2, gridy2,
  ! gridnx, gridny parameters and only traces the contour, lying inside
  ! the grid. If the traced contour reaches the grid boundary, then the
  ! grid boundary itself becomes part of the contour. The contour itself
  ! is not firmly restricted to the grid though and its points do not
  ! necessarily coincide with grid points. They do always lie on grid
  ! lines though. So, the grid size impacts the feature size, which can be
  ! resolved by the tracer (the number of points along the grid to be
  ! computed), but not the precision of the individual contour points.
  !
  ! The tracing procedure maintains an optional log file, which allows it
  ! to restart an interrupted computation and cache previously computed
  ! states. To use this feature, no action is required on caller's part
  ! besides specifying a (valid) log file name. If the valid log file
  ! exists, computation will automatically restart. If the state
  ! evolver is simple and fast, omitting the log file argument will make
  ! the tracing run entirely in memory without caching (in this case
  ! the load and write members of the state type will never be called
  ! and their implementation is not necessary).
  !
  ! Paramerers:
  !    Istate  -- initial state of the tracer. The closest grid point
  !               to the coordinates of this state must be inside
  !               the traced region of the phase diagram.
  !    gridx1  -- X of the lower left point of the phase diagram grid
  !    gridy1  -- Y of the lower left point of the phase diagram grid
  !    gridx2  -- X of the upper right point of the phase diagram grid
  !    gridy2  -- Y of the upper right point of the phase diagram grid
  !    gridnx  -- number of grid points along the X axis
  !    gridnx  -- number of grid points along the Y axis
  !    CCW     -- if .true. the contour is walked counterclockwise, if
  !               .false. -- clockwise
  !    NMAX    -- the maximum number of contour points to store in the
  !               resulting array
  !    N       -- the actual number of computed contour points
  !    X       -- the X coordinates of the computed contour points
  !    Y       -- the Y coordinates of the computed contour points
  !    IERROR  -- 0 if the computation completed successfully, otherwise
  !               signals the type of the encountered error:
  !               10 -- the grid point, closest to the initial state, lies
  !                     outside of the specified grid.
  !               11 -- the grid point, closest to the initial state, lies
  !                     outside of the traced region of the phase diagram;
  !                     evolving the state towards it crossed the phase boundary
  !               12 -- supplied log unit is not usable for log.
  !               13 -- the state advancement (evolve) function returned
  !                     a value, which is out of range [0.0, 1.0]
  !               14 -- not enough space in X,Y arrays for all the computed
  !                     contour points
  !               15 -- LOG file header mismatch, in case the grid was changed
  !                     please delete the old LOG file to start from scratch
  !    logunit -- (optional) opened file unit to use for log storage;
  !               must be opened with access="stream" and form="unformatted"
  !    NEWCPCB -- (optional) a call back function to be called when a
  !               new contour point is computed
  subroutine TRACE(Istate, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, CCW,&
       NMAX, N, X, Y, IERROR, logunit, NEWCPCB)
    implicit none
    class(state), intent(INOUT) :: Istate
    double precision, intent(IN) :: gridx1, gridy1, gridx2, gridy2
    integer, intent(IN) :: gridnx, gridny
    logical, intent(IN) :: CCW
    integer, intent(IN) :: NMAX
    integer, intent(OUT) :: N
    double precision, intent(OUT) :: X(NMAX), Y(NMAX)
    integer, intent(OUT) :: IERROR
    interface
       ! Optional callback, executed when a new grid line is traversed.
       !
       ! Parameters:
       !     X, Y    -- coordinates of the currently computed point on the
       !                phase diagram
       !     done    -- .FALSE. means that the above specified grid point is
       !                only being evaluated with respect to its belonging
       !                to the contour; .TRUE. means the evaluation is complete
       !                and Xi,Yi specify a newly found contour point.
       subroutine NEWCPCALLBACK_func(X, Y, done)
         double precision, intent(IN):: X, Y
         logical, intent(IN):: done
       end subroutine NEWCPCALLBACK_func
    end interface
    integer, intent(IN), optional :: logunit
    procedure (NEWCPCALLBACK_func), optional :: NEWCPCB

    logical have_callback
    double precision xcbc,ycbc,xcbn,ycbn
    double precision xc, yc, v
    integer ic, jc, ld, nd, neighbors

    double precision gsx, gsy

    integer, parameter :: diri(4) = [ 0, 1, 0,-1]
    integer, parameter :: dirj(4) = [-1, 0, 1, 0]

    !   1 2 3 4  d             direction
    !   3 4 1 2  revd    reverse direction
    integer, parameter :: revd(4) = [3, 4, 1, 2]

!    integer, parameter :: dirnext(16,5) = reshape( &
!         [ 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 0, &
!           4, 1, 4, 2, 4, 1, 4, 1, 4, 2, 4, 2, 4, 3, 4, 0, &
!           1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 0, &
!           2, 2, 2, 2, 3, 3, 4, 1, 2, 2, 2, 2, 3, 3, 4, 0, &
!           3, 3, 4, 1, 3, 1, 1, 1, 3, 3, 4, 2, 3, 3, 4, 0 ], &
!         [16, 5])

    integer, parameter :: dirnextCCW(16,4) = reshape( &
         [ 4, 1, 4, 2, 4, 1, 4, 1, 4, 2, 4, 2, 4, 3, 4, 0, &
           1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 0, &
           2, 2, 2, 2, 3, 3, 4, 1, 2, 2, 2, 2, 3, 3, 4, 0, &
           3, 3, 4, 1, 3, 1, 1, 1, 3, 3, 4, 2, 3, 3, 4, 0 ], &
           [16, 4])
    integer, parameter :: ld0CCW = 2

    integer, parameter :: dirnextCW(16,4) = reshape( &
         [ 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 4, 3, 4, 0, &
           3, 3, 2, 2, 3, 3, 1, 1, 3, 3, 2, 2, 3, 3, 4, 0, &
           4, 3, 4, 2, 4, 3, 4, 1, 4, 3, 4, 2, 4, 3, 4, 0, &
           1, 1, 1, 1, 1, 1, 1, 1, 4, 3, 4, 2, 4, 3, 4, 0 ], &
           [16, 4])
    integer, parameter :: ld0CW = 3

    integer :: dirnext(16,4)
    integer :: ld0

    ! search log in memory
    integer logMaxSize ! approx max number of contour pts
    integer, parameter :: logMagic = -87114001 ! Z'FACEBEEF'
    integer, dimension(:), allocatable :: log_oi, log_oj, log_d, log_seek
    double precision, dimension(:), allocatable :: log_v
    integer logEntrySize ! determined via inquire
    integer logSize

    ! log file on disk
    character(len=12) :: u_form, u_act, u_access
    logical ok_form, ok_act, ok_access
    integer ios

    ! map of points
    ! the map array is going to be very sparse. If memory is a concern, it
    ! is best to replace it with a hash table. But for now let's use an array,
    ! which is easier to handle in Fortran.
    double precision,allocatable :: map (:,:,:)

    integer i, j, k, id, jd
    logical out_of_the_grid

    integer :: startpos, startdatapos, magic
    logical fullyread

    logical, parameter :: DEBUG_TRAVERSAL=.false.
    logical, parameter :: DEBUG_LOGFILE=.false.

    if (CCW) then
       dirnext = dirnextCCW
       ld0 = ld0CCW
    else
       dirnext = dirnextCW
       ld0 = ld0CW
    end if

    have_callback=PRESENT(NEWCPCB)

    ! check the supplied log file
    if (PRESENT(logunit)) then
       inquire(logunit, access=u_access, action=u_act, form=u_form, iostat=ios)
       if (ios.NE.0) then
          IERROR = 12
          return
       end if
       ok_form =trim(u_form) == "UNFORMATTED"
       ok_act = trim(u_act)=="READWRITE"
       ok_access = trim(u_access).eq."STREAM"
       if ((.NOT.ok_form).OR.(.NOT.ok_act).OR.(.NOT.ok_access)) then
          IERROR = 12
          return
       end if
    end if

    ios = 0   ! assume no I/O error in case log is disabled
    N = 0     ! no contour points computed yet

    gsx=(gridx2-gridx1)/gridnx
    gsy=(gridy2-gridy1)/gridny

    ! allocate and initialize the map
    allocate(map(1:4,0:gridnx,0:gridny))
    do i=0,gridnx
       do j=0,gridny
          do k=1,4
             map(k,i,j)=-1.0D0
          end do
       end do
    end do

    ! allocate and initialize log
    logSize=0
    logMaxSize = 3*NMAX+gridny  ! roughly the expected log size
    allocate(log_oi(logMaxSize))
    allocate(log_oj(logMaxSize))
    allocate(log_d(logMaxSize))
    allocate(log_v(logMaxSize))
    allocate(log_seek(logMaxSize))

    ! if supplied, read disk log into memory
    if (PRESENT(logunit)) then
       rewind(logunit, iostat=ios)
       if (ios.NE.0) then
          IERROR = 12
          return
       end if
       if (DEBUG_LOGFILE) &
            write (*, '("Loading previously computed LOG: ")')

       inquire(iolength=logEntrySize) log_oi(1), log_oj(1), &
            log_d(1), log_v(1), log_seek(1)
       do
          if (logSize.EQ.0) then ! check header
             if (.NOT.log_header_matches_or_none(ios)) then
                ! old log with different grid size is present
                ! bail out to give user a chance to save data
                IERROR = 15
                return
             end if
             if (ios.ne.0) then
                rewind(logunit, iostat=ios)
                CALL ASSERT(ios.EQ.0)
                exit  ! out of the loop, log will be overwritten
             end if
          end if

          if (logSize.GE.logMaxSize) then
             IERROR = 14
             return
          end if

          fullyread=.false.
          inquire(unit=logunit, pos=startpos)
          read (logunit, iostat=ios) log_oi(logSize+1), log_oj(logSize+1), &
               log_d(logSize+1), log_v(logSize+1), log_seek(logSize+1)
          if (DEBUG_LOGFILE.and.(ios.eq.0)) then
             print '("{",I0,", ",I0,", ",I0,", ",G0.5,", ",I0,"}")', &
                  log_oi(logSize+1), log_oj(logSize+1), &
                  log_d(logSize+1), log_v(logSize+1), log_seek(logSize+1)
          end if
          if ((ios.eq.0).and.&
               (log_oi(logSize+1).ge.0).and.(log_oi(logSize+1).le.gridnx).and.&
               (log_oj(logSize+1).ge.0).and.(log_oj(logSize+1).le.gridny).and.&
               (log_d(logSize+1).ge.1).and.(log_d(logSize+1).le.4).and.&
               (log_v(logSize+1).ge.0.0D0).and.(log_v(logSize+1).le.1.0D0)) then
             ! nothing outrageous in the log index, let's read the data
             inquire(unit=logunit, pos=startdatapos)
             CALL ASSERT((startdatapos-startpos).EQ.logEntrySize)
             if (log_seek(logSize+1).ge.0) then
                if (log_seek(logSize+1).ne.startdatapos) &
                     read(logunit, pos=log_seek(logSize+1))
                call Istate%load(logunit, ios) ! load the data
                if (log_seek(logSize+1).ne.startdatapos) &
                     read(unit=logunit, pos=startdatapos)
             end if
             if (ios.eq.0) then
                read(logunit, iostat=ios) magic
                fullyread = ((ios.eq.0).and.(magic.eq.logMagic))
             end if
          end if
          if (fullyread) then
             logSize = logSize + 1
             ic = log_oi(logSize)
             jc = log_oj(logSize)
             nd = log_d (logSize)
             v  = log_v (logSize)
             id = ic+diri(nd)
             jd = jc+dirj(nd)
             out_of_the_grid = (id.LT.0).OR.(id.GT.gridnx).OR. &
                  (jd.LT.0).OR.(jd.GT.gridny) ! if the next step is out
             if (have_callback) then
                xcbc=xg(ic)
                ycbc=yg(jc)
                if (.not.out_of_the_grid) then
                   xcbn=xg(id)
                   ycbn=yg(jd)
                   ! callback "trying"
                   ! notify that the point is "being computed"
                   CALL NEWCPCB(xcbn, ycbn , .false.)
                   if (v.LT.1.0D0) then
                      ! notify that the contour point is found
                      CALL NEWCPCB(xcbc+v*(xcbn-xcbc),ycbc+v*(ycbn-ycbc),.true.)
                   end if
                else
                   CALL NEWCPCB(xcbc, ycbc, .true.)
                end if
             end if
             ! mark the transition on the map
             map(nd,ic,jc)=v
             if ((v.lt.1.0D0).AND..NOT.out_of_the_grid) then
                ! mark the reverse jump as well to match the "check"
                map(revd(nd),id,jd) = 1-v
             end if
          else
             ! the log entry is broken, let's set it for rewrite
             read(logunit, pos=startpos)
             exit
          end if
       end do
       if (DEBUG_LOGFILE) &
            write (*, '("Loaded ",I0," entries. ")') logSize
       if (logSize.GT.0) then
          ! find last considered grid point, which is inside the boundary
          ! this will be the starting point for continued calculation
          do k=logSize,1,-1
             if (log_v(k).eq.1.0D0) then
                ld=log_d(k)
                ic=log_oi(k)+diri(ld)
                jc=log_oj(k)+dirj(ld)
                inquire(unit=logunit, pos=startdatapos)
                CALL ASSERT(log_seek(k).ge.0)
                read(logunit, pos=log_seek(k))
                call Istate%load(logunit, ios) ! load the data
                CALL ASSERT(ios.eq.0)
                read(logunit, pos=startdatapos)
                exit
             end if
          end do
       end if
    end if

    if (logSize.EQ.0) then         ! just starting or disk log is disabled
       if (PRESENT(logunit)) then
          ! if log is enabled, write the header
          rewind(logunit)
          call write_log_header
       end if
       call Istate%getxy(xc,yc)
       ic=nint((xc-gridx1)/gsx)
       jc=nint((yc-gridy1)/gsy)
       if ((ic.LT.0).OR.(ic.GT.gridnx).OR.(jc.LT.0).OR.(jc.GT.gridny)) then
          IERROR = 10
          return
       end if
       v = Istate%evolve(xg(ic), yg(jc))
       if (v.NE.1.0D0) then
          IERROR = 11
          return
       end if

       ! now we go down until we reach the boundary
       do
          if (logSize.GE.logMaxSize) then
             IERROR = 14
             return
          end if
          v = check(ic, jc, 1)
          if (v.lt.0.0D0.or.v.GT.1.0D0) then
             IERROR = 13
             return
          end if
          if (v.lt.1.0D0) exit
          ic = ic + diri(1)
          jc = jc + dirj(1)
       end do
       ld = ld0
    end if
    CALL ASSERT(logSize.GE.0)

    ! here is the main loop
    do
       neighbors=0
       if (map(1,ic,jc).ne.-1.0D0) neighbors = neighbors + 8
       if (map(2,ic,jc).ne.-1.0D0) neighbors = neighbors + 4
       if (map(3,ic,jc).ne.-1.0D0) neighbors = neighbors + 2
       if (map(4,ic,jc).ne.-1.0D0) neighbors = neighbors + 1
       nd = dirnext(1+neighbors,ld)

       if (nd.eq.0) exit

       if (logSize.GE.logMaxSize) then
          IERROR = 14
          return
       end if
       v = check(ic, jc, nd)
       if (v.lt.0.0D0.or.v.GT.1.0D0) then
          IERROR = 13
          return
       end if
       if (DEBUG_TRAVERSAL) then
          write (*,'("check at {",i0,",",i0,"} towards ",i0," yields ",g0.5)', &
               ADVANCE='NO') ic, jc, nd, v
       end if

       if (v.eq.1.0D0) then
          ic = ic + diri(nd)
          jc = jc + dirj(nd)
          ld = nd
       end if

       if (DEBUG_TRAVERSAL) then
          if (v.eq.1.0D0) then
             print '(" we go to {",i0,",",i0,"}.")', ic, jc
          else
             print '(" we stay at {",i0,",",i0,"}.")', ic, jc
          end if
       end if

    end do

    ! the contour is traced, let's copy it to the output array
    N = 0
    do k=1, logSize
       if (log_v(k).lt.1.0) then  ! a countour point
          N = N + 1
          if (N.gt.NMAX) then
             IERROR = 14
             return
          end if
          v=log_v(k)
          ic=log_oi(k)
          jc=log_oj(k)
          X(N) = xg(ic)+v*(xg(log_oi(k)+diri(log_d(k)))-xg(ic))
          Y(N) = yg(jc)+v*(yg(log_oj(k)+dirj(log_d(k)))-yg(jc))
       end if
    end do

    ! successfully finished
    IERROR = 0

  contains
    ! x coordinate of the vertical grid line i
    double precision function xg(i)
      integer, intent(IN) :: i
      CALL ASSERT((i.GE.0).AND.(i.LE.gridnx))
      xg = gridx1 + i * gsx
    end function xg

    ! y coordinate of the horizontal grid line j
    double precision function yg(j)
      integer, intent(IN) :: j
      CALL ASSERT((j.GE.0).AND.(j.LE.gridny))
      yg = gridy1 + j * gsy
    end function yg

    ! try to move from the grid point i, j to the direction d
    double precision function check(i, j, d)
      integer, intent(IN) :: i, j, d
      double precision xn, yn
      integer id, jd, k, logseek, currpos

      logical out_of_the_grid
      double precision xc, yc, xcn, ycn

      double precision xl, yl

      CALL ASSERT((i.GE.0).AND.(i.LE.gridnx))
      CALL ASSERT((j.GE.0).AND.(j.LE.gridny))
      CALL ASSERT((d.GE.1).AND.(d.LE.4))

      call Istate%getxy(xc,yc)
      ! check that the state is actually at the specified grid point
      ! (to a certain precision, of course)
      CALL ASSERT((DABS(xg(i)-xc).LE.1D-6).AND.(DABS(yg(j)-yc).LE.1D-6))

      id = i + diri(d)
      jd = j + dirj(d)
      out_of_the_grid = (id.LT.0).OR.(id.GT.gridnx).OR. &
           (jd.LT.0).OR.(jd.GT.gridny) ! if the next step is out of the grid
      logseek=-1
      if (PRESENT(logunit).AND..NOT.out_of_the_grid) then
         ! if log is enabled, let's first check if we have already traversed
         ! the destination point
         if (map(1,id,jd).ne.-1.0D0 .or. map(2,id,jd).ne.-1.0D0 .or. &
              map(3,id,jd).ne.-1.0D0 .or. map(4,id,jd).ne.-1.0D0 ) then
            ! if yes, let's find the stored state in the log
            if (DEBUG_LOGFILE) &
                 write (*,'("Looking for {",I0,",",I0,"} in the LOG. ")', &
                 ADVANCE="NO") id, jd
            do k=logSize, 1, -1
               if ((log_oi(k)+diri(log_d(k)).eq.id) .and. &
                    (log_oj(k)+dirj(log_d(k)).eq.jd) .and. &
                    (log_v(k).eq.1.0d0)) then
                  logseek=log_seek(k)
                  if (DEBUG_LOGFILE) &
                       write (*,'("Found at seek =", I0)', &
                       ADVANCE="NO") logseek
                  CALL ASSERT(logseek.gt.0)
                  ! load the state from the log
                  inquire(unit=logunit, pos=currpos)
                  read(logunit, pos=logseek)
                  call Istate%load(logunit, ios)
                  CALL ASSERT(ios.eq.0)
                  call Istate%getxy(xl, yl)
                  if (DEBUG_LOGFILE) &
                       write (*,'(" the with {x,y}={", G0.5,",",G0.5,"}.")',&
                       ADVANCE="NO") xl, yl
                  read(logunit, pos=currpos) ! restore the last log position
                  exit
               end if
            end do
            if (logseek.gt.0) then
               check = 1.0D0
            else
               if (DEBUG_LOGFILE) &
                    write (*,'(" NOT FOUND.")', ADVANCE="NO")
            end if
            if (DEBUG_LOGFILE) &
                 write (*,*)
         end if
      end if

      if (logseek.lt.0) then
         if (.not.out_of_the_grid) then
            ! then we try to evolve
            xn = xg(id)
            yn = yg(jd)
            if (have_callback) then
               ! notify that the point is being computed
               CALL NEWCPCB(xn, yn, .false.)
            end if
            check = Istate%evolve(xn,yn)
            if (check.lt.1.0D0) then  ! reached a boundary
               if (ASSERTS_ENABLED) then
                  ! check that computation, which did not reach the target,
                  ! also did not advance the state accidentally
                  call Istate%getxy(xcn,ycn)
                  CALL ASSERT((DABS(xc-xcn).LE.1D-6).AND.(DABS(yc-ycn).LE.1D-6))
               end if
               if (have_callback) then
                  ! notify that the boundary point is found
                  CALL NEWCPCB(xc+check*(xn-xc), yc+check*(yn-yc), .true.)
               end if
            end if
         else
            ! out of the grid point is not computed, but is a boundary
            check = 0.0D0
            if (have_callback) then
               ! notify that the boundary point is found
               CALL NEWCPCB(xc, yc, .true.)
            end if
         end if
      end if

      map(d,i,j)=check
      if ((check.lt.1.0D0).AND..NOT.out_of_the_grid) then
         ! in principle the topology must guarantee this, but
         ! numerical ODE solving does not always uphold this identity
         ! due to numerical errors
         map(revd(d),i+diri(d),j+dirj(d)) = 1-check
      end if

      ! update log in memory
      log_oi(logSize+1)=i
      log_oj(logSize+1)=j
      log_d (logSize+1)=d
      log_v(logSize+1)=check
      log_seek(logSize+1)=logseek
      logSize = logSize +1

      ! write log to disk, if enabled and if the final state was reached
      if (PRESENT(logunit)) then
         if ((logseek.lt.0).and.(check.eq.1.0D0)) then
            inquire(unit=logunit, pos=currpos)
            log_seek(logSize) = currpos + logEntrySize
         end if
         CALL ASSERT(((check.eq.1.0D0).and.(log_seek(logSize).ge.0)).or. &
              ((check.lt.1.0D0).and.(log_seek(logSize).lt.0)))
         if (DEBUG_LOGFILE) &
              print '("writing LOG entry ",I0, ":  {",I0,", ",I0,", ",I0,", ", &
              & G0.5,", ",I0,"}")', logSize, log_oi(logSize), log_oj(logSize), &
              log_d(logSize), log_v(logSize), log_seek(logSize)
         write (logunit, iostat=ios) log_oi(logSize), log_oj(logSize), &
              log_d(logSize), log_v(logSize), log_seek(logSize)
         if ((logseek.lt.0).and.(check.eq.1.0D0)) then
            call Istate%write(logunit, ios)
         end if
         write (logunit, iostat=ios) logMagic
         flush(logunit)
      end if
    end function check

    subroutine write_log_header()
      implicit none
      write(logunit, iostat=ios) logMagic
      CALL ASSERT(ios.EQ.0)
      write(logunit, iostat=ios) gridx1, gridy1, gridx2, gridy2, &
           gridnx, gridny, CCW, logMagic
      CALL ASSERT(ios.EQ.0)
    end subroutine write_log_header

    logical function log_header_matches_or_none(ios)
      implicit none
      integer, intent(OUT):: ios
      double precision :: Lgridx1, Lgridy1, Lgridx2, Lgridy2
      integer :: Lgridnx, Lgridny
      logical :: LCCW
      integer :: magic1, magic2

      log_header_matches_or_none = .FALSE.
      read(logunit, iostat=ios) magic1, Lgridx1, Lgridy1, Lgridx2, Lgridy2, &
           Lgridnx, Lgridny, LCCW, magic2
      if (ios.NE.0) then
         ! actually an I/O error is OK here, it would just mean that
         ! the LOG is going to be overwritten
         log_header_matches_or_none = .TRUE.
         return
      end if
      if ( &
           (gridx1.ne.Lgridx1).or.(gridx2.ne.Lgridx2).or. &
           (gridy1.ne.Lgridy1).or.(gridy2.ne.Lgridy2).or. &
           (gridnx.ne.Lgridnx).or.(gridny.ne.Lgridny).or. &
           (CCW.neqv.LCCW) &
           ) return
      if ((magic1.ne.logMagic).or.(magic2.ne.logMagic)) return
      log_header_matches_or_none = .TRUE.

    end function log_header_matches_or_none
  end subroutine TRACE

end module tracer
