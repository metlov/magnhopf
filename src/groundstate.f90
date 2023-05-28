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
!     Functions in this module compute the energy of the simple 1D and 0D
!     states of helimagnet -- such as uniform, conical and helicoidal.
!     There is also function interfacing with skyrmion energy minimization
!     routines, for computing the classical phase diagram of the
!     helimagnet (including the skyrmions) and the corresponding ground
!     state energy an every point of the phase diagram.
module groundstate
  use debug
  implicit none

  public :: hasStableHelices

contains

  ! the critical q for small 0 < |h| < Pi^2/16
  double precision function QCRHSMALL(h)
    implicit none
    double precision, intent(IN) :: h
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    double precision :: c
    double precision :: zeroin
    CALL ASSERT((DABS(h).GE.0.0D0).AND.(DABS(h).LE.(pi*pi/16.0D0)))
    ! evaluate analytically at endpoints
    if (DABS(h).LE.2.0D0*EPSILON(h)) THEN
       QCRHSMALL = pi*pi/4.0D0
       RETURN
    end if
    if (DABS(h-pi*pi/16.0D0).LE.2.0D0*EPSILON(h)) THEN
       QCRHSMALL = 0.0D0
       RETURN
    end if
    c = zeroin(2.D0*EPSILON(c),1.0d0-2.D0*EPSILON(c), QCRHSMALL_equ, 1d-6)
    QCRHSMALL = (c**2*DABS(h))/(1.0D0 - c**2)

  contains

    double precision function DATANH(X)
      implicit none
      double precision,intent(IN) :: X
      DATANH = 0.5D0 * DLOG ( (1.0D0+X) / (1.0D0-X) )
    end function DATANH

    double precision function QCRHSMALL_equ(c)
      implicit none
      double precision,intent(IN) :: c
      QCRHSMALL_equ= -DABS(h) + ((c**2 - c**4)*pi**2)/ &
           (4.0D0*(-c - DATANH(c) + c**2*DATANH(c))**2)
    end function QCRHSMALL_equ
  end function QCRHSMALL

   ! the critical q for large Pi^2/16 < h < 1
  double precision function QCRHLARGE(h)
    implicit none
    double precision, intent(IN) :: h
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    double precision :: c
    double precision :: zeroin
    CALL ASSERT((DABS(h).GE.(pi*pi/16.0D0)).AND.(DABS(h).LE.1.0D0))
    ! evaluate analytically at endpoints
    if (DABS(h-pi*pi/16.0D0).LE.2.0D0*EPSILON(h)) THEN
       QCRHLARGE = 0.0D0
       RETURN
    end if
    if (DABS(DABS(h)-1.0D0).LE.2.0D0*EPSILON(h)) THEN
       QCRHLARGE = -1.0D0
       RETURN
    end if
    c = zeroin(2.D0*EPSILON(c),1.0d0-2.D0*EPSILON(c), QCRHLARGE_equ, 1d-6)
    QCRHLARGE = (-1.0D0 + c**2)*DABS(h)

  contains

    double precision function QCRHLARGE_equ(c)
      implicit none
      double precision,intent(IN) :: c
      QCRHLARGE_equ=-DABS(h) - ((-1.0D0 + c**2)*pi**2)/ &
           (4.0D0*(c**2 - c**4 + DASIN(DSQRT(1.0D0 - c**2))**2 + &
           2.0D0*c*DASIN(DSQRT(1.0D0 - c**2))*DSQRT(1.0D0 - c**2)))
    end function QCRHLARGE_equ

  end function QCRHLARGE

  logical function hasStableHelices(q, h)
    implicit none
    double precision, intent(IN) :: q, h
    double precision, parameter :: eps = 1D-4
    double precision, parameter :: pisqby4 = 2.46740110027233965470862274997d0
    double precision, parameter :: pisqby16 = 0.616850275068084913677155687492d0

    if ((q.LT.-1.0D0).OR.(q.GT.pisqby4).OR.(h.GT.1.0D0).OR.(h.LT.-1.0D0)) then
       ! definitely no helices
       hasStableHelices = .false.
    else
       if (((q.eq.-1.0D0).AND.(DABS(h).EQ.1.0D0)).OR. &
            ((q.eq.0.0D0).AND.(h.EQ.0.0D0)) .OR. &
            ((q.eq.pisqby4).AND.(h.EQ.0.0D0)) .OR. &
            ((q.eq.0.0D0).AND.(DABS(h).EQ.pisqby16)) .OR. &
            ((DABS(h).LE.eps).AND.((q.GE.0.0D0).AND.(q.LE.pisqby4))) ) then
          ! some corner points
          hasStableHelices = .true.
       else
          if ((q.LT.0.0D0).AND.(DABS(h).LT.DSQRT(DABS(q)))) then
             ! inside the horizontal parabola on the left
             hasStableHelices = .false.
          else
             if (DABS(h).GT.pisqby16) then
                if (q.GT.0.0D0) then
                   ! definitely no helices, no need to solve the eq
                   hasStableHelices = .false.
                else
                   hasStableHelices = (q.LT.QCRHLARGE(DABS(h)))
                end if
             else
                if (q.LT.0.0D0) then
                   ! definitely there are helices, no need to solve the eq
                   hasStableHelices = .true.
                else
                   hasStableHelices = (q.LT.QCRHSMALL(DABS(h)))
                end if
             end if
          end if
       end if
    end if
  end function hasStableHelices


  ! The equilibrium energy of the helical state at given q and h
  ! The q and h values must correspond to hasStableHelices(q, h) being .true.
  double precision function energyHelix(q, h)
    implicit none
    double precision, intent(IN) :: q, h
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    integer, parameter :: NFPTSMAX = 1000

    ! parameters to the MUSN call
    integer, parameter :: N = 2
    double precision A, B, ER(5), AMP
    integer NRTI
    integer, parameter :: ITLIM = 30
    double precision TI(NFPTSMAX), Y(N,NFPTSMAX), QQ(N,N,NFPTSMAX)
    integer, parameter :: NUDIM = 3 ! N*(N+1)/2
    double precision U(NUDIM,NFPTSMAX), PHIREC(NUDIM,NFPTSMAX)
    double precision D(N,NFPTSMAX)
    integer KPART, LW, LIW, LWG
    double precision W(7*N + 3*N*NFPTSMAX + 4*N*N  )
    integer IW(3*N + NFPTSMAX)
    double precision WGR(NFPTSMAX*6)
    integer IERROR

    integer I
    double precision C1SUM

    CALL ASSERT(hasStableHelices(q, h))

    ! some corner values
    if ((q.EQ.0.0D0).AND.(h.EQ.0.0D0)) then
       energyHelix = -0.5D0
       return
    end if

    A      = 0.0D0
    B      = 2.0D0*pi
    ER(1)  = 1.0D-6
    ER(2)  = 1.0D-6
    ER(3)  = 1.1D-15
    NRTI   = 120
    AMP    = 0.80D0
    LW     = size(W)
    LIW    = size(IW)
    LWG    = size(WGR)
    IERROR = 0     ! set to 1 for making MUSL print diagnostics
    call MUSN(EQUATIONS, INITIAL_APPROX, BOUNDARY, N, &
         A, B, ER, TI, NFPTSMAX, NRTI, AMP, ITLIM, &
         Y, QQ, U, NUDIM, D, &
         PHIREC, KPART, W, LW, IW, LIW, WGR, LWG, IERROR)
    if (IERROR.eq.0) then
       C1SUM=0.0D0
       do I=1,NRTI
          C1SUM=C1SUM+Y(2,I)
       end do
       energyHelix = -C1SUM/NRTI
    else
       call xerror("MUSN RETURNED IERROR",IERROR,0)
       CALL ASSERT(.false.)
       energyHelix = 0.0D0
    end if
  contains
    subroutine EQUATIONS(T, Y, FUN)
      implicit none
      double precision T, Y(2), FUN(2)
      double precision rootArg

      rootArg = -2*h*DCOS(t) - q*DCOS(t)**2 + 2*Y(2)
      if (rootArg.GE.0.0D0) then
         FUN(1) = DSQRT(rootArg)
      else
         FUN(1) = 0.0D0
      end if
      FUN(2) = 0.0D0
    end subroutine EQUATIONS
  subroutine BOUNDARY(N,YA,YB,FG,DGA,DGB)
    implicit none
    integer N
    double precision YA(N),YB(N),FG(N),DGA(N,N),DGB(N,N)

    CALL ASSERT(N.eq.2)
    FG(1) = YA(1)                   ! f(A) = 0
    FG(2) = YB(1) - 2.0D0*pi        ! f(B) = 2*pi

    DGA=0.0D0
    DGB=0.0D0

    DGA(1,1) = 1.D0
    DGB(2,1) = 1.D0
    return
  end subroutine BOUNDARY
  subroutine INITIAL_APPROX(T,Y)
    implicit none
    double precision, intent(IN) :: T
    double precision, intent(OUT) :: Y(2)
    double precision v1, v2, v3
    integer np

    Y(1) = 2.0D0*pi*DSIN(T/4.0D0)**2

    ! some ugly black magic for setting the initial C1
    np = 3
    v1 =  h + q/2
    call check_v(v1,np)
    v2 = -h + q/2
    call check_v(v2,np)
    if (ABS(q).GT.0.0D0) then
       v3 = -h**2/(2.0D0*q)
    else
       v3 = 1.5D0  ! an out-of-range value
    end if
    call check_v(v3,np)
    if (np.eq.0) then
       Y(2) = (0.5D0 + pi * pi/ 8.0D0)*0.5D0
    else
       if (q.GT.0.0D0) then
          Y(2) = MAX (v1, v2, v3)
       else
          Y(2) = MIN (v1, v2, v3)
       end if
    end if
    ! print '("Y1(",G0.5,")=",G0.5," Y2(",G0.5,")=",G0.5)', T, Y(1), T, Y(2)
  end subroutine INITIAL_APPROX
  subroutine check_v(v, np)
    double precision, intent(INOUT) :: v
    integer, intent(INOUT) :: np

    if (.NOT.((v.GE.0.5D0).AND.(v.LE.(pi*pi/8.0D0)))) then
       ! assign a value, which will definitely be excluded in
       ! a later selection
       if (q.GT.0.0D0) then
          v = 0.0D0
       else
          v = 1.5D0
       end if
       ! note that this value was out of range
       np = np - 1
    end if
  end subroutine check_v
  end function energyHelix

  logical function isSkXGroundState(q, h, mu)
    implicit none
    double precision, intent(IN) :: q, h, mu
    include "skyrmion_ground.inc"  ! load pre-computed contour
    integer INOUT

    CALL ASSERT(SIZE(data_q).EQ.SIZE(data_h))
    CALL PNPOLY(q, h, data_q, data_h, SIZE(data_q), INOUT)
    isSkXGroundState = (INOUT.EQ.1)
  end function isSkXGroundState

  double precision function energySkX(q, h, mu)
    use skyrmions
    implicit none
    double precision, intent(IN) :: q, h, mu
    integer, parameter :: NFPTSMAX = 1000
    integer :: NFPTS
    double precision :: nu, FPTS(NFPTSMAX,3)

    CALL ASSERT(isSkXGroundState(q, h, mu))
    nu    =-0.13D0
    NFPTS = 0
    call SOLVE_EULER_SK(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    CALL ASSERT(NFPTS.GE.2)

    energySkX = EtotSk(h, q, nu, NFPTSMAX, NFPTS, FPTS)
  end function energySkX

  ! determines the ground state type at a particular point of the phase diagram
  ! q, h, mu        -- coordinates of the point
  ! state           -- the type of the ground state, takes the values:
  !                   -1 -- unknown (or all masked, in which case stateE=1D10)
  !                    1 -- uniform
  !                    2 -- cones
  !                    3 -- helices
  !                    4 -- skyrmions
  ! stateE           -- the energy of the ground state
  ! statesToConsider -- an optional logical array, allowing to mask particular
  !                     states from consideration (by default all states
  !                     are considered)
  subroutine getGroundState(q, h, mu, state, stateE, statesToConsider)
    double precision, intent(IN) :: q, h, mu
    integer, intent(OUT) :: state
    double precision, intent(OUT) :: stateE
    logical, intent(IN), optional :: statesToConsider(4)
    logical sTC(4)
    double precision currE

    sTC=.true.
    if (PRESENT(statesToConsider)) sTC=statesToConsider

    state  = -1
    stateE = 1D10

    if (sTC(1)) then
       state = 1
       if ((h**2).LT.(q**2)) then
          stateE = MIN(-h - q/2.0D0, h - q/2.0D0, h**2/(2.0D0 * q))
       else
          stateE = MIN(-h - q/2.0D0, h - q/2.0D0)
       end if
    end if

    if (sTC(2)) then
       if ((h.GT.MIN(-1.0D0 + q, 1.0D0 - q)) .AND. &
            (h.LT.MAX(-1.0D0 + q, 1.0D0 - q))) then
          currE = (-1.0D0 - h**2 + q)/(2.0D0*(1.0D0 - q))
          if (currE.LT.stateE) then
             state = 2
             stateE = currE
          end if
       end if
    end if

    if (sTC(3)) then
       if (hasStableHelices(q, h)) then
          currE = energyHelix(q, h)
          if (currE.LT.stateE) then
             state = 3
             stateE = currE
          end if
       end if
    end if

    if (sTC(4)) then
       if (isSkXGroundState(q, h, mu)) then
          currE = energySkX(q, h, mu)
          if (currE.LT.stateE) then
             state = 4
             stateE = currE
          end if
       end if
    end if
  end subroutine getGroundState

  subroutine getStateSymbol(state, symb)
    integer, intent(IN) :: state
    character(len=3), intent(OUT) :: symb

    select case(state)
    case (1)
       symb="UNI"
    case (2)
       symb="CON"
    case (3)
       symb="HEL"
    case (4)
       symb="SkX"
    case default
       CALL ASSERT(.false.)
    end select
  end subroutine getStateSymbol

end module groundstate
