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
!     HERE A SIMPLE STATE INSTANCE IS DEFINED FOR TESTING THE TRACER

module test_state
  use tracer
  implicit none

  type, extends(state) :: testing_state
   contains
     procedure :: init => ts_create
     procedure :: evolve => ts_evolve
  end type testing_state

contains
  subroutine ts_create(this, xs, ys)
    class(testing_state), intent(OUT) :: this
    double precision, intent(IN) :: xs, ys
    call state_create(this, xs, ys)
  end subroutine ts_create

  double precision function ts_evolve(this, xt, yt)
    implicit none
    class(testing_state), intent(INOUT) :: this
    double precision, intent(IN) :: xt, yt
    double precision :: xo, yo
    double precision zeroin

    if (f(xt, yt).LE.0.0D0) then
       ts_evolve = state_evolve(this, xt, yt)
    else
       call this % getxy(xo, yo)
       ts_evolve = zeroin(0.0d0,1.0d0, equ, 1d-6)
    end if
    ! print '("evolve(",g0.5,",",g0.5,")=",g0.5)', xt, yt, ts_evolve

  contains

    double precision function equ(t)
      implicit none
      double precision t
      equ=f( xo + t*(xt - xo), yo + t*(yt - yo) )
    end function equ

    double precision function f(x, y)
      implicit none
      double precision, intent(IN) :: x, y
      f = &
           (DSQRT((2.0D0 + x)**2/4.0D0 + (-3.0D0 + y)**2)/2.0D0 + &
           DSQRT(3.0D0*(-1.0D0+x)**2 + (7.0D0*(-3.0D0 + y)**2)/3.0D0)/4.0D0)* &
           (DSQRT((2.0D0 + x)**2/4.0D0 + ( 3.0D0 + y)**2)/2.0D0 + &
           DSQRT(3.0D0*(-1.0D0+x)**2 + (7.0D0*( 3.0D0 + y)**2)/3.0D0)/4.0D0) &
           -8.0D0
    end function f

  end function ts_evolve

end module test_state
