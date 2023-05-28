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
!     A STATE CLASS IMPLEMENTATION FOR HOPFION

module hopf_state
  use skyrme_state
  use tracer
  implicit none

  type, extends(skyrmion_state) :: hopfion_state
     integer :: HTYPE
     double precision :: g
     logical :: gEqu
   contains
     ! superclass methods overrides
     procedure :: init => hs_create
     procedure :: initHopfion => hs_initHopfion
     procedure :: bisectSaved => hs_bisectSaved
     procedure :: write => hs_write
     procedure :: load => hs_load
     ! setting and querying predicates
     procedure :: predicate => hs_predicate
     procedure :: predicateMax => hs_predicateMax
     procedure :: predicateInfo => hs_predicateInfo
     ! energy computation
     procedure :: energy => hs_energy
     ! misc
     procedure :: write_nml => hs_write_nml
     procedure :: load_nml => hs_load_nml
     procedure :: get_short_code => hs_get_short_code
     procedure :: print => hs_print
  end type hopfion_state

contains
  ! a separate assignment operator
  subroutine hs_set_from(to_obj,from_obj)
    implicit none
    class(hopfion_state), intent(out) :: to_obj
    class(hopfion_state), intent(in) :: from_obj
    call sk_set_from(to_obj, from_obj)
    to_obj%HTYPE=from_obj%HTYPE
    to_obj%g=from_obj%g
    to_obj%gEqu=from_obj%gEqu
  end subroutine hs_set_from

  ! this would create a T1 hopfion state at q=xs, h=ys, gamma=1 and mu=0
  subroutine hs_create(this, xs, ys)
    use profile
    implicit none
    class(hopfion_state), intent(OUT) :: this
    double precision, intent(IN) :: xs, ys
    call this%initHopfion(xs, ys, 1, .false.)
  end subroutine hs_create

  ! This would create a T1/T2 hopfion state at q, h, gamma=1 and mu=0
  ! note it is only guaranteed to succeed in case q=h=0, otherwise
  ! it might fail assertion. The suggested approach is to create an
  ! q=h=0 hopfion and then evolve it towards the desired point on the
  ! phase diagram.
  subroutine hs_initHopfion(this, q, h, type, gEqu)
    use debug
    use profile
    implicit none
    class(hopfion_state), intent(OUT) :: this
    double precision, intent(IN) :: h, q
    integer, intent(IN) :: type
    logical, intent(IN) :: gEqu

    CALL ASSERT(.not. gEqu)    ! not yet implemented
    CALL ASSERT(type.eq.1.or.type.eq.2)

    call state_create(this, q, h)

    this%HTYPE = type
    this%mu    = 0.0D0
    this%g     = 1.0D0
    this%gEqu  = gEqu
    this%selectedPredicate = 0
    this%NFPTS = 0
    select case(type)
    case (1)
       this%nu    =-0.13D0
       call SOLVE_EULER(this%HTYPE, this%g, h, q, this%nu, this%mu, &
            NFPTSMAX, this%NFPTS, this%FPTS)
    case (2)
       this%nu    = 0.1D0
       call SOLVE_EULER(this%HTYPE, this%g, h, q, this%nu, this%mu, &
            NFPTSMAX, this%NFPTS, this%FPTS)
    case default
       CALL ASSERT(.false.)
    end select
    call ASSERT(this%NFPTS.GE.2)
  end subroutine hs_initHopfion

  ! Does one bisection step or a full bisection.
  recursive double precision function hs_bisectSaved(this, q, h, mu, &
       report, onestep)
    use profile
    class(hopfion_state), intent(INOUT) :: this
    double precision, intent(IN) :: q, h, mu
    logical, intent(IN) :: report, onestep
    ! update the type of the "old" variable in subclasses
    type(hopfion_state) old
    integer INFO

    call hs_set_from(old, this)
    if (onestep) then
       ! --------------------------------------------------------------------
       ! subclasses need to rewrite this part with the rest mostly unchanged
       !
       hs_bisectSaved = 0.0D0 ! fail by default
       INFO = 0
       call ITERATE_MS(this%HTYPE, this%g, h, q, this%nu, &
            mu, NFPTSMAX, this%NFPTS, this%FPTS, INFO)
       if (this%NFPTS.GE.2 .and. this%predicate()) hs_bisectSaved = 1.0D0
       !
       !---------------------------------------------------------------------
    else
       hs_bisectSaved = this%bisectTowards(q, h, mu, report)
    end if
    if (hs_bisectSaved .eq. 1.0D0) then
       ! successfully evolved
       this%mu=mu
       hs_bisectSaved = state_evolve(this, q, h)
    else
       ! attempt failed, restore the old profile
       call hs_set_from(this, old)
    end if
  end function hs_bisectSaved

  subroutine hs_write(this, unit, ios)
    implicit none
    class(hopfion_state), intent(IN) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios
    call sk_write(this, unit, ios)
    if (ios /= 0) return
    write(unit,iostat=ios) this%HTYPE, this%g, this%gEqu
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine hs_write

  ! state loader
  subroutine hs_load(this, unit, ios)
    class(hopfion_state), intent(OUT) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios
    call sk_load(this, unit, ios)
    if (ios /= 0) return
    read(unit,iostat=ios) this%HTYPE, this%g, this%gEqu
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine hs_load

  logical function hs_predicate(this)
    use debug
    implicit none
    class(hopfion_state), intent(IN) :: this
    integer i
    select case(this%selectedPredicate)
    case (0)
       hs_predicate = .true.
    case (1)
       hs_predicate = .true.
       do i=2, this%NFPTS
          if (this%FPTS(i,2)-this%FPTS(i-1,2).LT.-1D-5) then
             hs_predicate = .false.
             exit
          end if
       end do
    case default
       CALL ASSERT(.false.)
    end select
  end function hs_predicate

  integer function hs_predicateMax(this)
    implicit none
    class(hopfion_state), intent(IN) :: this
    hs_predicateMax = 1
  end function hs_predicateMax

  subroutine hs_predicateInfo(this, pred, buflen, buf, bufShort)
    use debug
    implicit none
    class(hopfion_state), intent(IN) :: this
    integer, intent(IN):: pred, buflen
    character(len=buflen), intent(out) :: buf, bufShort
    integer p
    p = pred
    if (pred.eq.-1) p = this%selectedPredicate
    select case(p)
    case (0)
       buf = "No predicate."
       bufShort = "NO PREDICATE"
    case (1)
       buf = "Monotonous f(r)."
       bufShort = "MONOTONOUS F"
    case default
       CALL ASSERT(.false.)
    end select
  end subroutine hs_predicateInfo

  double precision function hs_energy(this)
    use energy
    implicit none
    class(hopfion_state), intent(IN) :: this
    double precision q, h
    call this%getxy(q,h)
    hs_energy = Etot(this%HTYPE, 4.0D0*DSQRT(2.0D0) , this%g, h, q, &
         this%nu, this%mu, NFPTSMAX, this%NFPTS, this%FPTS)
  end function hs_energy

  subroutine hs_write_nml(this, unit, ios)
    implicit none
    class(hopfion_state), intent(IN) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    integer HTYPE
    double precision g
    logical gEqu

    namelist /HOPFION/ HTYPE, g, gEqu

    HTYPE = this%HTYPE
    g  = this%g
    gEqu = this%gEqu
    write (unit=unit, nml=HOPFION, iostat=ios)
    if (ios /= 0) return
    call sk_write_nml(this, unit, ios)
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine hs_write_nml

  subroutine hs_load_nml(this, unit, ios)
    class(hopfion_state), intent(OUT) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    integer HTYPE
    double precision g
    logical gEqu
    namelist /HOPFION/ HTYPE, g, gEqu
    read (unit=unit, nml=HOPFION, iostat=ios)
    if (ios /= 0) return
    call sk_load_nml(this, unit, ios)
    this%HTYPE = HTYPE
    this%g = g
    this%gEqu = gEqu
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine hs_load_nml

  subroutine hs_get_short_code(this, code)
    class(hopfion_state), intent(IN) :: this
    character (len=2), intent(OUT) :: code

    select case(this%HTYPE)
    case (1)
       code = "H1"
    case (2)
       code = "H2"
    case default
       CALL ASSERT(.false.)
    end select
  end subroutine hs_get_short_code

  subroutine hs_print(this)
    implicit none
    class(hopfion_state), intent(IN) :: this
    print '("HTYPE=",I0,"; g=",G0.5,"; gEqu=",L1)', &
         this%HTYPE, this%g, this%gEqu
    call sk_print(this)
  end subroutine hs_print

end module hopf_state
