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
!     A STATE CLASS IMPLEMENTATION FOR SKYRMION
module skyrme_state
  use tracer
  use skyrmions
  implicit none

  integer, parameter:: NFPTSMAX = 1000

  type, extends(state) :: skyrmion_state
     integer :: NFPTS
     double precision :: FPTS(NFPTSMAX,3)
     double precision :: nu, mu
     integer :: selectedPredicate
   contains
     ! superclass methods overrides
     procedure :: init => sk_init
     procedure :: initHopfion => sk_initHopfion
     procedure :: evolve => sk_evolve
     procedure :: write => sk_write
     procedure :: load => sk_load
     ! state evolution on the h-q-mu plane
     procedure :: evolve3D => sk_evolve3D
     procedure :: bisectTowards => sk_bisectTowards
     procedure :: bisectSaved => sk_bisectSaved
     ! setting and querying predicates
     procedure :: predicate => sk_predicate
     procedure :: predicateMax => sk_predicateMax
     procedure :: predicateSet => sk_predicateSet
     procedure :: predicateInfo => sk_predicateInfo
     ! energy computation
     procedure :: energy => sk_energy
     ! misc
     procedure :: write_nml => sk_write_nml
     procedure :: load_nml => sk_load_nml
     procedure :: get_short_code => sk_get_short_code
     procedure :: print => sk_print
  end type skyrmion_state

contains
  ! a separate assignment operator
  subroutine sk_set_from(to_obj,from_obj)
    implicit none
    class(skyrmion_state), intent(out) :: to_obj
    class(skyrmion_state), intent(in) :: from_obj
    double precision xc, yc, dummy
    integer I
    call from_obj%getxy(xc,yc)
    dummy = state_evolve(to_obj,xc,yc)
    to_obj%NFPTS=from_obj%NFPTS
    do I=1,to_obj%NFPTS
       to_obj%FPTS(I,1)=from_obj%FPTS(I,1)
       to_obj%FPTS(I,2)=from_obj%FPTS(I,2)
       to_obj%FPTS(I,3)=from_obj%FPTS(I,3)
    end do
    to_obj%mu=from_obj%mu
    to_obj%nu=from_obj%nu
    to_obj%selectedPredicate=from_obj%selectedPredicate
  end subroutine sk_set_from

  ! this would create a T1 hopfion state at q=xs, h=ys, gamma=1 and mu=0
  subroutine sk_init(this, xs, ys)
    use profile
    implicit none
    class(skyrmion_state), intent(OUT) :: this
    double precision, intent(IN) :: xs, ys
    double precision  h, q
    call state_create(this, xs, ys)
    q     = xs
    h     = ys
    this%nu    =-0.13D0
    this%mu    = 0.0D0
    this%NFPTS = 0
    call SOLVE_EULER_SK(h, q, this%nu, NFPTSMAX, this%NFPTS, this%FPTS)
    this%selectedPredicate = 0
  end subroutine sk_init

  ! an entry point to be able to initialize the hopfion via the skyrmion
  ! interface (have to do it due to the absence of casts)
  subroutine sk_initHopfion(this, q, h, type, gEqu)
    use debug
    implicit none
    class(skyrmion_state), intent(OUT) :: this
    double precision, intent(IN) :: h, q
    integer, intent(IN) :: type
    logical, intent(IN) :: gEqu
    call ASSERT(.false.)
  end subroutine sk_initHopfion

  ! this computes evolution in the q-h plane at a fixed mu
  ! sk_evolve3D allows motion in the whole q-h-my space.
  double precision function sk_evolve(this, xt, yt)
    implicit none
    class(skyrmion_state), intent(INOUT) :: this
    double precision, intent(IN) :: xt, yt
    sk_evolve = this%bisectSaved(xt, yt, this%mu, .false., .false.)
  end function sk_evolve

  subroutine sk_write(this, unit, ios)
    implicit none
    class(skyrmion_state), intent(IN) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios
    integer I

    call state_write(this, unit, ios)
    if (ios /= 0) return
    write(unit,iostat=ios) this%mu, this%nu, this%selectedPredicate, this%NFPTS
    if (ios /= 0) return
    write(unit,iostat=ios) (this%FPTS(I,1), I=1, this%NFPTS)
    if (ios /= 0) return
    write(unit,iostat=ios) (this%FPTS(I,2), I=1, this%NFPTS)
    if (ios /= 0) return
    write(unit,iostat=ios) (this%FPTS(I,3), I=1, this%NFPTS)
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine sk_write

  ! state loader
  subroutine sk_load(this, unit, ios)
    class(skyrmion_state), intent(OUT) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios
    integer I

    call state_load(this, unit, ios)
    if (ios /= 0) return
    read(unit,iostat=ios) this%mu, this%nu, this%selectedPredicate, this%NFPTS
    if (ios /= 0) return
    read(unit,iostat=ios) (this%FPTS(I,1), I=1, this%NFPTS)
    if (ios /= 0) return
    read(unit,iostat=ios) (this%FPTS(I,2), I=1, this%NFPTS)
    if (ios /= 0) return
    read(unit,iostat=ios) (this%FPTS(I,3), I=1, this%NFPTS)
    if (ios /= 0) return ! just in case something will be added after this stmt
  end subroutine sk_load

  ! Computes evolution in the q-h-mu plane.
  ! It needs to be overridden in each sublass to manage the state save.
  ! The actual bisection in sk_bisectTowards is done polymorphically.
  double precision function sk_evolve3D(this, q, h, mu, report)
    implicit none
    class(skyrmion_state), intent(INOUT) :: this
    double precision, intent(IN) :: q, h, mu
    logical, intent(IN) :: report
    sk_evolve3D = this%bisectSaved(q, h, mu, report, .false.)
  end function sk_evolve3D

  ! This function is polymorphic and should work in subclasses as well
  double precision function sk_bisectTowards(this, qt, ht, mut, report)
    class(skyrmion_state), intent(INOUT) :: this
    double precision, intent(IN) :: qt, ht, mut
    logical, intent(IN) :: report

    ! accuracy of the critical line probing
    double precision, parameter :: minstep = 1E-6
    double precision t, step , qs, hs, mus, q, h, mu

    mus=this%mu
    call this%getxy(qs,hs)
    step = 1.0D0
    t = 0.0D0
    do while((t.LT.1.0D0).AND.(step.GE.minstep))
       if (t + step.GT.1.0D0) step=1.0D0-t
       q  =  qs + (t + step)*( qt -  qs)
       h  =  hs + (t + step)*( ht -  hs)
       mu = mus + (t + step)*(mut - mus)
       if (report) write(6, '("trying {",G0.5,",", G0.5, ",", G0.5,"}")', &
            ADVANCE='no')  q, h, mu
       if (this%bisectSaved(q, h, mu, report, .true. ) .eq. 1.0D0) then
          if (report) write(6, '(" got nu=",G0.5," .")') this%nu
          t = t + step
       else
          if (report) write(6, '(" missed.")')
          step=step/2.0D0
       end if
    end do
    sk_bisectTowards = t
  end function sk_bisectTowards

  ! Does one bisection step or a full bisection.
  !
  ! It also saves/restores the state (which can't be done polymorphically),
  ! so this function must be overridden in subclasses.
  !
  ! The idea behind the onestep parameter is to have a single non-polymorphic
  ! function in the class, which needs a certain override in subclasses
  ! (even though some small amount of code will be repeated), but to
  ! have the rest of the evolution code generic.
  recursive double precision function sk_bisectSaved(this, q, h, mu, &
       report, onestep)
    class(skyrmion_state), intent(INOUT) :: this
    double precision, intent(IN) :: q, h, mu
    logical, intent(IN) :: report, onestep
    ! update the type of the "old" variable in subclasses
    type(skyrmion_state) old

    call sk_set_from(old, this)
    if (onestep) then
       ! --------------------------------------------------------------------
       ! subclasses need to rewrite this part with the rest mostly unchanged
       !
       sk_bisectSaved = 0.0D0 ! fail by default
       call SOLVE_EULER_SK(h, q, this%nu, NFPTSMAX, this%NFPTS, this%FPTS)
       if (this%NFPTS.GE.2 .and. this%predicate()) sk_bisectSaved = 1.0D0
       !
       !---------------------------------------------------------------------
    else
       sk_bisectSaved = this%bisectTowards(q, h, mu, report)
    end if
    if (sk_bisectSaved .eq. 1.0D0) then
       ! successfully evolved
       this%mu=mu
       sk_bisectSaved = state_evolve(this, q, h)
    else
       ! attempt failed, restore the old profile
       call sk_set_from(this, old)
    end if
  end function sk_bisectSaved

  logical function sk_predicate(this)
    use debug
    use groundstate
    implicit none
    class(skyrmion_state), intent(IN) :: this
    logical statesToConsider(4)
    integer stateT
    double precision h, q, stateE

    select case(this%selectedPredicate)
    case (0)
       sk_predicate = .true.
    case (1)
       statesToConsider=.true.
       statesToConsider(4)=.false. ! except skyrmions themselves
       call this%getxy(q,h)
       call getGroundState(q, h, this%mu, stateT, stateE, statesToConsider)
       sk_predicate = (this%energy().LT.stateE)
    case default
       CALL ASSERT(.false.)
    end select
  end function sk_predicate

  integer function sk_predicateMax(this)
    implicit none
    class(skyrmion_state), intent(IN) :: this

    sk_predicateMax = 1
  end function sk_predicateMax

  subroutine sk_predicateInfo(this, pred, buflen, buf, bufShort)
    use debug
    implicit none
    class(skyrmion_state), intent(IN) :: this
    integer, intent(IN):: pred, buflen
    character(len=buflen), intent(out) :: buf, bufShort
    integer p
    if (pred.eq.-1) then
       p = this%selectedPredicate
    else
       p = pred
    end if
    select case(p)
    case (0)
       buf = "No predicate."
       bufShort = "NO PREDICATE"
    case (1)
       buf = "Skyrmion ground state."
       bufShort = "SkX GROUND S"
    case default
       CALL ASSERT(.false.)
    end select
  end subroutine sk_predicateInfo

  ! sets the predicate on the current object
  ! this function is generic and should not require override in subclasses
  logical function sk_predicateSet(this, pred)
    implicit none
    class(skyrmion_state), intent(INOUT) :: this
    integer, intent(IN):: pred
    integer oldpred
    if (pred.GE.0.AND.pred.LE.this%predicateMax()) then
       oldpred = this%selectedPredicate
       this%selectedPredicate=pred
       sk_predicateSet = this%predicate()
       if (.not.sk_predicateSet) this%selectedPredicate=oldpred
    else
       sk_predicateSet = .false.
    end if
  end function sk_predicateSet

  double precision function sk_energy(this)
    use energy
    implicit none
    class(skyrmion_state), intent(IN) :: this
    double precision q, h
    call this%getxy(q,h)
    sk_energy = EtotSk(h, q, this%nu, NFPTSMAX, this%NFPTS, this%FPTS)
  end function sk_energy

  subroutine sk_write_nml(this, unit, ios)
    implicit none
    class(skyrmion_state), intent(IN) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    double precision q, h, nu, mu
    integer selectedPredicate
    integer NFPTS
    double precision, allocatable :: FPTS1(:), FPTS2(:), FPTS3(:)

    namelist /SKYRMION/ q, h, nu, mu, selectedPredicate, &
         NFPTS, FPTS1, FPTS2, FPTS3

    call this%getxy(q,h)
    mu = this%mu
    nu = this%nu
    selectedPredicate = this%selectedPredicate
    NFPTS = this%NFPTS
    allocate(FPTS1(NFPTS))
    allocate(FPTS2(NFPTS))
    allocate(FPTS3(NFPTS))
    FPTS1 = this%FPTS(1:NFPTS,1)
    FPTS2 = this%FPTS(1:NFPTS,2)
    FPTS3 = this%FPTS(1:NFPTS,3)

    write (unit=unit, nml=SKYRMION, iostat=ios)

    deallocate(FPTS1)
    deallocate(FPTS2)
    deallocate(FPTS3)
  end subroutine sk_write_nml

  subroutine sk_load_nml(this, unit, ios)
    class(skyrmion_state), intent(OUT) :: this
    integer, intent(IN) :: unit
    integer, intent(OUT) :: ios

    double precision q, h, nu, mu, dummy
    integer selectedPredicate
    integer NFPTS
    double precision FPTS1(NFPTSMAX),FPTS2(NFPTSMAX), FPTS3(NFPTSMAX)

    namelist /SKYRMION/ q, h, nu, mu, selectedPredicate, &
         NFPTS, FPTS1, FPTS2, FPTS3

    selectedPredicate=0

    read (unit=unit, nml=SKYRMION, iostat=ios)

    dummy = state_evolve(this, q, h)
    this%mu = mu
    this%nu = nu
    this%selectedPredicate = selectedPredicate
    this%NFPTS = NFPTS
    this%FPTS(1:NFPTS,1) = FPTS1(1:NFPTS)
    this%FPTS(1:NFPTS,2) = FPTS2(1:NFPTS)
    this%FPTS(1:NFPTS,3) = FPTS3(1:NFPTS)
  end subroutine sk_load_nml

  subroutine sk_get_short_code(this, code)
    class(skyrmion_state), intent(IN) :: this
    character (len=2), intent(OUT) :: code

    code="Sk"
  end subroutine sk_get_short_code

  subroutine sk_print(this)
    implicit none
    class(skyrmion_state), intent(IN) :: this
    double precision q, h
    integer I
    call this%getxy(q,h)
    print '("q=",G0.5,"; h=",G0.5, "nu=", G0.5,"; mu=",G0.5)', &
         q, h, this%nu, this%mu
    print '("NFPTS=", I0)', this%NFPTS
    do I=1, this%NFPTS
       print '(G0.5, " ", G0.5)', this%FPTS(I,1), this%FPTS(I,2)
    end do
  end subroutine sk_print

end module skyrme_state
