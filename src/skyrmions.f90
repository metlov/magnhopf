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
!     This file contains functions for computing equilibrium skyrmion profiles.
module skyrmions
  use debug
  implicit none

  public :: PexSk, PzSk, PaSk, PdmSk
  public :: EtermsSk, EtotSk, SOLVE_EULER_SK

contains
  double precision function PexSk(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: FPTS(NFPTSMAX,3)

    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    PexSk=PROFILE_INTEGRATE(PexSk_INTEGRAND,0.0D0,1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PexSk_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PexSk_INTEGRAND = fp**2*pi**3*r + (pi*DSIN(f*pi)**2)/r
    end function PexSk_INTEGRAND
  end function PexSk

  double precision function PzSk(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: FPTS(NFPTSMAX,3)

    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    PzSk=PROFILE_INTEGRATE(PzSk_INTEGRAND,0.0D0,1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PzSk_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PzSk_INTEGRAND = 2.0D0*pi*(r - r*DCOS(f*pi))
    end function PzSk_INTEGRAND
  end function PzSk

  double precision function PaSk(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: FPTS(NFPTSMAX,3)

    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    PaSk=PROFILE_INTEGRATE(PaSk_INTEGRAND,0.0D0,1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PaSk_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PaSk_INTEGRAND = 2.0D0*pi*r*DSIN(f*pi)**2
    end function PaSk_INTEGRAND
  end function PaSk

  double precision function PdmSk(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: FPTS(NFPTSMAX,3)

    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    PdmSk=PROFILE_INTEGRATE(PdmSk_INTEGRAND,0.0D0,1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PdmSk_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PdmSk_INTEGRAND = pi*(-2.0D0*fp*pi*r - DSIN(2.0D0*f*pi))
    end function PdmSk_INTEGRAND
  end function PdmSk

  ! Evaluates the totals for all the energy terms of the skyrmions with a
  ! given profile, assuming they are packed into an HCP lattice
  !
  ! The output array TERMS must be defined as
  ! DOUBLE PRECISION TERMS(4)
  ! on exit it would hold the energy terms in the following order:
  ! TERMS(1)  --  Eex, total exchange energy per unit volume
  ! TERMS(2)  --  Edm, total Dzyaloshinskii-Moriys energy per unit volume
  ! TERMS(3)  --  Ez , total Zeeman energy per unit volume
  ! TERMS(4)  --  Ea , total anisotropy energy per unit volume
  subroutine EtermsSk (h, q, nu, NFPTSMAX, NFPTS, FPTS, TERMS)
    implicit none
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: h, q, nu, FPTS(NFPTSMAX,3)
    double precision, intent(OUT) :: TERMS(4)
    double precision V

    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    V = 2.0D0*DSQRT(3.0D0)

    TERMS(1) = nu**2 * PexSk(NFPTSMAX, NFPTS, FPTS)/V                      ! Eex
    TERMS(2) = nu * PdmSk(NFPTSMAX, NFPTS, FPTS)/V                         ! Edm
    TERMS(3) = 0.0D0
    if (h .ne.0.0D0) TERMS(3) = h * (PzSk(NFPTSMAX, NFPTS, FPTS) - V)/V    ! Ez
    TERMS(4) = 0.0D0
    if (q .ne.0.0D0) TERMS(4) = q/2.0D0*(PaSk(NFPTSMAX, NFPTS, FPTS)-V)/V  ! Ea
  end subroutine EtermsSk

  ! Evaluates the total energy of hopfions with the given profile, assuming
  ! the hopfions are packed into a 3D lattice with V R^3 volume per hopfion
  !
  ! For HCP and FCC lattices V = 4.0D0*DSQRT(2.0D0)
  double precision function EtotSk(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    implicit none
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(IN) :: h, q, nu, FPTS(NFPTSMAX,3)
    double precision TERMS(4)

    call EtermsSk(h, q, nu, NFPTSMAX, NFPTS, FPTS, TERMS)
    EtotSk=sum(TERMS)
  end function EtotSk

! Computes the equilibrium skyrmion profile, possibly using a specified "old"
! profile as an initial approximation.
!
! Parameters:
! IN:
!     h        - the external field
!     q        - the anisotropy quality factor
!     NFPTSMAX - dimension of FPTS(NPTSMAX,3) array (must be sufficient
!                to store the resulting profile points)
! IN/OUT:
!     nu     - the dimensionless skyrmion size parameter, takes initial
!              approximation on input, final equilibrium value on output
!     NFPTS  - ON INPUT: the number of points in the "old" profile to use
!              as initial aproximation or 0 if a built-in initial profile
!              must be used
!              ON OUTPUT: the number of points in the resulting equilibrium
!              profile or 0 if the equaton can not be solved.
!     FPTS   - ON INPUT: an "old" hopfion profile
!              ON OUTPUT: the resulting equilibrium hopfion profile.
!
subroutine SOLVE_EULER_SK(h, q, nu, NFPTSMAX, NFPTS, FPTS)
  implicit none
  integer, intent(IN) :: NFPTSMAX
  integer, intent(INOUT) ::  NFPTS
  double precision, intent(IN) :: h, q
  double precision, intent(INOUT) :: nu, FPTS(NFPTSMAX,3)

! copy of the initial hofion profile
  integer NFPTSOLD
  double precision FPTSOLD(NFPTSMAX,3)

! parameters to the MUSN call
  integer, parameter :: N = 4
  double precision A, B, ER(5), AMP
  integer NRTI
  integer, parameter :: ITLIM = 30
  double precision TI(NFPTSMAX), Y(N,NFPTSMAX), QQ(N,N,NFPTSMAX)
  integer, parameter :: NUDIM = 10 ! N*(N+1)/2
  double precision U(NUDIM,NFPTSMAX), PHIREC(NUDIM,NFPTSMAX)
  double precision D(N,NFPTSMAX)
  integer KPART, LW, LIW, LWG
  double precision W(7*N + 3*N*NFPTSMAX + 4*N*N  )
  integer IW(3*N + NFPTSMAX)
  double precision WGR(NFPTSMAX*6)
  integer IERROR

  double precision, parameter :: pi = 3.1415926535897932384626433832795d0

  integer I
  double precision left
  double precision NUSUM

  CALL ASSERT(NFPTS.le.NFPTSMAX)

  ! copy the old profile
  NFPTSOLD=NFPTS
  do I=1,NFPTSOLD
     FPTSOLD(I,1)=FPTS(I,1)
     FPTSOLD(I,2)=FPTS(I,2)
     FPTSOLD(I,3)=FPTS(I,3)
  end do

  left=0.001D0

  A      = left
  B      = 1.0D0
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
     CALL ASSERT(NRTI.lt.NFPTSMAX)
     FPTS(1,1)=0.0D0
     FPTS(1,2)=1.0D0
     FPTS(1,3)=(Y(2,I)-FPTS(1,2))/TI(1)    ! approximation to the derivative
     NUSUM=0.0D0
     do I=1,NRTI
        FPTS(1+I,1)=TI(I)
        FPTS(1+I,2)=Y(2,I)
        FPTS(1+I,3)=Y(1,I)
        NUSUM=NUSUM+Y(3,I)
     end do
     NFPTS=NRTI+1
     nu=NUSUM/NRTI
  else
     call xerror("MUSN RETURNED IERROR",IERROR,0)
     NFPTS=0
  end if
contains
  subroutine EQUATIONS(r,Y,FUN)
    implicit none
    double precision r, Y(4),FUN(4)
    double precision f, fp, nu, dEdnu
    f     = Y(2)
    fp    = Y(1)
    nu    = Y(3)
    dEdnu = Y(4)
    FUN(1) = (nu*r - fp*nu**2*pi*r - nu*r*DCOS(2.0D0*f*pi) + &
         h*r**2*DSIN(f*pi) + nu**2*DCOS(f*pi)*DSIN(f*pi) + &
         q*r**2*DCOS(f*pi)*DSIN(f*pi))/(nu**2*pi*r**2)
    FUN(2) = Y(1)
    FUN(3) = 0.0D0
    FUN(4) = (fp**2*nu*pi**3*r**2 + fp*pi**2*(-2.0D0 + fp*nu*pi)*r**2 + &
         2.0D0*nu*pi*DSIN(f*pi)**2 - pi*r*DSIN(2.0D0*f*pi))/ &
         (2.0D0*DSQRT(3.0D0)*r)
  end subroutine EQUATIONS
  subroutine BOUNDARY(N,YA,YB,FG,DGA,DGB)
    implicit none
    integer N
    double precision YA(N),YB(N),FG(N),DGA(N,N),DGB(N,N)

    CALL ASSERT(N.eq.4)
    FG(1) = YA(2) - (1.0D0-left)    ! f(A) = 1 - left
    FG(2) = YB(2)                   ! f(B) = 0
    FG(3) = YA(4)                   ! dEdnu(A) = 0
    FG(4) = YB(4)                   ! dEdnu(B) = 0

    DGA=0.0D0
    DGB=0.0D0

    DGA(1,2) = 1.D0
    DGB(2,2) = 1.D0
    DGA(3,4) = 1.D0
    DGB(4,4) = 1.D0
    return
  end subroutine BOUNDARY
  subroutine INITIAL_APPROX(T,Y)
    implicit none
    double precision, intent(IN) :: T
    double precision, intent(OUT) :: Y(4)
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    logical SKIP
    integer IERR
    if (NFPTS.GE.2) then
       SKIP=.false.
       call DPCHFD(NFPTS, FPTS(1,1), FPTS(1,2), FPTS(1,3), 1, SKIP, &
            1, T, Y(2), Y(1), IERR)
       CALL ASSERT(IERR.EQ.0)
    else
       Y(1) = -T
       Y(2) = 1.0D0 - T
    end if
    Y(3) = nu
    Y(4) = (pi*(2.0D0*nu*DSIN(pi*Y(2))**2 - T*DSIN(2.0D0*pi*Y(2)) + &
         2.0D0*pi*T**2*Y(1)*(-1.0D0 + nu*pi*Y(1))))/(2.0D0*DSQRT(3.0D0)*T)
  end subroutine INITIAL_APPROX
end subroutine SOLVE_EULER_SK

end module skyrmions
