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

module quadratures
  use debug
  implicit none
  public :: PROFILE_INTEGRATE

contains
  ! PROFILE_INTEGRATE
  !
  ! INTEGRATES A FUNCTION OVER A HOPFION PROFILE, SPECIFIED NUMERICALLY
  !
  ! PARAMETERS:
  ! INTEGRAND -- Function to integrate. Must have signature
  !              DOUBLE PRECISION FI(r, f, fp, g)
  !              WHERE r is the current radial coordinate, f is the value
  !              of the profile function at r (f=f[r]) and fp is the value
  !              of the profile's function derivative at r (fp=f'[r])
  ! G         -- The hopfion's aspect ratio
  ! A         -- Lower limit of the integration
  ! B         -- Upper limit of the integration
  ! NFPTS     -- Number of the profile function points in FPTS array
  ! FPTS      -- The array (NFPTS,3) containing the numerical profile function
  !              abscissae, values and the first derivatives at these abscissae.
  double precision recursive function PROFILE_INTEGRATE(INTEGRAND, A, B, &
       NFPTSMAX, NFPTS, FPTS)
    implicit none
    external INTEGRAND
    integer NFPTSMAX, NFPTS
    double precision INTEGRAND, A, B, FPTS(NFPTSMAX,3)

    double precision epsabs, epsrel, result, abserr
    integer I, neval, ier
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(A.ge.FPTS(1,1))
    CALL ASSERT(B.le.FPTS(NFPTS,1))
    CALL ASSERT(B.ge.A)

    epsabs=1.0d-7
    epsrel=1.0d-7
    call VDQNG(PROFILE_INTEGRAND, a, b, epsabs, epsrel, &
         result, abserr, neval, ier)
    PROFILE_INTEGRATE=result
  contains
    recursive subroutine PROFILE_INTEGRAND(N, X, F)
      integer N
      double precision X(N),F(N)
      double precision FVALS(N), FPVALS(N)
      logical SKIP
      integer IERR
      CALL ASSERT(N.gt.0)
      SKIP=.false.
      call DPCHFD(NFPTS, FPTS(1,1), FPTS(1,2), FPTS(1,3), 1, SKIP, &
           N, X, FVALS, FPVALS, IERR)
      CALL ASSERT(IERR.EQ.0)
      do I=1,N
         F(I)=INTEGRAND(X(I), FVALS(I), FPVALS(I))
      end do
    end subroutine PROFILE_INTEGRAND
  end function PROFILE_INTEGRATE
end module quadratures
