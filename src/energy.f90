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
!     FUNCTIONS IN THIS FILE ARE USED TO COMPUTE THE HOPFION ENERGY
!     AND THE EULER'S EQUATION RESOLVED FOR THE SECOND DERIVATIVE
!
!     PARAMETERS, COMMON TO MOST FUNCTIONS:
!
!     HTYPE   AN INTEGER 1 OR 2, SPECIFYING THE HOPFION TYPE
!     G       THE HOPFION ASPECT RATIO (G<1 CORRESPONDS TO PROLATE
!             SPHEROIDS, G>1 TO OBLATE ONES, G=1 IS SPHERE)
!     NFPTS   AN INTEGER, SPECIFYING THE NUMBER OF CURRENT POINTS IN
!             THE FPTS ARRAY
!     FPTS    A DOUBLE ARRAY FPTS(NPTS,3) CONTAINING THE VALUES OF
!             THE HOPFION PROFILE ABSCISSAE FPTS(NPTS,1), THE CORRESPONDING
!             VALUES OF F IN FPTS(NPTS,2) AND THE VALUES OF ITS FIRST
!             DERIVATIVE IN FPTS(NPTS,3). THE ABSCISSAE ARRAY IS ASSUMED TO
!             BE SORTED IN ASCENDING ORDER
module energy
  use debug
  implicit none

  public :: Pex, Pz, Pa, Pdm, T1W1, Pms
  public :: Eterms, Etot, EulerFPP, EulerTaylor

contains
  !     EVALUATES THE EXCHANGE ENERGY
  double precision function Pex(HTYPE, g, NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer HTYPE, NFPTSMAX, NFPTS
    double precision FPTS(NFPTSMAX,3),g
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    Pex=PROFILE_INTEGRATE(PexT1_INTEGRAND, 0.0D0, 1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PexT1_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PexT1_INTEGRAND = &
           (64.0D0*pi*r**2*(fp**2*(4.0D0 + g**2)*r**2* &
           (1.0D0 + (-2.0D0 + f)*f + r**2)**2 - &
           2.0D0*(-1.0D0 + f)*fp*r*(1.0D0 + (-2.0D0 + f)*f + r**2)* &
           (5.0D0 + 5.0D0*(-2.0D0 + f)*f + (3.0D0 + 2.0D0*g**2)*r**2) + &
           (-1.0D0 + f)**2*(15.0D0 + 5.0D0*(-2.0D0 + f)*f* &
           (6.0D0 + 3.0D0*(-2.0D0 + f)*f + (2.0D0 + 4.0D0*g**2)*r**2) + &
           r**2*(10.0D0 + 11.0D0*r**2 + 4.0D0*g**2*(5.0D0 + r**2)))))/ &
           (15.0D0*(1.0D0 + (-2.0D0 + f)*f + r**2)**4)
    end function PexT1_INTEGRAND
  end function Pex

  !     EVALUATES THE ZEEMAN ENERGY (does not depend on g or hopfion type)
  double precision function Pz(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer NFPTSMAX, NFPTS
    double precision FPTS(NFPTSMAX,3)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    Pz=PROFILE_INTEGRATE(PzT1_INTEGRAND, 0.0D0, 1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PzT1_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision, intent(IN) :: r, f, fp
      PzT1_INTEGRAND = &
           (64.0D0*(-1.0D0 + f)**2*pi*r**4)/ &
           (3.0D0*(1.0D0 + (-2.0D0 + f)*f + r**2)**2)
    end function PzT1_INTEGRAND
  end function Pz

  !     EVALUATES THE ANISOTROPY ENERGY (does not depend on g or hopfion type)
  double precision function Pa(NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer NFPTSMAX, NFPTS
    double precision FPTS(NFPTSMAX,3)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    Pa=PROFILE_INTEGRATE(PaT1_INTEGRAND, 0.0D0, 1.0d0, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function PaT1_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision r, f, fp
      PaT1_INTEGRAND = &
           (128.0D0*(-1.0D0 + f)**2*pi*r**4* &
           (5.0D0 - 6.0D0*r**2 + 5.0D0*r**4 + &
           (-2.0D0 + f)*f*(10.0D0 + 5.0D0*(-2.0D0 + f)*f - 6.0D0*r**2)))/ &
           (15.0D0*(1.0D0 + (-2.0D0 + f)*f + r**2)**4)
    end function PaT1_INTEGRAND
  end function Pa

  !     EVALUATES THE DZYALOSHINSKII-MORIYA CHIRAL ENERGY
  double precision function Pdm(HTYPE, gg, NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer, intent(IN) :: HTYPE, NFPTSMAX, NFPTS
    double precision, intent(IN) :: FPTS(NFPTSMAX,3),gg
    double precision g
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    g = gg
    if (HTYPE.eq.2) g = -gg

    Pdm=PROFILE_INTEGRATE(PdmT1_INTEGRAND, 0.0D0, 1.0d0, NFPTSMAX, NFPTS, FPTS)

    if (HTYPE.eq.2) Pdm = -Pdm
  contains
    double precision function PdmT1_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision r, f, fp
      PdmT1_INTEGRAND = &
           (32.0D0*pi*r**2*(fp*r*(1.0D0 + (-2.0D0 + f)*f + r**2)* &
           (5.0D0 + (2.0D0 + 8.0D0*g)*r**2 + 5.0D0*r**4 + &
           (-2.0D0+f)*f*(10.0D0+5.0D0*(-2.0D0+f)*f+2.0D0*(1.0D0+4.0D0*g)*r**2))+&
           (-1.0D0 + f)*(-15.0D0 + 5.0D0*(3.0D0 - 8.0D0*g)*r**2 + &
           (-29.0D0 + 24.0D0*g)*r**4 + 5.0D0*r**6 + &
           (-2.0D0+f)*f*(-45.0D0+(30.0D0-80.0D0*g)*r**2+(-29.0D0+24.0D0*g)*r**4-&
           5.0D0*(-2.0D0 + f)*f*(9.0D0 + 3.0D0*(-2.0D0 + f)*f + &
           (-3.0D0 + 8.0D0*g)*r**2)))))/ &
           (15.0D0*(1.0D0 + (-2.0D0 + f)*f + r**2)**4)
    end function PdmT1_INTEGRAND
  end function Pdm

  !     EVALUATES THE W1 MAGNETOSTATIC FUNCTION FOR T1 HOPFIONS
  double precision function T1W1(XX, NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer NFPTSMAX, NFPTS
    double precision XX,FPTS(NFPTSMAX,3)
    double precision X
    CALL ASSERT(XX.ge.-0.01D0.and.XX.le.1.01D0)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    X=XX
    if (XX.lt.0.0D0) X=0.0D0
    if (XX.gt.1.0D0) X=1.0D0

    T1W1=PROFILE_INTEGRATE(T1W1_INTEGRAND, 0.0D0, X, NFPTSMAX, NFPTS, FPTS)
  contains
    double precision function T1W1_INTEGRAND(r, f, fp)
      double precision r, f, fp
      T1W1_INTEGRAND = &
           (16.0d0*(-1.0d0 + f)**2*r**4)/(1.0d0 + (-2.0d0 + f)*f + r**2)**2
    end function T1W1_INTEGRAND
  end function T1W1

  !     EVALUATES THE W1 MAGNETOSTATIC FUNCTION FOR T1 HOPFIONS
  double precision function T2Vi(I, XX, NFPTSMAX, NFPTS, FPTS)
    use debug
    use quadratures
    implicit none
    integer, intent(IN) :: I, NFPTSMAX, NFPTS
    double precision, intent(IN) :: XX,FPTS(NFPTSMAX,3)
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    double precision X
    CALL ASSERT(XX.ge.-0.01D0.and.XX.le.1.01D0)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    X=XX
    if (XX.lt.0.0D0) X=0.0D0
    if (XX.gt.1.0D0) X=1.0D0

    select case(I)
    case (1)
       T2Vi=PROFILE_INTEGRATE(T2V1_INTEGRAND, 0.0D0, X, NFPTSMAX, NFPTS, FPTS)
    case (2)
       T2Vi=PROFILE_INTEGRATE(T2V2_INTEGRAND, 0.0D0, X, NFPTSMAX, NFPTS, FPTS)
    case (3)
       T2Vi=PROFILE_INTEGRATE(T2V3_INTEGRAND, 0.0D0, X, NFPTSMAX, NFPTS, FPTS)
    case (4)    ! V4 = 2*V2 - V3 , but evaluated in a single integration sweep
       T2Vi=PROFILE_INTEGRATE(T2V4_INTEGRAND, 0.0D0, X, NFPTSMAX, NFPTS, FPTS)
    case default
       CALL ASSERT(.false.)
    end select
  contains
    double precision function T2V1_INTEGRAND(r, f, fp)
      double precision r, f, fp
      T2V1_INTEGRAND = &
           (256.0D0*(-1.0D0 + f)*pi*r**7* &
           (-2.0D0*(-1.0D0 + f)*r + fp*(1.0D0 - f + r)*(-1.0D0 + f + r)))/ &
           (1225.0D0*((-1.0D0 + f)**2 + r**2)**3)
    end function T2V1_INTEGRAND
    double precision function T2V2_INTEGRAND(r, f, fp)
      double precision r, f, fp
      T2V2_INTEGRAND = &
           (32.0D0*(-1.0D0 + f)*pi*r**4* &
           (-5.0D0*(-1.0D0 + f)**3 + 4.0D0*(-1.0D0 + f)**2*fp*r + &
           3.0D0*(-1.0D0 + f)*r**2 - 4.0D0*fp*r**3))/ &
           (75.0D0*((-1.0D0 + f)**2 + r**2)**3)
    end function T2V2_INTEGRAND
    double precision function T2V3_INTEGRAND(r, f, fp)
      double precision r, f, fp
      T2V3_INTEGRAND = &
           (32.0D0*(-1.0D0 + f)*pi*r**4* &
           (-5.0D0*(-1.0D0 + f)**3 + 4.0D0*(-1.0D0 + f)**2*fp*r + &
           3.0D0*(-1.0D0 + f)*r**2 - 4.0D0*fp*r**3))/ &
           (45.0D0*((-1.0D0 + f)**2 + r**2)**3)
    end function T2V3_INTEGRAND
    double precision function T2V4_INTEGRAND(r, f, fp)
      double precision r, f, fp
      T2V4_INTEGRAND = &
           (32.0D0*(-1.0D0 + f)*pi*r**4* &
           (-5.0D0*(-1.0D0 + f)**3 + 4.0D0*(-1.0D0 + f)**2*fp*r + &
           3.0D0*(-1.0D0 + f)*r**2 - 4.0D0*fp*r**3))/ &
           (225.0D0*((-1.0D0 + f)**2 + r**2)**3)
    end function T2V4_INTEGRAND
  end function T2Vi

  !     EVALUATES THE MAGNETOSTATIC ENERGY
  double precision function Pms(HTYPE, g, NFPTSMAX, NFPTS, FPTS)
    use quadratures
    implicit none
    integer HTYPE, NFPTSMAX, NFPTS
    double precision FPTS(NFPTSMAX,3),g
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(g.eq.1.0d0)     ! TODO: currently magnetostatics is only
    !                             computed for spherical hopfions
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    select case(HTYPE)
    case (1)
       Pms=2.0D0*PROFILE_INTEGRATE(PmsT1_INTEGRAND, 0.0D0, 1.0d0, &
            NFPTSMAX, NFPTS, FPTS)
    case (2)
       Pms=2.0D0*PROFILE_INTEGRATE(PmsT2_INTEGRAND, 0.0D0, 1.0d0, &
            NFPTSMAX, NFPTS, FPTS)
    case default
       CALL ASSERT(.false.)
    end select
    ! The factor 2 above is there because we actually are integrating only
    ! over the triangle 0<r1<1, 0<r2<r1
  contains
    double precision function PmsT1_INTEGRAND(r, f, fp)
      double precision, parameter :: pi = 3.1415926535897932384626433832795d0
      double precision r, f, fp
      PmsT1_INTEGRAND = &
           (32.0D0*(-1.0D0 + f)**2*pi*r*T1W1(r,NFPTSMAX,NFPTS,FPTS)) / &
           (9.0D0*(1.0D0 + (-2.0D0 + f)*f + r**2)**2)
    end function PmsT1_INTEGRAND
    double precision function PmsT2_INTEGRAND(r, f, fp)
      double precision r, f, fp
      double precision V1, V2, V3
      V1 = T2Vi(1, r, NFPTSMAX, NFPTS, FPTS)
      V2 = T2Vi(2, r, NFPTSMAX, NFPTS, FPTS)
      V3 = T2Vi(3, r, NFPTSMAX, NFPTS, FPTS)
      PmsT2_INTEGRAND = &
           (16.0D0*(-1.0D0 + f)* &
           (2.0D0*fp*(1.0D0 - f + r)*(-1.0D0 + f + r)* &
           (V1 + r**2*(V2 - V3)) + &
           (-1.0D0 + f)*r*(-4.0D0*V1 - (-1.0D0 + f)**2*V3 + &
           r**2*(-4.0D0*V2 + 3.0D0*V3))))/ &
           ((-1.0D0 + f)**2 + r**2)**3
    end function PmsT2_INTEGRAND
  end function Pms

  ! Evaluates the totals for all the energy terms of the hopfion with a
  ! given profile, assuming that the hopfions are packed into a 3D lattice
  ! with V R^3 volume per hopfion
  !
  ! For HCP and FCC lattices V = 4.0D0*DSQRT(2.0D0)
  !
  ! The output array TERMS must be defined as
  ! DOUBLE PRECISION TERMS(5)
  ! on exit it would hold the energy terms in the following order:
  ! TERMS(1)  --  Eex, total exchange energy per unit volume
  ! TERMS(2)  --  Edm, total Dzyaloshinskii-Moriys energy per unit volume
  ! TERMS(3)  --  Ez , total Zeeman energy per unit volume
  ! TERMS(4)  --  Ea , total anisotropy energy per unit volume
  ! TERMS(5)  --  Ems, total magnetostatic energy per unit volume
  subroutine Eterms (HTYPE, V, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS, TERMS)
    implicit none
    integer, intent(IN) :: HTYPE, NFPTSMAX, NFPTS
    double precision, intent(IN) :: V, g, h, q, nu, mu, FPTS(NFPTSMAX,3)
    double precision, intent(OUT) :: TERMS(5)
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(g.gt.0.0D0)
    CALL ASSERT(NFPTS.le.NFPTSMAX)
    CALL ASSERT(NFPTS.ge.2)

    TERMS(1) = nu**2 * Pex(HTYPE, g, NFPTSMAX, NFPTS, FPTS)/V             ! Eex
    TERMS(2) = nu * Pdm(HTYPE, g, NFPTSMAX, NFPTS, FPTS)/V                ! Edm
    TERMS(3) = 0
    if (h .ne.0.0D0) TERMS(3) = h * (Pz(NFPTSMAX, NFPTS, FPTS) - V)/V     ! Ez
    TERMS(4) = 0
    if (q .ne.0.0D0) TERMS(4) = q/2.0D0*(Pa(NFPTSMAX, NFPTS, FPTS)-V)/V   ! Ea
    TERMS(5) = 0
    if (mu.gt.0.0D0) TERMS(5) = mu**2*Pms(HTYPE,g,NFPTSMAX,NFPTS,FPTS)/V  ! Ems
  end subroutine Eterms

  ! Evaluates the total energy of hopfions with the given profile, assuming
  ! the hopfions are packed into a 3D lattice with V R^3 volume per hopfion
  !
  ! For HCP and FCC lattices V = 4.0D0*DSQRT(2.0D0)
  double precision function Etot(HTYPE, V, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    implicit none
    integer, intent(IN) :: HTYPE, NFPTSMAX, NFPTS
    double precision, intent(IN) :: V, g, h, q, nu, mu, FPTS(NFPTSMAX,3)
    double precision TERMS(5)

    call Eterms(HTYPE, V, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS, TERMS)
    Etot=sum(TERMS)
  end function Etot

  !     EVALUATES THE R.H.S. OF THE EULER'S EQUATION, RESOLVED FOR
  !     THE SECOND DERIVATIVE F''(r)
  double precision function EulerFPP(HTYPE, g, r, f, fp, nu, h, &
       q, mu, NFPTSMAX, NFPTS, FPTS)
    implicit none
    integer, intent(IN) :: HTYPE, NFPTSMAX, NFPTS
    double precision, intent(IN) :: g, r, f, fp, nu, h, q, mu
    double precision, intent(IN) :: FPTS(NFPTSMAX,3)
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    double precision W1, V4
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(NFPTS.eq.0.or.NFPTS.ge.2)
    CALL ASSERT(.not.(mu.gt.0.0D0.and.NFPTS.lt.2))
    CALL ASSERT(g.eq.1.0d0.or.NFPTS.eq.0) ! TODO: currently magnetostatics is
    !                                       only computed for spherical hopfions
    CALL ASSERT(NFPTS.le.NFPTSMAX)


    !if (f.LT.0.0D0) write(0,'("f(",G0.5,")=", G0.5)') r,f
    !if (DABS(f).GT.10.0D0) write(0,'("f(",G0.5,")=", G0.5)') r,f
    if (DABS(f).GT.10.0D0) then
       ! if the solution starts to break down, let it break down to the end
       EulerFPP=10.0D0
       return
    end if
    select case (HTYPE)
    case (1)
       W1=0.0D0
       if (mu.gt.0.0D0.and.NFPTS.ge.2) then
          W1=T1W1(r, NFPTSMAX, NFPTS, FPTS)
       end if
       EulerFPP = &
               -0.3333333333333333333333333D0* &
       (-15.0D0*f**6*(4.0D0*(-2.0D0 + g)*nu + 7.0D0*(h + q))*r**3 + &
          12.0D0*(-1.0D0 + f)**2*fp*(4.0D0 + g**2)*nu**2*r**2* &
           (1.0D0 + (-2.0D0 + f)*f + r**2)**2 - &
          6.0D0*(-1.0D0 + f)*fp**2*(4.0D0 + g**2)*nu**2*r**3* &
           (1.0D0 + (-2.0D0 + f)*f + r**2)**2 + &
          3.0D0*f**5*r**3*(-40.0D0*nu*(6.0D0+nu) + 30.0D0*g*nu*(4.0D0+g*nu) + &
             3.0D0*q*(35.0D0 - 9.0D0*r**2) + 5.0D0*h*(21.0D0 + r**2)) - &
          3.0D0*f**4*r**3*(5.0D0*q*(35.0D0-27.0D0*r**2)+25.0D0*h*(7.0D0+r**2)+ &
             2.0D0*nu*(-100.0D0*(3.0D0 + nu) + 75.0D0*g*(2.0D0 + g*nu) + &
                (24.0D0 - 44.0D0*g)*r**2)) - &
          3.0D0*f**3*r**3*(4.0D0*nu* &
              (-25.0D0*(-4.0D0*(2.0D0 + nu) + g*(4.0D0 + 3.0D0*g*nu)) + &
                (-4.0D0*(12.0D0 + nu) + g*(88.0D0 + 9.0D0*g*nu))*r**2) + &
             q*(-175.0D0 + 270.0D0*r**2 - 27.0D0*r**4) + &
             5.0D0*h*(-35.0D0 - 10.0D0*r**2 + r**4)) + &
          3.0D0*f**2*r**3*(5.0D0*h*(-21.0D0 - 10.0D0*r**2 + 3.0D0*r**4) - &
             3.0D0*q*(35.0D0 - 90.0D0*r**2 + 27.0D0*r**4) + &
             4.0D0*nu*(-25.0D0*(-6.0D0 - 4.0D0*nu + 3.0D0*g*(1.0D0 + g*nu)) + &
                3.0D0*(-4.0D0*(6.0D0 + nu) + g*(44.0D0 + 9.0D0*g*nu))*r**2 - &
                5.0D0*(-2.0D0 + g)*r**4)) - &
          3.0D0*f*r**3*(5.0D0*h*(1.0D0 + r**2)*(-7.0D0 + 2.0D0*r**2 + r**4) + &
             2.0D0*nu*(20.0D0*(6.0D0+5.0D0*nu) - 15.0D0*g*(4.0D0+5.0D0*g*nu) + &
                2.0D0*(-12.0D0*(4.0D0 + nu) + g*(88.0D0 + 27.0D0*g*nu))*r**2 + &
                (40.0D0 - 20.0D0*g + 4.0D0*nu + g**2*nu)*r**4) + &
             q*(-35.0D0 + 135.0D0*r**2 - 81.0D0*r**4 + 5.0D0*r**6)) + &
          3.0D0*r**3*(5.0D0*h*(-1.0D0 + r**2)*(1.0D0 + r**2)**2 + &
             2.0D0*nu*(20.0D0*(1.0D0 + nu) - 5.0D0*g*(2.0D0 + 3.0D0*g*nu) + &
                2.0D0*(-4.0D0*(3.0D0 + nu) + g*(22.0D0 + 9.0D0*g*nu))*r**2 + &
                (4.0D0*(5.0D0 + nu) + g*(-10.0D0 + g*nu))*r**4) + &
             q*(-5.0D0 + 27.0D0*r**2 - 27.0D0*r**4 + 5.0D0*r**6)) - &
          5.0D0*mu**2*W1 - 35.0D0*f**6*mu**2*W1 - &
          25.0D0*f**4*mu**2*(7.0D0 + r**2)*W1 + &
          5.0D0*f**5*mu**2*(21.0D0 + r**2)*W1 - &
          5.0D0*f**3*mu**2*(-35.0D0 - 10.0D0*r**2 + r**4)*W1 + &
          5.0D0*mu**2*r**2*(-1.0D0 + r**2 + r**4)*W1 - &
          5.0D0*f*mu**2*(1.0D0 + r**2)*(-7.0D0 + 2*r**2 + r**4)*W1 + &
          5.0D0*f**2*mu**2*(-21.0D0 - 10.0D0*r**2 + 3.0D0*r**4)*W1 + &
          5.0D0*f**7*(3.0D0*(h + q)*r**3 + mu**2*W1))/ &
          ((4.0D0 + g**2)*nu**2*(r + (-2.0D0 + f)*f*r + r**3)**3)
    case (2)
       V4=0.0D0
       if (mu.gt.0.0D0) then
          V4=T2Vi(4, r, NFPTSMAX, NFPTS, FPTS)
       end if
       EulerFPP = &
            -0.5D0*(2*pi*r**4*(-4*nu*(5*(2 + g) - 2*(6 + 11*g)*r**2 + &
            5*(2 + g)*r**4) + &
            2*nu**2*(4*(5 - 2*r**2 + r**4) + g**2*(-15 + 18*r**2 + r**4)) + &
            (-1 + r**2)*(5*h*(1 + r**2)**2 + q*(5 - 22*r**2 + 5*r**4))) + &
            (r*(2940*fp**2*(4 + g**2)*nu**2*pi*(r + r**3)**3 + &
            8*fp*pi*r**2*(1856*mu**2*r**4*(-1 + r**2)**2 + &
            735*(4 + g**2)*nu**2*(1 + r**2)**3) + &
            3675*f**9*(2*pi*(h + q)*r**3 - 3*mu**2*V4) + &
            735*f**8*(2*pi*r**2*(4*fp*(4 + g**2)*nu**2 - &
            5*(4*(2 + g)*nu + 9*(h + q))*r) + &
            135*mu**2*V4) + mu**2*(-1 + r**2)* &
            (64*pi*r**5*(-245 + 219*r**2) - 11025*(1 + r**2)**3*V4) + &
            2*f**5*(2*pi*r**2*(3*fp*(-7424*mu**2*r**4 - &
            735*fp*(4 + g**2)*nu**2*r*(7 + r**2) - &
            980*(4 + g**2)*nu**2*(28 + 9*r**2)) + &
            r*(77175*h*(3 + r**2) - 15435*q*(-15 + 11*r**2) + &
            16*mu**2*r**2*(-5145 + 464*r**2) - &
            2940*nu*(-140*(2 + g) + (6 + 51*g)*r**2) - &
            2205*nu**2*(35*(4 - 3*g**2) + (4 + g**2)*r**2))) - &
            231525*mu**2*(3 + r**2)*V4) + &
            2*f**6*(2*pi*r**2*(245*r*(-840*(2 + g)*nu + &
            105*(4 - 3*g**2)*nu**2 - 630*q + &
            6*(2 + 17*g)*nu*r**2 + 7*(16*mu**2 + 33*q)*r**2 - &
            105*h*(6 + r**2)) + &
            fp*(5145*fp*(4 + g**2)*nu**2*r + 3712*mu**2*r**4 + &
            1470*(4 + g**2)*nu**2*(28 + 3*r**2))) + &
            77175*mu**2*(6 + r**2)*V4) + &
            490*f**7*(2*pi*r**2*(-3*fp*(4 + g**2)*nu**2*(16 + fp*r) + &
            r*(15*nu*(32 + 16*g - 4*nu + 3*g**2*nu) + 270*q - &
            16*mu**2*r**2 - 33*q*r**2 + &
            15*h*(18 + r**2))) - 45*mu**2*(18 + r**2)*V4) + &
            2*f**4*(2*pi*r**2*(fp*(3712*mu**2*r**4*(15 - 2*r**2) + &
            3675*fp*(4 + g**2)*nu**2*r*(7 + 3*r**2) + &
            1470*(4 + g**2)*nu**2*(70 + 45*r**2 + 3*r**4)) + &
            5*r*(16*mu**2*r**2*(1715 - 464*r**2) - 5145*h*(9 + 5*r**2) + &
            5145*q*(-9 + 11*r**2) + 735*nu**2*(35*(4 - 3*g**2) + &
            3*(4 + g**2)*r**2) + &
            294*nu*(-350*(2 + g) + (2 + 17*g)*r**2*(15 + r**2)))) + &
            77175*mu**2*(9 + 5*r**2)*V4) + &
            2*f**3*(2*pi*r**2*(fp*(14848*mu**2*r**4*(-5 + 2*r**2) - &
            735*fp*(4 + g**2)*nu**2*r*(35 + 30*r**2 + 3*r**4) - &
            5880*(4 + g**2)*nu**2*(14 + 3*r**2*(5 + r**2))) + &
            r*(-16*mu**2*r**2*(8575 - 4640*r**2 + 219*r**4) - &
            735*nu**2*(700 - 525*g**2 + 30*(4 + g**2)*r**2 + &
            (-4 + 19*g**2)*r**4) - &
            3675*h*(-42 - 35*r**2 + r**6) - &
            5880*nu*(-70*(2 + g) + (2 + 17*g)*r**2*(5 + r**2)) + &
            735*q*(210 + 11*r**2*(-35 + r**4)))) + &
            11025*mu**2*(-42 - 35*r**2 + r**6)*V4) + &
            2*f**2*(2*fp*pi*r**2*(2205*fp*(4+g**2)*nu**2*r*(1+r**2)*(7+3*r**2)+&
            3712*mu**2*r**4*(15 - 12*r**2 + r**4) + &
            1470*(4 + g**2)*nu**2*(1 + r**2)*(28 + 17*r**2 + r**4)) + &
            pi*r**3*(3675*h*(1 + r**2)*(-35 - 6*r**2 + 5*r**4) + &
            32*mu**2*r**2*(5145 - 4640*r**2 + 657*r**4) - &
            735*q*(175 - 435*r**2 - 27*r**4 + 71*r**6) + &
            5880*nu**2*(-4*(-25 - 8*r**2 + r**4) + &
            g**2*(-75 + 3*r**2 + 14*r**4)) - &
            2940*nu*(135*(2 + g) + r**2*(-18 - 233*g + &
            r**2*(-22 - 107*g + 5*(2 + g)*r**2)))) - &
            33075*mu**2*(-6 - 7*r**2 + r**6)*V4) + &
            f*(2*pi*r**2*(-32*mu**2*r**3*(1715 - 2320*r**2 + 657*r**4) + &
            5880*nu*(r + r**3)*(15*(2 + g) - 4*(6 + 11*g)*r**2 + &
            5*(2 + g)*r**4) - &
            735*(r + r**3)*(-35*(h + q) + 5*(-5*h + 27*q)*r**2 + &
            3*(5*h - 27*q)*r**4 + &
            5*(h + q)*r**6) + &
            2*fp*(-2940*(4 + g**2)*nu**2*(1 + r**2)**2*(4 + r**2) - &
            735*fp*(4 + g**2)*nu**2*r*(1 + r**2)**2*(7 + r**2) - &
            7424*mu**2*r**4*(3 - 4*r**2 + r**4)) - &
            1470*nu**2*(r + r**3)*(4*(25 - 6*r**2 + r**4) + &
            g**2*(-75 + 54*r**2 + r**4))) + &
            11025*mu**2*(-9 - 14*r**2 + 6*r**6 + r**8)*V4)))/&
            (735.*(1 + (-2 + f)*f + r**2)))/&
            ((4 + g**2)*nu**2*pi*r**4*(1 + (-2 + f)*f + r**2)**3)
    case default
       CALL ASSERT(.false.)
    end select
  end function EulerFPP

  ! Computes the analytical Euler's equation solution at r<<1, based
  ! on Taylor expansion. This is used to start integration from some small
  ! 0<r0<<1 avoiding potential indeterminacies in the explicit formula for
  ! f''[r] in the equation's r.h.s. at r=0.
  !
  ! On output fl=f[r] and fpl=f'[r], other parameters are for input.
  subroutine EulerTaylor(HTYPE, g, r, nu, h, &
       q, mu, NFPTSMAX, NFPTSOLD, FPTSOLD, fl, fpl)
    implicit none
    integer, intent(IN) :: HTYPE, NFPTSMAX, NFPTSOLD
    double precision g, r, nu, h, q, mu, FPTSOLD(NFPTSMAX,3), fl, fpl

    double precision famp
    CALL ASSERT(HTYPE.ge.1.and.HTYPE.le.2)
    CALL ASSERT(g.gt.0.0D0)
    CALL ASSERT(.not.(mu.gt.0.0D0.and.NFPTSOLD.lt.2))
    CALL ASSERT(NFPTSOLD.le.NFPTSMAX)

    famp = (h-8.0D0*nu + 4.0D0*g*nu - 8.0D0*nu**2 + 6.0D0*g**2*nu**2 + q)/ &
         ((4.0D0 + g**2)*nu**2)
    if (HTYPE.EQ.2) then
       famp = famp + 16.0D0/(4.0D0*nu + g**2*nu)
    end if
    fl=famp*r**2/2.0D0
    fpl=famp*r
  end subroutine EulerTaylor
end module energy
