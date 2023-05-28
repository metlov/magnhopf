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
!     THIS FILE DEFINES THE ENTRY POINT TO THE TESTING SUITE.
#define UNIT_TEST(t) IF(.NOT.(t)) CALL TF(__LINE__); CALL TP()
#define TOL_EQ(A, B) ABS((A)-(B))>0.000001
#define UNIT_TEST_EQ(val, exp) IF(TOL_EQ(val,exp)) CALL TFVE(__LINE__, val, exp); CALL TP()
#define TOL_EQ_WEAK(A, B) ABS((A)-(B))>0.002
#define UNIT_TEST_EQ_WEAK(val, exp) IF(TOL_EQ_WEAK(val,exp)) CALL TFVE(__LINE__, val, exp); CALL TP()

subroutine TESTINTEGRAND(N, X, F)
  use debug
  implicit none
  integer I,N
  double precision X(*), F(*)
  double precision PREV
  do I=1,N
     ! check the ordering of the supplied abscissae
     if (I>2) then
        CALL ASSERT(X(I)>PREV)
     end if
     PREV=X(I)
     ! evaluate the function
     F(I)=X(I)**2
  end do
end subroutine TESTINTEGRAND


!     Test driver program
program TEST
  use debug
  call TEST_INITIALIZE

  call TEST_VDQNG
  call TEST_T1W1
  call TEST_T2Vi
  call TEST_PFUNCS
  call TEST_EULERFPP
  call TEST_TAYLOR
  call TEST_SOLVE_EULER
  call TEST_ITERATE_MS
  call TEST_TEST_STATE
  call TEST_TRACER
  call TEST_HOPFION_STATE
  call TEST_SKYRMIONS
  call TEST_GROUNDSTATE

  IF (TEST_REPORT().GT.0) STOP 1
contains

  subroutine TEST_INITIALIZE
    integer TESTS_TOTAL, TESTS_FAILED
    common /utest/ TESTS_TOTAL, TESTS_FAILED
    TESTS_TOTAL = 0
    TESTS_FAILED = 0
  end subroutine TEST_INITIALIZE

  integer function TEST_REPORT ()
    integer TESTS_TOTAL, TESTS_FAILED
    common /utest/ TESTS_TOTAL, TESTS_FAILED
    write (*,*)"Test run complete"
    write (*,*)"Tests passed:", TESTS_TOTAL - TESTS_FAILED
    write (*,*)"Tests failed:", TESTS_FAILED
    TEST_REPORT = TESTS_FAILED
  end function TEST_REPORT

  subroutine TP()
    integer TESTS_TOTAL, TESTS_FAILED
    common /utest/ TESTS_TOTAL, TESTS_FAILED
    TESTS_TOTAL = TESTS_TOTAL + 1
  end subroutine TP

  subroutine TF(line)
    integer,intent(IN):: line
    integer TESTS_TOTAL, TESTS_FAILED
    common /utest/ TESTS_TOTAL, TESTS_FAILED
    write (*,'(A,I0,A)') "TEST FAILED on line ", LINE, " ."
    TESTS_FAILED = TESTS_FAILED + 1
  end subroutine TF

  subroutine TFVE(LINE, VAL, EXPECTED)
    integer TESTS_TOTAL, TESTS_FAILED
    common /utest/ TESTS_TOTAL, TESTS_FAILED
    integer line
    double precision VAL, EXPECTED
    write (*,'(A,I0,A)',ADVANCE="NO") &
         "TEST FAILED on line ", LINE, ": expected "
    if (EXPECTED.EQ.DBLE(INT(EXPECTED))) then
       if (VAL.EQ.DBLE(INT(VAL))) then
          write (*,'(I0,A,I0,A)') INT(EXPECTED), " got ", INT(VAL), " ."
       else
          write (*,'(I0,A,G0.10,A)') INT(EXPECTED), " got ", VAL, " ."
       end if
    else
       write (*,'(G0.10,A,G0.10,A)') EXPECTED, " got ", VAL
    end if
    TESTS_FAILED = TESTS_FAILED + 1
  end subroutine TFVE

  subroutine TEST_VDQNG
    implicit none
    external TESTINTEGRAND
    double precision a,b, epsabs, epsrel, result, abserr
    integer neval, ier
    a=0.0d0
    b=1.0d0
    epsabs=1.0d-7
    epsrel=1.0d-7
    call VDQNG(TESTINTEGRAND, a, b, epsabs, epsrel, result, abserr, neval, ier)
    !  print *, "neval=", neval, " abserr=", abserr
    UNIT_TEST_EQ(result,0.33333333333D0)
  end subroutine TEST_VDQNG


  ! fill FPTS array with f(r)=r^2
  subroutine FRSQUARE(NFPTS,FPTS)
    implicit none
    integer NFPTS
    double precision FPTS(NFPTS,3)
    integer I
    do I=1,NFPTS
       FPTS(I,1)=(I-1.0D0)/(NFPTS-1.0D0)
       FPTS(I,2)=FPTS(I,1)**2
       FPTS(I,3)=2*FPTS(I,1)
    end do
  end subroutine FRSQUARE

  subroutine TEST_T1W1
    use energy
    implicit none
    integer,parameter :: NFPTS=20
    double precision FPTS(NFPTS,3)
    double precision res

    call FRSQUARE(NFPTS,FPTS)

    res=T1W1(0.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    res=T1W1(0.321D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.010767384837782204867D0)

    res=T1W1(0.673D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.32339958783908859908D0)

    res=T1W1(1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.65660542509739933908D0)

  end subroutine TEST_T1W1

  subroutine TEST_T2Vi
    use energy
    implicit none
    integer,parameter :: NFPTS=20
    double precision FPTS(NFPTS,3)
    double precision res

    call FRSQUARE(NFPTS,FPTS)

    res=T2Vi(1, 0.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    res=T2Vi(1, 0.321D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-9.9627914118210910899D-7)

    res=T2Vi(1, 0.673D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.0033740177670486638820D0)

    res=T2Vi(1, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.043398550824031225570D0)

    res=T2Vi(2, 0.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    res=T2Vi(2, 0.321D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.0044140341476204913103D0)

    res=T2Vi(2, 0.673D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.060304872217170728769D0)

    res=T2Vi(2, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.27503823730575906845D0)

    res=T2Vi(3, 0.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    res=T2Vi(3, 0.321D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.0073567235793674855171D0)

    res=T2Vi(3, 0.673D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.10050812036195121462D0)

    res=T2Vi(3, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.45839706217626511408D0)

    res=T2Vi(4, 0.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    res=T2Vi(4, 0.321D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.0014713447158734971034D0)

    res=T2Vi(4, 0.673D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-0.020101624072390242923D0)

    res=T2Vi(4, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.091679412435253022817D0)

  end subroutine TEST_T2VI

  subroutine TEST_PFUNCS
    use energy
    implicit none
    integer,parameter :: NFPTS=20
    double precision FPTS(NFPTS,3)
    double precision res

    call FRSQUARE(NFPTS,FPTS)

    res=Pex(1, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,100.87076474568137151D0)

    res=Pex(1, 21.0D0/37.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,84.111004423340882119D0)

    res=Pex(1, 81.0D0/29.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,269.03100086768466128D0)

    res=Pz(NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,2.7503823730575906842D0)

    res=Pa(NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,2.1935798194235852631D0)

    res=Pdm(1, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,20.538738350437707216D0)

    res=Pdm(1, 21.0D0/37.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,18.867221882960053221D0)

    res=Pdm(1, 81.0D0/29.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,27.469767840582117310D0)

    res=Pdm(2, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-12.807974688353557493D0)

    res=Pdm(2, 21.0D0/37.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-14.479491155831211487D0)

    res=Pdm(2, 81.0D0/29.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,-5.8769451982091473980D0)

    res=Pms(1, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.85410163304684600048D0)

    res=Pms(2, 1.0D0, NFPTS, NFPTS, FPTS)
    UNIT_TEST_EQ(res,0.77408965244537254445D0)

  end subroutine TEST_PFUNCS

  subroutine TEST_EULERFPP
    use energy
    implicit none
    integer,parameter :: NFPTS=20
    integer TNFPTS
    double precision FPTS(NFPTS,3)
    double precision g, r, f, fp, nu, h, q, mu, res

    call FRSQUARE(NFPTS,FPTS)

    ! T1 hopfions
    TNFPTS=NFPTS

    g=1.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r

    nu=0.324D0
    h=0.4378D0
    q=0.1254D0
    mu=1.0D0/0.13D0
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,6.295589448652320879D0)

    r=0.001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,-16.980053345276892435D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,-16.980643194634963592D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    g=21.0D0/37.0D0
    TNFPTS = 0
    mu=0.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,-27.044550907220708196D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,-29.270189715027807913D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    g=81.0D0/29.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,9.0421857715993366485D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,14.863759709788049781D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(1, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)


    ! T2 hopfions
    TNFPTS=NFPTS

    g=1.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r

    nu=0.324D0
    h=0.4378D0
    q=0.1254D0
    mu=1.0D0/0.13D0
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,-19.429196669140286372D0)

    r=0.001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,32.400771813838703787D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,32.402072854747750566D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    g=21.0D0/37.0D0
    TNFPTS = 0
    mu=0.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,19.411330798093772062D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,27.857525572791135798D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)

    g=81.0D0/29.0D0

    r=0.214
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,26.056102725403032425D0)


    r=0.000000001D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,35.786109431466066922D0)

    r=1.0D0
    f=r**2
    fp=2.0D0*r
    res=EulerFPP(2, g, r, f, fp, nu, h, q, mu, NFPTS, TNFPTS, FPTS)
    UNIT_TEST_EQ(res,0.0D0)
  end subroutine TEST_EULERFPP

  subroutine TEST_TAYLOR
    use energy
    implicit none
    integer,parameter :: NFPTS=20
    double precision FPTS(NFPTS,3)
    double precision g, r, nu, h, q, mu, fl, fpl

    call FRSQUARE(NFPTS,FPTS)

    ! T1 hopfions
    g=1.0D0
    r=0.01D0
    nu=-0.324D0
    h=0.4378D0
    q=0.1254D0
    mu=1.0D0/0.13D0

    call EulerTaylor(1, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00015710714830056393842D0)
    UNIT_TEST_EQ(fpl,0.031421429660112787685D0)

    g=21.0D0/37.0D0
    call EulerTaylor(1, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00019645606369491165409D0)
    UNIT_TEST_EQ(fpl,0.039291212738982330819D0)

    g=81.0D0/29.0D0
    call EulerTaylor(1, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00014566965854639870665D0)
    UNIT_TEST_EQ(fpl,0.029133931709279741331D0)

    ! T2 hopfions
    g=1.0D0
    r=0.01D0
    nu=0.324D0
    h=0.4378D0
    q=0.1254D0
    mu=1.0D0/0.13D0

    call EulerTaylor(2, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00040402072854747751867D0)
    UNIT_TEST_EQ(fpl,0.080804145709495503734D0)

    g=21.0D0/37.0D0
    call EulerTaylor(2, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00035857525572791135990D0)
    UNIT_TEST_EQ(fpl,0.071715051145582271980D0)

    g=81.0D0/29.0D0
    call EulerTaylor(2, g, r, nu, h, q, mu, &
         NFPTS, NFPTS, FPTS, fl, fpl)
    UNIT_TEST_EQ(fl,0.00043786109431466067139D0)
    UNIT_TEST_EQ(fpl,0.087572218862932134278D0)
  end subroutine TEST_TAYLOR

  subroutine TEST_SOLVE_EULER
    use energy
    use profile
    implicit none
    integer,parameter :: NFPTSMAX=200
    integer NFPTS
    double precision FPTS(NFPTSMAX,3)
    double precision g, nu, h, q, mu
    double precision e

    ! T1 hopfions
    g     = 1.0D0
    h     = 0.0D0
    q     = 0.0D0
    mu    = 0.0D0

    nu    =-0.13D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.20165D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.0365643D0)

    ! and another quick test by computing the energy of this profile as
    ! if the magnetostatics was actually included
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, 1.0D0, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,0.0808402D0)

    ! For different aspect ratios
    g     = 1.3D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.191076D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.0050418D0)

    g     = 0.3D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.2072D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.0852203D0)

    g     = 0.1D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.201754D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.0870915D0)

    g     = 0.01D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.198104D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.086234D0)

    ! changing field
    g     = 1.0D0
    h     = 0.1D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.200977D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.100743D0)

    g     = 1.0D0
    h     =-0.1D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.200324D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,0.0200649D0)

    ! changing anisotropy
    g     = 1.0D0
    h     = 0.1D0
    q     = 0.1D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.20096D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.135024D0)

    g     = 1.0D0
    h     =-0.1D0
    q     =-0.1D0
    NFPTS = 0
    call SOLVE_EULER(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.198271D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,0.0469813D0)

    ! T2 hopfions
    g     = 1.0D0
    h     = 0.0D0
    q     = 0.0D0
    mu    = 0.0D0

    nu    = 0.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.107499D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.0626134D0)
    else
       UNIT_TEST(.false.)
    end if

    ! T2 for different aspect ratios
    g=1.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.0924209D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.0633397D0)
    else
       UNIT_TEST(.false.)
    end if

    g=0.8D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.132661D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.0635564D0)
    else
       UNIT_TEST(.false.)
    end if

    ! T2 changing field
    g=1.0D0
    h=0.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.100251D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.129043D0)
    else
       UNIT_TEST(.false.)
    end if

    h=-0.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.109913D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.00340384D0)
    else
       UNIT_TEST(.false.)
    end if

    ! T2 changing anisotropy
    h=0.1D0
    q=0.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.098581D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.163905D0)
    else
       UNIT_TEST(.false.)
    end if

    h=-0.1D0
    q=-0.1D0
    NFPTS = 0
    call SOLVE_EULER(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.107764D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,0.0251227D0)
    else
       UNIT_TEST(.false.)
    end if

  end subroutine TEST_SOLVE_EULER

  subroutine TEST_ITERATE_MS
    use energy
    use profile
    implicit none
    integer,parameter :: NFPTSMAX=200
    integer NFPTS
    double precision FPTS(NFPTSMAX,3)
    double precision g, nu, h, q, mu
    double precision e
    integer INFO

    ! T1 hopfions
    g     = 1.0D0
    h     = 0.0D0
    q     = 0.0D0
    mu    = 1.0D0

    nu    =-0.13D0
    NFPTS = 0
    INFO  = 0
    call ITERATE_MS(1, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS, INFO)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.179429D0)
    e = Etot(1, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,0.0205418D0)

    ! T2 hopfions
    g     = 1.0D0
    h     = 0.0D0
    q     = 0.0D0

    mu    = 1.0D0/3.0D0
    nu    = 0.1D0
    NFPTS = 0
    INFO  = 0
    call ITERATE_MS(2, g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS, INFO)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,0.0992666D0)
    if (NFPTS.gt.0) then
       e = Etot(2, 4.0D0*DSQRT(2.0D0), g, h, q, nu, mu, NFPTSMAX, NFPTS, FPTS)
       UNIT_TEST_EQ(e,-0.0563073D0)
    else
       CALL ASSERT(.false.)
    end if
  end subroutine TEST_ITERATE_MS

  subroutine TEST_TEST_STATE
    use test_state
    type(testing_state) :: s
    double precision xs,ys,t
    integer :: unit, ios
    call s%init(0.0D0,-3.0D0)
    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys,-3.0D0)

    open(newunit=unit, access="stream", form="unformatted", &
         status="scratch", iostat=ios)
    UNIT_TEST(ios.EQ.0)

    call s%write(unit, ios)
    UNIT_TEST(ios.EQ.0)

    call s%init(0.0D0,-3.8D0)
    t = s%evolve(0.0D0,-4.0D0)
    UNIT_TEST_EQ(t,0.980026D0)

    rewind(unit)
    call s%load(unit, ios)
    UNIT_TEST(ios.EQ.0)

    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys,-3.0D0)

    call s%init(0.0D0,-3.8D0)
    t = EVOLVE_TEST_STATE(s)
    UNIT_TEST_EQ(t,0.980026D0)

    close(unit)
  end subroutine TEST_TEST_STATE

  double precision function EVOLVE_TEST_STATE(Istate)
    use tracer
    implicit none
    class(state), intent(INOUT) :: Istate
    double precision xs,ys
    call Istate%getxy(xs,ys)
    EVOLVE_TEST_STATE = Istate%evolve(0.0D0,-4.0D0)
  end function EVOLVE_TEST_STATE

  subroutine TEST_TRACER
    use test_state
    use tracer
    implicit none

    type(testing_state) :: s

    ! grid definition
    double precision, parameter :: gridx1=-5, gridx2=5, gridy1=-5, gridy2=5
    integer, parameter :: gridnx=4, gridny=4

    integer, parameter :: NMAX = 50
    double precision :: X(NMAX), Y(NMAX)
    integer IERROR, i, N, lu, ios

    integer, parameter :: testContourN = 8
    double precision, parameter :: testContourX(testContourN) = &
         [0.0D0, 2.105655103D0, 1.049179031D0, 2.105655103D0, 0.0D0, &
         -1.937740143D0, -0.2399826271D0, -1.937740143D0]
    double precision, parameter :: testContourY(testContourN) = &
         [-3.996005279D0, -2.5D0, 0.0D0, 2.5D0, 3.996005279D0, 2.5D0, &
         0.0D0, -2.5D0]

    ! first do pure in-memory tracing
    call s%init(0.0D0,-3.0D0)
    call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, .true., &
         NMAX, N, X, Y, IERROR)
    UNIT_TEST(IERROR.EQ.0)
    UNIT_TEST_EQ(N*1.0D0,testContourN*1.0D0)
    do i = 1, min(N,testContourN)
       UNIT_TEST_EQ(X(I),testContourX(I))
       UNIT_TEST_EQ(Y(I),testContourY(I))
    end do

    ! then with log file writing (note that 1 is used instead of NMAX
    ! to make tracing fail and produce a truncated log)
    open(newunit=lu, access="stream", form="unformatted", &
         status="SCRATCH", iostat=ios)
    UNIT_TEST(ios.EQ.0)
    call s%init(0.0D0,-3.0D0)
    call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, .true., &
         1, N, X, Y, IERROR, lu)
    UNIT_TEST(IERROR.EQ.14)   ! this failure is expected

    ! should continue previous somputation
    call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, .true., &
         NMAX, N, X, Y, IERROR, lu)
    UNIT_TEST(IERROR.EQ.0)
    UNIT_TEST_EQ(N*1.0D0,testContourN*1.0D0)
    do i = 1, min(N,testContourN)
       UNIT_TEST_EQ(X(I),testContourX(I))
       UNIT_TEST_EQ(Y(I),testContourY(I))
    end do

    ! this time should just read the completed log
    call TRACE(s, gridx1, gridy1, gridx2, gridy2, gridnx, gridny, .true., &
         NMAX, N, X, Y, IERROR, lu)
    UNIT_TEST(IERROR.EQ.0)
    UNIT_TEST_EQ(N*1.0D0,testContourN*1.0D0)
    do i = 1, min(N,testContourN)
       UNIT_TEST_EQ(X(I),testContourX(I))
       UNIT_TEST_EQ(Y(I),testContourY(I))
    end do
  end subroutine TEST_TRACER

  subroutine TEST_HOPFION_STATE
    use hopf_state
    type(hopfion_state) :: s
    double precision xs,ys,t,e
    integer :: unit, ios
    call s%init(0.0D0,0.0D0)
    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys, 0.0D0)

    open(newunit=unit, access="stream", form="unformatted", &
         status="scratch", iostat=ios)
    UNIT_TEST(ios.EQ.0)

    call s%write(unit, ios)
    UNIT_TEST(ios.EQ.0)

    call s%init(0.0D0,0.0D0)
    t = s%evolve(0.0D0,-0.4D0)
    UNIT_TEST_EQ(t,1.0D0)

    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys,-0.4D0)

    t = s%evolve(0.0D0,-0.5D0)
    UNIT_TEST_EQ_WEAK(t,0.853394D0)

    rewind(unit)
    call s%load(unit, ios)
    UNIT_TEST(ios.EQ.0)

    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys, 0.0D0)

    rewind(unit)

    ! T2 hopfions
    call s%initHopfion(0.0D0,0.0D0, 2, .false.)
    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys, 0.0D0)

    e = s%energy()
    UNIT_TEST_EQ(e, -0.0626134D0)

    call s%write(unit, ios)
    UNIT_TEST(ios.EQ.0)

    t = s%evolve(0.0D0,-0.3D0)
    UNIT_TEST_EQ(t,1.0D0)

    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys,-0.3D0)

    t = s%evolve(0.0D0,-0.4D0)
    UNIT_TEST_EQ_WEAK(t,0.344864D0)

    rewind(unit)
    call s%load(unit, ios)
    UNIT_TEST(ios.EQ.0)

    call s%getxy(xs,ys)
    UNIT_TEST_EQ(xs, 0.0D0)
    UNIT_TEST_EQ(ys, 0.0D0)

    UNIT_TEST(s%HTYPE.eq.2)

    close(unit)
  end subroutine TEST_HOPFION_STATE

  ! fill FPTS array with f(r)=r^2
  subroutine SkTESTPROFILE(NFPTSMAX, NFPTS,FPTS)
    implicit none
    integer, intent(IN) :: NFPTSMAX, NFPTS
    double precision, intent(OUT) :: FPTS(NFPTSMAX,3)
    integer I
    do I=1,NFPTS
       FPTS(I,1)=(I-1.0D0)/(NFPTS-1.0D0)
       FPTS(I,2)=(1.0D0-FPTS(I,1)**2)
       FPTS(I,3)=-2.0D0*FPTS(I,1)
    end do
  end subroutine SkTESTPROFILE

  subroutine TEST_SKYRMIONS
    use skyrmions
    implicit none
    integer,parameter :: NFPTSMAX=1000
    integer,parameter :: NFPTSPFUNC=20
    integer NFPTS
    double precision FPTS(NFPTSMAX,3)
    double precision res

    double precision nu, h, q
    double precision e

    call SkTESTPROFILE(NFPTSMAX,NFPTSPFUNC,FPTS)

    res=PexSk(NFPTSMAX, NFPTSPFUNC, FPTS)
    UNIT_TEST_EQ(res,32.920805178206522403D0)

    res=PzSk(NFPTSMAX, NFPTSPFUNC, FPTS)
    UNIT_TEST_EQ(res,3.1415926535897932385D0)

    res=PaSk(NFPTSMAX, NFPTSPFUNC, FPTS)
    UNIT_TEST_EQ(res,1.5707963267948966192D0)

    res=PdmSk(NFPTSMAX, NFPTSPFUNC, FPTS)
    UNIT_TEST_EQ(res,13.698908620923286344D0)

    res=EtotSk(0.3D0,0.2D0,-0.2D0, NFPTSMAX, NFPTSPFUNC, FPTS)
    UNIT_TEST_EQ(res,-0.4933552546168910881D0)

    ! solve the Euler equation for the skyrmion profile
    h     = 0.3D0
    q     = 0.2D0
    nu    =-0.25D0
    NFPTS = 0
    call SOLVE_EULER_SK(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.254259D0)
    e = EtotSk(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.573822D0)

    h     = 0.0D0
    q     = 0.0D0
    nu    =-0.25D0
    NFPTS = 0
    call SOLVE_EULER_SK(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST(NFPTS.gt.0)
    UNIT_TEST_EQ(nu,-0.231411D0)
    e = EtotSk(h, q, nu, NFPTSMAX, NFPTS, FPTS)
    UNIT_TEST_EQ(e,-0.424442D0)

  end subroutine TEST_SKYRMIONS

  subroutine TEST_GROUNDSTATE
    use groundstate
    implicit none
    double precision, parameter :: pi = 3.1415926535897932384626433832795d0
    double precision Qcr, eHelix, stateE
    integer state

    Qcr = QCRHSMALL(0.0D0)
    UNIT_TEST_EQ(Qcr,pi*pi/4.0D0)

    Qcr = QCRHSMALL(0.3D0)
    UNIT_TEST_EQ(Qcr,1.0600373997949184685D0)

    Qcr = QCRHSMALL(pi**2/16.0D0-0.1D0)
    UNIT_TEST_EQ(Qcr,0.30995929458698632274D0)

    Qcr = QCRHSMALL(pi**2/16.0D0)
    UNIT_TEST_EQ(Qcr,0.0D0)

    Qcr = QCRHLARGE(0.71685027506808491368D0)
    UNIT_TEST_EQ(Qcr,-0.29040183377996124640D0)

    Qcr = QCRHLARGE(0.8D0)
    UNIT_TEST_EQ(Qcr,-0.51732793718587486799D0)

    Qcr = QCRHLARGE(1.0D0)
    UNIT_TEST_EQ(Qcr,-1.0D0)

    UNIT_TEST(hasStableHelices( 0.0D0 , 0.0D0))
    UNIT_TEST(hasStableHelices(-0.52D0, 0.77D0))
    UNIT_TEST(hasStableHelices(-0.4D0 ,-0.7D0))
    UNIT_TEST(hasStableHelices( 0.37D0, 0.22D0))
    UNIT_TEST(hasStableHelices( 2.04D0, 0.01D0))
    UNIT_TEST(.NOT.hasStableHelices( 0.96D0, 0.39D0))
    UNIT_TEST(.NOT.hasStableHelices(-0.24D0,-0.73D0))
    UNIT_TEST(.NOT.hasStableHelices(-0.34D0,-0.02D0))

    eHelix = energyHelix( 0.0D0, 0.0D0)
    UNIT_TEST_EQ(eHelix,-0.5D0)

    eHelix = energyHelix( 0.9D0, 0.2D0)
    UNIT_TEST_EQ(eHelix,-0.753905D0)

    eHelix = energyHelix( 2.25D0, 0.004D0)
    UNIT_TEST_EQ(eHelix,-1.15363D0)

    eHelix = energyHelix( 2.26D0, 0.004D0)
    UNIT_TEST_EQ(eHelix,-1.157131D0)

    eHelix = energyHelix( 2.35D0, 0.001D0)
    UNIT_TEST_EQ(eHelix,-1.18918D0)

    eHelix = energyHelix( -0.47D0, 0.75D0)
    UNIT_TEST_EQ(eHelix,-0.528709D0)

    eHelix = energyHelix( 1.06004D0, 0.299999D0)
    UNIT_TEST_EQ(eHelix,-0.830019D0)

    call getGroundState(0.45D0,0.006D0, 0.0D0, state, stateE)
    UNIT_TEST_EQ(state*1.0D0,3.0D0)
    UNIT_TEST_EQ(stateE,-0.615688D0)

    UNIT_TEST(isSkXGroundState(0.52D0, 0.43D0, 0.0D0))
    UNIT_TEST(.NOT.isSkXGroundState( 0.5D0, 0.2D0, 0.0D0))

  end subroutine TEST_GROUNDSTATE
  
end program TEST
