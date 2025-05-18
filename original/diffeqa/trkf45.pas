{**************************************************************************
*   Integrate a System of Ordinary Differential Equations By the          *
*   Runge-Kutta-Fehlberg method (double precision)                        *
* ----------------------------------------------------------------------- *
* REFERENCE:     H A Watts and L F Shampine,                              *
*                Sandia Laboratories,                                     *
*                Albuquerque, New Mexico.                                 *
* ----------------------------------------------------------------------- *
* SAMPLE RUN:                                                             *
*                                                                         *
* PROGRAM TRKF45                                                          *
* Demonstrate the RKF45 ODE integrator.                                   *
*                                                                         *
* TEST01                                                                  *
* Solve a scalar equation:                                                *
*                                                                         *
*  Y' = 0.25 * Y * ( 1 - Y / 20 )                                         *
*                                                                         *
*       T          Y       Y_Exact    Error                               *
*                                                                         *
*    0.00000    1.00000    1.00000    0.0000000                           *
*    4.00000    2.50322    2.50322   -0.0000008                           *
*    8.00000    5.60009    5.60009   -0.0000020                           *
*   12.00000   10.27773   10.27773   -0.0000018                           *
*   16.00000   14.83682   14.83683   -0.0000104                           *
*   20.00000   17.73016   17.73017   -0.0000029                           *
*                                                                         *
* TEST02                                                                  *
* Solve a vector equation:                                                *
*                                                                         *
*  Y'(1) =   Y(2)                                                         *
*  Y'(2) = - Y(1)                                                         *
*                                                                         *
*       T            Y1            Y2                                     *
*                                                                         *
*    0.00000      1.000000      0.000000                                  *
*    0.52360      0.866025     -0.500000                                  *
*    1.04720      0.500000     -0.866026                                  *
*    1.57080      0.000000     -1.000000                                  *
*    2.09440     -0.500000     -0.866026                                  *
*    2.61799     -0.866026     -0.500001                                  *
*    3.14159     -1.000001     -0.000000                                  *
*    3.66519     -0.866027      0.500000                                  *
*    4.18879     -0.500001      0.866026                                  *
*    4.71239     -0.000000      1.000001                                  *
*    5.23599      0.500001      0.866027                                  *
*    5.75959      0.866027      0.500001                                  *
*    6.28319      1.000002      0.000000                                  *
*                                                                         *
* TEST03                                                                  *
* Solve a vector equation:                                                *
*                                                                         *
*  Y'(1) = Y(2)                                                           *
*  Y'(2) = Y(3)                                                           *
*  Y'(3) = Y(4)                                                           *
*  Y'(4) = Y(5)                                                           *
*  Y'(5) = (45 * Y([3) * Y(4) * Y(5) - 40 * Y(4)^3 / (9 * Y(3)^2)         *
*                                                                         *
*       T          Y1          Y2          Y3          Y4          Y5     *
*                                                                         *
*    0.00000    1.000000    1.000000    1.000000    1.000000    1.000000  *
*    0.13636    1.146098    1.146091    1.145868    1.140683    1.056035  *
*    0.27273    1.313535    1.313400    1.311282    1.285316    1.052479  *
*    0.40909    1.505384    1.504599    1.496117    1.423333    0.952094  *
*    0.54545    1.725082    1.722230    1.698438    1.538585    0.711111  *
*    0.68182    1.976382    1.968402    1.913699    1.608970    0.288089  *
*    0.81818    2.263276    2.244383    2.133997    1.607806   -0.339181  *
*    0.95455    2.589841    2.550107    2.347698    1.508010   -1.150271  *
*    1.09091    2.960029    2.883687    2.539850    1.289462   -2.060935  *
*    1.22727    3.377386    3.241053    2.693755    0.948202   -2.920461  *
*    1.36364    3.844750    3.615889    2.793723    0.503786   -3.543020  *
*    1.50000    4.363961    4.000000    2.828427   -0.000000   -3.771236  *
*                                                                         *
* ----------------------------------------------------------------------- *
*  NOTE: Uses unit urkf45.pas.                                            *
*                                                                         *
*                                 TPW Release 1.1 By J-P Moreau, Paris.   *
*                                          (www.jpmoreau.fr)              *
* ----------------------------------------------------------------------- *
* Release 1.1 (03/15/05): added TEST #3.                                  *
**************************************************************************}
program trkf45;

Uses WinCrt, URKF45;


function y_exact_01(t:double): double;
{--------------------------------------------------------------------
!
!  Y_EXACT_01 evaluates the exact solution of the ODE (TEST #1).
!
!-------------------------------------------------------------------}
Begin
  y_exact_01 := 20.0 / ( 1.0 + 19.0 * exp(- 0.25 * t))
End;


Procedure test01;
{--------------------------------------------------------------------
!
!  TEST01 solves a scalar ODE in double precision.
!
!-------------------------------------------------------------------}
Var
    neqn: integer;
    abserr:double;
    i_step, iflag, n_step: integer;
    relerr, t, t_out, t_start, t_stop: double;
    y:VEC;
Begin
  num:=1;    {example #1}
  neqn:=1;   {one equation}
  writeln;
  writeln(' TEST01');
  writeln(' Solve a scalar equation:');
  writeln;
  writeln('  Y'' = 0.25 * Y * ( 1 - Y / 20 )');
  writeln;

  abserr := 0.000001e+00;
  relerr := 0.000001e+00;

  iflag := 1;

  t_start := 0.0;
  t_stop := 20.0;

  n_step := 5;

  t_out := 0.0;
  y[1] := 1.0;

  writeln('        T           Y        Y_Exact     Error');
  writeln;
  writeln(t_out:12:5, y[1]:12:5, y_exact_01(t_out):12:5, (y[1] - y_exact_01(t_out)):14:7);

  for i_step := 1 to n_step do
  begin

    t := ((n_step - i_step + 1) * t_start + (i_step - 1) * t_stop) / (1.0*n_step);
    t_out := ( ( n_step - i_step ) * t_start + (i_step) * t_stop) / (1.0*n_step);

    rkf45 (neqn, y, t, t_out, relerr, abserr, iflag);

    writeln(t_out:12:5, y[1]:12:5, y_exact_01(t_out):12:5, (y[1] - y_exact_01(t_out)):14:7)

  end
End;

Procedure test02;
{--------------------------------------------------------------------
!
!  TEST02 solves a vector ODE.
!
!-------------------------------------------------------------------}
Var
    neqn: integer;
    abserr:double;
    i_step, iflag, n_step: integer;
    relerr, t, t_out, t_start, t_stop: double;
    y:VEC;
Begin

  num:=2;       {Example #2}
  neqn := 2;    {2 equations}

  writeln;
  writeln(' TEST02');
  writeln(' Solve a vector equation:');
  writeln;
  writeln('  Y''(1) =   Y(2)');
  writeln('  Y''(2) = - Y(1)');

  abserr := 0.000001;
  relerr := 0.000001;

  iflag := 1;

  t_start := 0.0;
  t_stop := 2.0 * PI;

  n_step := 12;

  t_out := 0.0;

  y[1] := 1.0;
  y[2] := 0.0;

  writeln;
  writeln('        T          Y1          Y2');
  writeln;
  writeln(t_out:12:5, y[1]:12:6, y[2]:12:6);

  for i_step := 1 to n_step do
  begin

    t := ((n_step - i_step + 1) * t_start + (i_step - 1) * t_stop) / (1.0*n_step);
    t_out := ((n_step - i_step) * t_start + (i_step  * t_stop) / (1.0* n_step));

    rkf45 (neqn, y, t, t_out, relerr, abserr, iflag);

    writeln(t_out:12:5, y[1]:12:6, y[2]:12:6)

  end
End;  {test02}

Procedure test03;
{--------------------------------------------------------------------
!
!  TEST03 solves a vector ODE.
!
!-------------------------------------------------------------------}
Var
    neqn: integer;
    abserr:double;
    i_step, iflag, n_step: integer;
    relerr, t, t_out, t_start, t_stop: double;
    y:VEC;
Begin

  num:=3;       {Example #3}
  neqn := 5;    {5 equations}

  writeln;
  writeln(' TEST03');
  writeln(' Solve a vector equation:');
  writeln;
  writeln(' Y''(1) = Y(2)');
  writeln(' Y''(2) = Y(3)');
  writeln(' Y''(3) = Y(4)');
  writeln(' Y''(4) = Y(5)');
  writeln(' Y''(5) = (45 * Y([3) * Y(4) * Y(5) - 40 * Y(4)^3 / (9 * Y(3)^2)');

  abserr := 1e-8;
  relerr := 1e-10;

  iflag := 1;

  t_start := 0.0;
  t_stop  := 1.5;

  n_step := 11;

  t_out := 0.0;

  {initial conditions}
  y[1] := 1.0;
  y[2] := 1.0;
  y[3] := 1.0;
  y[4] := 1.0;
  y[5] := 1.0;

  writeln;
  writeln('        T          Y1          Y2          Y3          Y4          Y5');
  writeln;
  writeln(t_out:12:5, y[1]:12:6, y[2]:12:6, y[3]:12:6, y[4]:12:6, y[5]:12:6);

  for i_step := 1 to n_step do
  begin

    t := ((n_step - i_step + 1) * t_start + (i_step - 1) * t_stop) / (1.0*n_step);
    t_out := ((n_step - i_step) * t_start + (i_step  * t_stop) / (1.0* n_step));

    rkf45 (neqn, y, t, t_out, relerr, abserr, iflag);
    writeln(t_out:12:5, y[1]:12:6, y[2]:12:6, y[3]:12:6, y[4]:12:6, y[5]:12:6);

  end
End;  {test03}


{main program}
BEGIN

  writeln;
  writeln(' PROGRAM TRKF45');
  writeln(' Demonstrate the RKF45 ODE integrator.');

  test01;
  readkey;
  test02;
  readkey;
  test03;
  readkey;
  DoneWinCrt

END.

{end of file trkf45.pas (J-P Moreau Release 1.1, 03/15/05)