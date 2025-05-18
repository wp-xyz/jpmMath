{****************************************************************
*    Test program for cubature over triangles using 3-point     *
*    Newton-Cotes and Romberg-Richardson extrapolation          *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
* (Integrate function f(x,y)=xy over triangle defined by three  *
*  points (0, 0), (10, 0) and (0, 10) ).                        *
*                                                               *
*  -----------------------------------------------------------  * 
*  # Newton-  Error       Value              Error              *
*    Cotes    Code       Integral          Estimation           *
*  -----------------------------------------------------------  *
*     1         0       833.333333     0.00000000000000E+0000   *
*     2         0       472.222222    -2.25694444444445E+0001   *
*     3         0       422.839506    -7.71604938271594E-0001   *
*     4         0       417.417090    -2.11813120309898E-0002   *
*     5         0       416.759828    -6.41857940308910E-0004   *
*     6         0       416.678292    -1.99062169485842E-0005   *
*     7         0       416.668119    -6.20892194547196E-0007   *
*  -----------------------------------------------------------  *
*                                                               *
* ------------------------------------------------------------- *
* Reference:                                                    *
*   "Numerical Algorithms with C,  By Gisela Engeln-Muellges    *
*    and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].        *
*                                                               *
*                         TPW Release 1.0 By J-P Moreau, Paris. *
*                                  (www.jpmoreau.fr)            *       
****************************************************************}
Program TK3NEC;

Uses WinCrt, Kubnec;

Var
    W, F: Double;
    i, j: Integer;


{main program}
BEGIN

  writeln;
  writeln('-----------------------------------------------------------');
  writeln(' # Newton-  Error      Value              Error            ');
  writeln('   Cotes    Code      Integral          Estimation         ');
  writeln('-----------------------------------------------------------');

  for j := 1 to 7 do
  begin

    {              Px,  Py,  Qx,  Qy,   Rx,  Ry  }
    i := Kub3RoRi(0.0, 0.0, 10.0, 0.0, 0.0, 10.0, j, W, F);

    write('   ',j:2,'        ',i:2);
    writeln('    ',W:12:6,'    ', F)

  end;
  writeln('-----------------------------------------------------------');

  ReadKey;
  DoneWinCrt

END.

{end of file tk3nec.pas}