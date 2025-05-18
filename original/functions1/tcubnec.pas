{****************************************************************
* Test program for cubature over rectangles using Newton-Cotes  *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
* (Integrate function EXP(-(x*x + y*y)) in rectangle defined by *
*  x1=0, x2=10, y1=0, y2=10).                                   *
*                                                               *
* -----------------------------------------------------------   *
*   #    ERROR     RESULT        ERROR     NUMBER OF FUNCTION   *
* METHOD CODE                  ESTIMATE        EVALUATIONS      *
* -----------------------------------------------------------   *
* Without error estimation:                                     *
*   1    ( 0)   0.7855606650                        121         *
*   2    ( 0)   0.7853439999                        441         *
*   3    ( 0)   0.7853778519                        961         *
*   4    ( 0)   0.7854017744                       1681         *
*   5    ( 0)   0.7854001104                       2601         *
*   6    ( 0)   0.7853979700                       3721         *
*   7    ( 0)   0.7853980474                       5041         *
* With error estimation:                                        *
*   1    ( 0)   0.7853981634    5.4E-0005           562         *
*   2    ( 0)   0.7853981634   -3.6E-0006          2122         *
*   3    ( 0)   0.7853981634   -1.4E-0006          4682         *
*   4    ( 0)   0.7853981634    5.7E-0008          8242         * 
*   5    ( 0)   0.7853981634    3.1E-0008         12802         *
*   6    ( 0)   0.7853981634   -7.6E-0010         18362         *
*   7    ( 0)   0.7853981634   -4.5E-0010         24922         *
* -----------------------------------------------------------   *
*                                                               *
* Reference:                                                    *
* "Numerical Algorithms with C,  By Gisela Engeln-Muellges      *
*  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].          *
*                                                               *
*                         TPW Release 1.0 By J-P Moreau, Paris. *
*                                  (www.jpmoreau.fr)            *
****************************************************************}
Program TKUBNEC;

Uses WinCrt, Kubnec;


Var

  a, b, c, d, W, F: Double;
  i, Nx, Ny, method, nF: Integer; 
  mSCH: Integer;                        {with error estimation, if <>0}

                    {Function func(x,y) is declared in unit kubnec.pas}

BEGIN

  a := 0.0; c := 0.0;          {define limits of integration rectangle}
  b := 10.0; d := 10.0;

  Nx:=10; Ny:=10; 

  {print header}
  writeln('-----------------------------------------------------------');
  writeln('  #    ERROR     RESULT        ERROR     NUMBER OF FUNCTION');
  writeln('METHOD CODE                  ESTIMATE        EVALUATIONS   ');
  writeln('-----------------------------------------------------------');
  
  {main loop}
  for mSCH := 0 to 1 do
  begin
    if mSCH<>0 then
      writeln('With error estimation:')
    else
      writeln('Without error estimation:');
    {integrate with 7 methods}
    for method := 1 to 7 do
    begin

      i := Kub4NeCn (a, b, Nx, c, d, Ny, method, W, mSCH, F); {see kubnec.pas}

      write(method:3,'    (',i:2,')  ',W:13:10);

      if mSCH<>0 then write('   ',F:8)
      else write('             ');

      nF := (Nx * method + 1) * (Ny * method + 1);

      if mSCH<>0 then nF := nF + (2 * Nx * method + 1)*(2 * Ny * method + 1);
      if F>0.0 then writeln('       ',nF:6)
      else writeln('       ',nF:6)
    end
  end;
  writeln('-----------------------------------------------------------');

  ReadKey;
  DoneWinCrt

END.

{end of file tkubnec.pas}