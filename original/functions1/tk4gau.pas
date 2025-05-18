{****************************************************************
*    Test program for cubature over rectangles using  Gauss     *
*                                                               *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
* (Integrate function EXP(-(x*x + y*y)) in rectangle defined by *
*  x1=0, x2=10, y1=0, y2=10).                                   *
*                                                               *
* ---------------------------------------------------           * 
*   #      Error     Value       Error      Function            *
* Method   code     Integral   Estimation     calls             *
* ---------------------------------------------------           *
*   0      ( 0)    0.27460788                   16              *
*   1      ( 0)    0.94345609                   64              *
*   2      ( 0)    0.77390367                  144              *
*   3      ( 0)    0.78537787                  256              *
*   4      ( 0)    0.78547216                  400              *
*   5      ( 0)    0.78539063                  576              *
*   6      ( 0)    0.78539862                  784              *
*   7      ( 0)    0.78539814                 1024              *
*   0      ( 0)    0.77973416   1.684e-001      80              *
*   1      ( 0)    0.78676409  -5.223e-002     320              * 
*   2      ( 0)    0.78527076   3.789e-003     720              *
*   3      ( 0)    0.78540422   8.783e-006    1280              *
*   4      ( 0)    0.78539799  -2.472e-005    2000              *
*   5      ( 0)    0.78539817   2.513e-006    2880              * 
*   6      ( 0)    0.78539816  -1.514e-007    3920              *
*   7      ( 0)    0.78539816   6.261e-009    5120              *
* ---------------------------------------------------           *
*                                                               *
* ------------------------------------------------------------- *
* Reference:                                                    *
* "Numerical Algorithms with C,  By Gisela Engeln-Muellges      *
*  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].          *
*                                                               *
*                         TPW Release 1.0 By J-P Moreau, Paris. *
*                                  (www.jpmoreau.fr)            *
****************************************************************}
Uses WinCrt, KubGauss;

Var
  
  a, b, c, d: Double;
  W, F: Double;
  Nx, Ny: Integer;
  i, Verfahren: Integer;
  mSCH: Integer;                         {<>0 = with eror estimate}


BEGIN

  Nx := 4; Ny := 4;
  a := 0.0; c := 0.0;
  b := 10.0; d := 10.0;

  writeln(' -------------------------------------------------------');
  writeln('   #      Error     Value          Error       Function ');
  writeln(' Method   code     Integral      Estimation      calls  ');
  writeln(' -------------------------------------------------------');

  for mSCH := 0 to 1 do
    for Verfahren := 0 to 7 do
    begin
     
      fcalls := 0;

      i := Kub4GauE(a, b, Nx, c, d, Ny, Verfahren, W, mSCH, F);

      write('   ',Verfahren,'      (',i:2,')  ',W:12:8);

      if mSCH<>0 then
        write('    ',F:12)
      else
        write('                ');

      writeln('     ',fcalls:4)
    end;
    writeln(' -------------------------------------------------------');
    ReadKey;
    DoneWinCrt
  
END.

{end of file tk4gau.pas}