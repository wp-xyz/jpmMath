{******************************************************************
* Solve an ordinary system of first order differential equations  *
* using the predictor-corrector method of Adams-Bashforth-Moulton *
* --------------------------------------------------------------- *
* SAMPLE RUNS:                                                    *
* Example #1                                                      *
*                                                                 *
* # example (0 to 3) ..........: 0                                *
* Number of equations in system: 2                                *
* Starting x ..................: 0                                * 
* Ending x ....................: 2                                *
* Absolute precision ..........: 1e-8                             *
* Relative precision ..........: 1e-10                            * 
* Starting integration step ...: 0.1                              *
* Max. number of function calls: 1000                             *
*                                                                 *
* Y(0) = 0                                                        *
* Y(1) = 0                                                        *
*                                                                 *
* D.E. System to solve:                                           *
* y1' = y1 * y2 + cos(x) - 0.5 * sin(2.0*x)                       *
* y2' = y1 * y1 + y2 * y2 - (1 + sin(x))                          *
*                                                                 *
* x=  2.00000000000000E+0000                                      *
* Y(0) =  9.19316247082538E-0002                                  *
* Y(1) = -1.36385503717724E+0000                                  *
*                                                                 *
* Number of function calls: 566                                   *
* Final step size ........:  8.33333333333330E-0003               *
* Error code .............: 0                                     *
*                                                                 *
* Example #2                                                      *
*                                                                 *
* # example (0 to 3) ..........: 1                                *
* Number of equations in system: 1                                *
* Starting x ..................: 0                                * 
* Ending x ....................: 1                                *
* Absolute precision ..........: 1e-8                             *
* Relative precision ..........: 1e-10                            * 
* Starting integration step ...: 0.1                              *
* Max. number of function calls: 500                              *
*                                                                 *
* Y(0) = 1                                                        *
*                                                                 *
* D.E. System to solve:                                           *
* y' = -y + x/((1+x)*(1+x))                                       *
*                                                                 *
* x=  1.00000000000000E+0000                                      *
* Y(0) =  4.99999998141136E-0001                                  *
*                                                                 *
* Number of function calls: 149                                   *
* Final step size ........:  3.33333333333333E-0002               *
* Error code .............: 0                                     *
*                                                                 *
* Example #3                                                      *
*                                                                 *
* # example (0 to 3) ..........: 3                                *
* Number of equations in system: 5                                *
* Starting x ..................: 0                                * 
* Ending x ....................: 2                                *
* Absolute precision ..........: 1e-8                             *
* Relative precision ..........: 1e-10                            * 
* Starting integration step ...: 0.1                              *
* Max. number of function calls: 2000                             *
*                                                                 *
* Y(0) = 1                                                        *
* Y(1) = 1                                                        *
* Y(2) = 1                                                        *
* Y(3) = 1                                                        *
* Y(4) = 1                                                        *
*                                                                 *
* D.E. System to solve:                                           *
* y1' = y2                                                        *
* y2' = y3                                                        * 
* y3' = y4                                                        *
* y4' = y5                                                        * 
* y5' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)   *
*                                                                 *
* x=  2.00000000000000E+0000                                      *
* Y(0) =  6.70820393221835E+0000                                  *
* Y(1) =  5.34164078505024E+0000                                  *
* Y(2) =  2.41495341322099E+0000                                  *
* Y(3) = -1.44897204948409E+0000                                  *
* Y(4) = -1.44897204334884E+0000                                  *
*                                                                 *
* Number of function calls: 1075                                  *
* Final step size ........:  1.04166666666659E-0002               *
* Error code .............: 0                                     *
*                                                                 *
* --------------------------------------------------------------- *
*                          TPW Release 1.1 By J-P Moreau, Paris.  *
*                                    (www.jpmoreau.fr)            *
*                                                                 *
* Release 1.1 (12/15/06): added two more examples.                *                                                            
******************************************************************}
Program Test_Praekorr;

Uses WinCrt1,      {Borland Turbo-pascal Unit for writeln, etc.}
     Fgauss,       {for type pVec}
     Ab_mou,       {for procedure prae_korr [BIBLI 11]}
     t_dgls;       {for bspnummer}

Var
    xbegin,xend,epsabs,epsrel,h: double;
    i,n,fmax,ncalls,error: integer;
    Y: pVEC;

BEGIN

  New(Y);

  writeln;                                   {# example - see t_dgls.pas}
  write(' # example (0 to 3) ..........: '); readln(bspnummer);
  write(' Number of equations in system: '); readln(n);
  write(' Starting x ..................: '); readln(xbegin);
  write(' Ending x ....................: '); readln(xend);
  write(' Absolute precision ..........: '); readln(epsabs);
  write(' Relative precision ..........: '); readln(epsrel); 
  write(' Starting integration step ...: '); readln(h);
  write(' Max. number of function calls: '); readln(fmax);
  writeln;
  for i:=0 to n-1 do  {initial value of y}
  begin
    write(' Y(',i,') = '); readln(Y^[i])
  end;         

  {call integration procedure using the predictor-corrector method
   of Adams-Bashforth-Moulton }

  error := prae_korr(xbegin,Y,n,xend,h,epsabs,epsrel,fmax,ncalls,1.0,True);

  {write results}
  writeln;
  writeln(' D.E. System to solve:');
  for i:=0 to n-1 do writeln(dgltxt[i]);
  writeln;
  writeln(' x= ', xend);
  for i:=0 to n-1 do
    writeln(' Y(',i,') = ', Y^[i]);
  writeln;
  writeln(' Number of function calls: ', ncalls);
  writeln(' Final step size ........: ', h);
  writeln(' Error code .............: ', error);

  Dispose(Y);
  ReadKey;
  DoneWinCrt

END.

{end of file tabmou.pas}