{*******************************************************************
*             Differential equations of order N                    *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991     *
*             [BIBLI 03].                                          *
*                                                                  *
*             Demonstration program using procedure Equadiffn of   *
*             unit eqdif1.pas.                                     *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
* Example: integrate y"=(2y'y'+yy)/y from x=4 to x=6               *
*          with initial conditions: y(4)=2 and y'(4)=-2tan(1)      *
*                                                                  *
* Exact solution is:   y = 2cos(1)/cos(x-5)                        *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER N          *
*              of type y(n) = f(y(n-1),...,y'',y,x)                *
*                                                                  *
*   order of equation: 2                                           *
*   begin value x    : 4                                           *
*   end value x      : 6                                           *
*   y0 value at x0   : 2                                           *
*   y1 value at x0   : -3.114815                                   *
*   number of points : 11                                          *
*   finesse          : 20                                          *
*                                                                  *
*       X            Y                                             *
* ---------------------------                                      *
*   4.000000     2.000000                                          *
*   4.200000     1.551018                                          *
*   4.400000     1.309291                                          *
*   4.600000     1.173217                                          *
*   4.800000     1.102583                                          *
*   5.000000     1.080605                                          *
*   5.200000     1.102583                                          *
*   5.400000     1.173217                                          *
*   5.600000     1.309291                                          *
*   5.800000     1.551018                                          *
*   6.000000     2.000000                                          *
*                                                                  *
*******************************************************************}
Program Test_Eqdifn;
Uses WinCrtMy, Type_def, Eqdif1;

VAR

  fi,i,ndata,ordre: INTEGER;
  xi,xf : REAL_AR;
  yi    : Table;
  vect  : RV;

  { Example: y"=(2y'y'+yy)/y }
  FUNCTION fp(x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    if abs(y[0])<Macheps then y[0]:=1e-12;
    fp:=(2*y[1]*y[1]+y[0]*y[0])/y[0]
  End;

BEGIN

  New(vect);

  Clrscr;
  writeln;
  writeln('    DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER N');
  writeln('          of type y(n) = f(y(n-1),...,y'',y,x)   ');
  writeln;
  write('    order of equation: '); readln(ordre);
  writeln;
  write('    begin value x    : '); readln(xi);
  write('    end value x      : '); readln(xf);
  for i:=0 to ordre-1 do
  begin
    write('    y',i,' value at x0   : '); readln(yi[i])
  end;
  write('    number of points : '); readln(ndata);
  write('    finesse          : '); readln(fi);

  equadiffn(fp,vect,xi,xf,yi,ndata,ordre,fi);
  Readkey;

  Dispose(vect);
  DoneWinCrt

END.

{End of file teqdifn.pas}