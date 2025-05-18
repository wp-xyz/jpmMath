{*******************************************************************
*        Differential equations with p variables of order 1        *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0 By Marc  *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991"    *
*             [BIBLI 03].                                          *
*                                                                  *
*             Demonstration program using procedure Equadiffp of   *
*             unit eqdif1.pas.                                     *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
* Example #1: integrate system of equations from x=0 to x=3:       *
*              y1' = y2 + y3 - 3*y1                                *
*              y2' = y1 + y3 - 3*y2                                *
*              y3' = y1 + y2 - 3*y3                                *
* with initial conditions: y1(0)=1, y2(0)=2 and y3(0)=-1           *
*                                                                  *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH P VARIABLE OF ORDER 1          *
*              of type yi' = f(y1,y2,...,yn), i=1..n               *
*                                                                  *
*   number of variables: 3                                         *
*   begin value x      : 0                                         *
*   end value x        : 3                                         *
*   y1 value at x0     : 1                                         *
*   y2 value at x0     : 2                                         *
*   y3 value at x0     : -1                                        *
*   number of points   : 7                                         *
*   finesse            : 30                                        *
*                                                                  *
*       X        Y1         Y2         Y3                          *
* --------------------------------------------------------         *
*   0.000000   1.000000   2.000000  -1.000000                      *
*   0.500000   0.449466   0.584801   0.178795                      *
*   1.000000   0.251358   0.269674   0.214727                      *
*   1.500000   0.149580   0.152058   0.144622                      *
*   2.000000   0.090335   0.090671   0.089664                      *
*   2.500000   0.054738   0.054784   0.054648                      *
*   3.000000   0.033193   0.033200   0.033181                      *
* --------------------------------------------------------         *
*                                                                  *
* Example #2: integrate system of equations from x=0 to PI*Sqrt(2) *
*              y1' = y2                                            *
*              y2' = -4*y1 - 3*y3                                  *
*              y3' = y4                                            *
*              y4' = -8*y1 - 2*y3                                  *
* with initial conditions: y1(0)=3, y2(0)=0, y3(0)=4, y4(0)=0      *
*                                                                  *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH P VARIABLE OF ORDER 1          *
*              of type yi' = f(y1,y2,...,yn), i=1..n               *
*                                                                  *
*   number of variables: 4                                         *
*   begin value x      : 0                                         *
*   end value x        : 4.442883                                  *
*   y1 value at x0     : 3                                         *
*   y2 value at x0     : 0                                         *
*   y3 value at x0     : 4                                         *
*   y4 value at x0     : 0                                         *
*   number of points   : 9                                         *
*   finesse            : 30                                        *
*                                                                  *
*       X        Y1         Y2         Y3         Y4               *
* --------------------------------------------------------         *
*   0.000000   3.000000   0.000000   4.000000   0.000000           *
*   0.555360   0.000000  -8.485281   0.000000 -11.313708           *
*   1.110721  -3.000000  -0.000001  -4.000000  -0.000002           *
*   1.666081  -0.000001   8.485281  -0.000001  11.313708           *
*   2.221442   3.000000   0.000003   4.000000   0.000003           *
*   2.776802   0.000001  -8.485281   0.000002 -11.313708           *
*   3.332162  -3.000000  -0.000004  -4.000000  -0.000005           *
*   3.887523  -0.000002   8.485281  -0.000002  11.313708           *
*   4.442883   3.000000   0.000005   4.000000   0.000007           *
* --------------------------------------------------------         *
*                                                                  *
* Release 1.1: Added example #2 (Sept. 2008).                      *
*                                                                  *
*                                               J-P Moreau, Paris. *
*                                               (www.jpmoreau.fr)  *
*******************************************************************}
Program Test_Eqdifp;
Uses WinCrtMy, Type_def, Eqdif1;

VAR

  fi,i,ndata,p : INTEGER;
  xi,xf : REAL_AR;
  yi    : Table;
  v1,v2,v3,v4: RV;

  { Exemple #1: y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3
  FUNCTION fp(k:INTEGER;x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    Case k of
      0: fp:=y[1]+y[2]-3*y[0];
      1: fp:=y[0]+y[2]-3*y[1];
      2: fp:=y[0]+y[1]-3*y[2];
      else fp:=0
    End
  End;    }

  { Exemple #2: y1'=y2, y2'=-4y1-3y3, y'3=y4, y'4=-8y1-2y3 }
  FUNCTION fp(k:INTEGER;x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    Case k of
      0: fp:=y[1];
      1: fp:=-4*y[0]-3*y[2];
      2: fp:=y[3];
      3: fp:=-8*y[0]-2*y[2]
      else fp:=0.0
    End
  End;

BEGIN

  New(v1); new(v2); new(v3); new(v4);

  Clrscr;
  writeln;
  writeln('    DIFFERENTIAL EQUATIONS WITH P VARIABLES OF ORDER 1');
  writeln('        of type yi'' = f(y1,y2,...,yn), i=1..n');
  writeln;
  write('    number of variables : '); readln(p);
  writeln;
  write('    begin value x       : '); readln(xi);
  write('    end value x         : '); readln(xf);
  for i:=0 to p-1 do
  begin
    write('    y',i+1,' value at x0      : '); readln(yi[i])
  end;
  write('    number of points    : '); readln(ndata);
  write('    finesse             : '); readln(fi);

  equadiffp(fp,v1,v2,v3,v4,xi,xf,yi,ndata,p,fi);
  Readkey;

  Dispose(v1); dispose(v2); Dispose(v3); dispose(v4);
  DoneWinCrt

END.

{End of file teqdifp.pas}