{*******************************************************************
*        Differential equations with p variables of order 1        *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0 By Marc  *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991"    *
*             [BIBLI 03].                                          *
*                                                                  *
*             Demonstration program using procedure Equadiffp of   *
*             unit eqdif1.pas with graph option.                   *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
* Example: integrate system of equations from x=0 to x=3:          *
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
*   finesse            : 100                                       *
*                                                                  *
*       X            Y1           Y2                               *
* --------------------------------------                           *
*   0.000000     1.000000     2.000000                             *
*   0.500000     0.449466     0.584801                             *
*   1.000000     0.251358     0.269674                             *
*   1.500000     0.149580     0.152058                             *
*   2.000000     0.090335     0.090671                             *
*   2.500000     0.054738     0.054784                             *
*   3.000000     0.033193     0.033200                             *
*                                                                  *
*******************************************************************}
Program Test_Eqdifp;
Uses WinCrtMy, Winprocs, Type_def, Eqdif1, Graph_2D;

VAR

  fi,i,ndata,p : INTEGER;
  xi,xf : REAL_AR;
  yi    : Table;
  v1,v2 : RV;

  { Exemple : y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3 }
  FUNCTION fp(k:INTEGER;x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    Case k of
      0: fp:=y[1]+y[2]-3*y[0];
      1: fp:=y[0]+y[2]-3*y[1];
      2: fp:=y[0]+y[1]-3*y[2];
      else fp:=0
    End
  End;

BEGIN

  New(v1); new(v2);

  WinCrtInit(' DIFFERENTIAL EQUATION WITH P VARIABLES OF ORDER 1');

  Repeat
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

    equadiffp(fp,v1,v2,xi,xf,yi,ndata,p,fi);

    Write(' Do you want a graph of solution (y/n) ? ');
    if Readkey='y' then
    begin
      ClrScr;
      CourbeXY(CrtDC,ndata,5,v1,xi,xf);
      Legendes(CrtDC,'SOLUTION Y1 = F(X)','X','Y');
      CourbeXY(CrtDC,ndata,3,v2,xi,xf);
      Legendes(CrtDC,'SOLUTION Y2 = F(X)','X','Y');
      TextOut(CrtDC,MaxX-250,50,'y1'' = y2 + y3 - 3*y1',20);
      TextOut(CrtDC,MaxX-250,65,'y2'' = y1 + y3 - 3*y2',20);
      TextOut(CrtDC,MaxX-250,80,'y3'' = y1 + y2 - 3*y3',20);
      TextOut(CrtDC,MaxX-240,95,'from x=0 to x=3',15);
      SortieGraphique
    end
  Until rep='n';

  Dispose(v1); dispose(v2);
  DoneWinCrt

END.

{End of file teqdifp.pas}