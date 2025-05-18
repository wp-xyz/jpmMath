{*******************************************************************
*             Differential equations of order 1                    *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991     *
*             [BIBLI 03].                                          *
*                                                                  *
*             Demonstration program using procedure Equadiff1 of   *
*             unit eqdif1.pas with graph option.                   *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
* Example: integrate y'=4x*(y+sqrt(y))/1+x² from x=0 to x=1        *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1          *
*                     of type y' = f(x,y)                          *
*                                                                  *
*   begin value x   : 0                                            *
*   end value x     : 1                                            *
*   y value at x0   : 1                                            *
*   number of points: 11                                           *
*   finesse         : 10                                           *
*                                                                  *
*       X            Y                                             *
* ---------------------------                                      *
*   0.000000     1.000000                                          *
*   0.100000     1.040400                                          *
*   0.200000     1.166400                                          *
*   0.300000     1.392400                                          *
*   0.400000     1.742400                                          *
*   0.500000     2.250000                                          *
*   0.600000     2.958400                                          *
*   0.700000     3.920400                                          *
*   0.800000     5.198400                                          *
*   0.900000     6.864400                                          *
*   1.000000     9.000000                                          *
*                                                                  *
*******************************************************************}
Program Test_Eqdif1;
Uses WinCrtMy, WinProcs, Type_def, Eqdif1, Graph_2D;

VAR

  fi,ndata:INTEGER;
  xi,xf,yi:REAL_AR;
  vect    :RV;

  {Example: y'=4x(y+root(y))/(1+x²) }
  FUNCTION fp(x,y:REAL_AR):REAL_AR; FAR;
  Begin
    fp:=4*x*(y+sqrt(y))/(1+x*x)
  End;

{main program}
BEGIN

  New(vect);

  WinCrtInit(' DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1');

  Repeat
    Clrscr;
    writeln;
    writeln('    DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER 1');
    writeln('             of type y'' = f(x,y)');
    writeln;
    write('    begin value x   : '); readln(xi);
    write('    end value x     : '); readln(xf);
    write('    y value at x0   : '); readln(yi);
    write('    number of points: '); readln(ndata);
    write('    finesse         : '); readln(fi);

    equadiff1(fp,vect,xi,xf,yi,ndata,fi);

    Write(' Do you want a graph of solution (y/n) ? ');
    if Readkey='y' then
    begin
      ClrScr;
      CourbeXY(CrtDC,ndata,10,vect,xi,xf);
      Legendes(CrtDC,'SOLUTION Y = F(X)','X','Y');
      TextOut(CrtDC,MaxX-250,50,'Equation: y''=4x(y+root(y))/(1+x²)',33);
      TextOut(CrtDC,MaxX-200,70,'from x=0 to x=1',15);
      SortieGraphique
    end
  Until rep='n';

  Dispose(vect);
  DoneWinCrt

END.

{End of file teqdif1.pas}