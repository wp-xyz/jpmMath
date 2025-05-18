{*******************************************************************
*        Differential equations with p variables of order 1        *
*             by Runge-Kutta method of order 4                     *
*                (Case of a Mass Pendulum)                         *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0 By Marc  *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991"    *
*             [BIBLI 03].                                          *
*                                                                  *
*             Demonstration program using procedure Equadiffp of   *
*             unit eqdifp.pas.                                     *
* ---------------------------------------------------------------- *
* EXPLANATION:                                                     *
* Let us consider the following Mass Pendulum:                     *
*                                                                  *
*                            . O                                   *
*                           /|                                     *
*                        d /_|theta                                * 
*                         /  |                                     *
*                        /   |                                     *
*                     <Mg>   |                                     *
*                                                                  *
* The deviation angle theta is determined by the differential      *
* equation:                                                        *
*             theta" = -K sin(theta) - K' theta                    *
* where:                                                           *
*             K is a coefficient depending on the pendulum         *
*             characteristics (K=Mgd/J),                           *
*             K' is another coefficient for fluid resistance.      *
*                                                                  *
* Moreover, one must know the initial conditions theta(0) and      *
* theta'(0).                                                       *
*                                                                  *
* This differential equation will be replaced by the system with   *
* two variables:                                                   *
*                  y1' = y2                                        *
*                  y'2 = -K*sin(y1) -K'*y2                         *
*                                                                  *
* (this allows calculating theta and theta' at each time step).    *
* We use a Runge-Kutta method of 4th order to solve the system     *
* with respect to time.  Here, K=(pi^2/4) and K'=0.1  (see unit    *
* Eqdifp.pas).                                                     *                   
* ---------------------------------------------------------------- *
* SAMPLE RUNS:                                                     *
*                                                                  *
* Example1: integrate system of equations from t=0 to t=12:        *
*              y1' = y2                                            *
*              y2' = -(pi^2/4)*sin(y1) - 0.1*y2                    *
* with initial conditions: y1(0)=0.25, y2(0)=0                     *
*                                                                  *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH P VARIABLE OF ORDER 1          *
*              of type yi' = f(y1,y2,...,yn), i=1..n               *
*                    (Case of a Mass Pendulum)                     *
*                                                                  *
*   Begin value time   : 0                                         *
*   End value time     : 12                                        *
*   theta value at t0  : 0.25                                      *
*   theta' value at t0 : 0                                         *
*   Number of points   : 25                                        *
*   Finesse            : 10                                        *
*                                                                  *
*     Time         Theta        Theta'                             *
* --------------------------------------                           *
*   0.000000     0.250000     0.000000                             *
*   0.500000     0.178615    -0.268789                             *
*   1.000000     0.009250    -0.372792                             *
*   1.500000    -0.157126    -0.259240                             *
*   2.000000    -0.226032    -0.004518                             *
*   2.500000    -0.163373     0.240320                             *
*   3.000000     0.010873     0.337319                             *
*   3.500000     0.140388     0.237242                             *
*   4.000000     0.204362     0.007531                             *
*   4.500000     0.149147    -0.215220                             *
*   5.000000     0.011745    -0.305188                             *
*   5.500000    -0.125654    -0.216683                             *
*   6.000000    -0.184778    -0.009451                             *
*   6.500000    -0.135956     0.193012                             *
*   7.000000    -0.012089     0.276100                             *
*   7.500000     0.112634     0.197596                             *
*   8.000000     0.167078     0.010581                             *
*   8.500000     0.123784    -0.173302                             *
*   9.000000     0.012066    -0.249776                             *
*   9.500000    -0.101090    -0.179966                             *
*  10.000000    -0.151083    -0.011147                             *
*  10.500000    -0.112595     0.155759                             *
*  11.000000    -0.011793     0.225958                             *
*  11.500000     0.090824     0.163747                             *
*  12.000000     0.136627     0.011313                             *
*                                                                  *
* Example2: integrate system of equations from t=0 to t=12:        *
*              y1' = y2                                            *
*              y2' = -(pi^2/4)*sin(y1) - 0.1*y2                    *
* with initial conditions: y1(0)=3.14, y2(0)=0                     *
*                                                                  *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH P VARIABLE OF ORDER 1          *
*              of type yi' = f(y1,y2,...,yn), i=1..n               *
*                    (Case of a Mass Pendulum)                     *
*                                                                  *
*   Begin value time   : 0                                         *
*   End value time     : 40                                        *
*   theta value at t0  : 3.14                                      *
*   theta' value at t0 : 0                                         *
*   Number of points   : 41                                        *
*   Finesse            : 10                                        *
*                                                                  *
*     Time     Theta     Theta'                                    *
* --------------------------------------                           *
*  0.000000  3.140000  0.000000                                    *
*  1.000000  3.137678 -0.005478                                    *
*  2.000000  3.124331 -0.026170                                    *
*  3.000000  3.062686 -0.120013                                    *
*  4.000000  2.781242 -0.545197                                    *
*  5.000000  1.576348 -2.134046                                    *
*  6.000000 -1.105676 -2.312060                                    *
*  7.000000 -2.175457  0.053644                                    *
*  8.000000 -1.057882  2.213567                                    *
*  9.000000  1.294094  1.644802                                    *
* 10.000000  1.684059 -0.817529                                    *
* 11.000000 -0.130730 -2.326526                                    *
* 12.000000 -1.563901 -0.249465                                    *
* 13.000000 -0.636818  1.933854                                    *
* 14.000000  1.166974  1.019504                                    *
* 15.000000  0.988226 -1.309951                                    *
* 16.000000 -0.727610 -1.446374                                    *
* 17.000000 -1.088264  0.765154                                    *
* 18.000000  0.357682  1.576103                                    *
* 19.000000  1.058234 -0.342549                                    *
* 20.000000 -0.086719 -1.524313                                    *
* 21.000000 -0.968792  0.032876                                    *
* 22.000000 -0.095485  1.389859                                    *
* 23.000000  0.858794  0.182214                                    *
* 24.000000  0.210044 -1.230693                                    *
* 25.000000 -0.748336 -0.322701                                    *
* 26.000000 -0.276785  1.074664                                    *
* 27.000000  0.646728  0.407293                                    *
* 28.000000  0.311067 -0.933079                                    *
* 29.000000 -0.557369 -0.451786                                    *
* 30.000000 -0.323904  0.809220                                    *
* 31.000000  0.480640  0.468557                                    *
* 32.000000  0.322919 -0.702785                                    *
* 33.000000 -0.415543 -0.466803                                    *
* 34.000000 -0.313287  0.612046                                    *
* 35.000000  0.360579  0.453126                                    *
* 36.000000  0.298489 -0.534862                                    *
* 37.000000 -0.314185 -0.432143                                    *
* 38.000000 -0.280842  0.469135                                    *
* 39.000000  0.274922  0.407030                                    *
* 40.000000  0.261876 -0.412985                                    *
*                                                                  *
* You may try other examples, such as:                             *
*    theta0=4, theta'0=0 from t=0 to t=20                          *
*    theta0=0, theta'0=20 from t=0 to t=10                         *
*    theta0=0, theta'0=5  from t=0 to t=10                         *
*    theta0=0, theta'0=3.15 from t=0 to t=20                       *
*                                                                  *
* You may also use output files theta.asc and theta1.asc to draw   *
* graphs of theta and theta' (not shown here), use 500 points and  *
* finesse = 1 to 4 to have a nice graph.                           *
*                                                                  *
*                            Pascal Release By J-P Moreau, Paris.  *
*                                     (www.jpmoreau.fr)            *
*******************************************************************}
Program Test_Pendulum;
Uses WinCrt, Type_def, Eqdifp;

VAR

  fi,i,ndata,p : INTEGER;
  xi,xf : REAL_AR;
  yi    : Table;
  v1,v2 : RV;

  { Exemple : y1'=y2, y2'=(pi^2/4)*sin(y1)-0.1*y2 }
  FUNCTION fp(k:INTEGER;x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    Case k of
      0: fp:=y[1];
      1: fp:=-(pi*pi/4)*sin(y[0])-0.1*y[1]
      else fp:=0
    End
  End;


{main program}
BEGIN

  New(v1); new(v2);

  Clrscr;
  writeln;
  writeln('    DIFFERENTIAL EQUATIONS WITH P VARIABLES OF ORDER 1');
  writeln('        of type yi'' = f(y1,y2,...,yn), i=1..n');
  writeln('              (Case of a Mass Pendulum)');
  writeln;

  p:=2;  {two variables}

  writeln;
  write('    Begin value time    : '); readln(xi);
  write('    End value time      : '); readln(xf);
  for i:=0 to p-1 do
  begin
    if i=0 then
    begin
      write('    theta value at t0   : ');
      readln(yi[i])
    end
    else
    begin
      write('    theta'' value at t0  : ');
      readln(yi[i])
    end
  end;
  write('    Number of points    : '); readln(ndata);
  write('    Finesse             : '); readln(fi);

  equadiffp(fp,v1,v2,xi,xf,yi,ndata,p,fi);
  Readkey;

  Dispose(v1); dispose(v2);
  DoneWinCrt

END.

{End of file pendulum.pas}