{*************************************************************************
*  Solve a boundary value problem for a second order DE of the form:     *
*  y" - p(x)*y = f(x) with conditions: y(x0)=a and y(x1)=b (x0<=x<=x1)   *
*  using a Runge-Kutta integration method.                               *
* ---------------------------------------------------------------------- *
* Explanations: if y0(x) is a particular solution of y"-p(x)y = f(x) (1) *
* and y1(x), y2(x) two lineary independant solutions of y"-p(x)y = 0 (2) *
* (right hand side removed), then the general solution of (1) has the    *
* form:  y(x) = y0(x) + c1*y1(x) + c2*y2(x), the constants c1, c2 are    *
* obtained from boundary conditions.                                     *
* A more efficient way to proceed is: 1) calculate a particular solution *
* y1(x) of (2) such as y1(x0)=0 and 2) calculate a particular solution   *
* y2(x) of (1) such as y2(x0)=a, then the general solution of (1) has the*
* form:    y(x) = y1(x) + C * y2(x), constant C is determined by:        *
*          C * y2(x1) + y1(x1) = b                                       *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
* Integrate y" - (1/x)y' - (3/x²)y = 5x² from x0=1 to x1=2 with boundary *
* conditions: y(x0)=11 and y(x1)=36. The exact solution is:              *
*           y(x) = 8/x + 2x² + x4                                        *
* The initial 2nd order problem is replaced by the following first order *
* DE systems:                                                            *
*                 1st system                2nd system                   *
*                 y' = z                    y' = z                       *
*                 z' = z/x² + 3y/x²         z' = z/x² + 3y/x² + 5x²      *
*                 with z'(x0)=0             with z'(x0)=11               *
*                                                                        *
*     X     Y estimated    Y true    Absolute Error                      *
* --------------------------------------------------                     *
*   1.0000  11.0000000   11.0000000    0.0E+0000                         *
*   1.0500  11.1498044   11.1498039    5.7E-0007                         *
*   1.1000  11.3988283   11.3988273    1.0E-0006                         *
*   1.1500  11.7472793   11.7472780    1.3E-0006                         *
*   1.2000  12.1962682   12.1962667    1.5E-0006                         *
*   1.2500  12.7476579   12.7476563    1.7E-0006                         *
*   1.3000  13.4039479   13.4039462    1.8E-0006                         *
*   1.3500  14.1681840   14.1681822    1.8E-0006                         *
*   1.4000  15.0438876   15.0438857    1.8E-0006                         *
*   1.4500  16.0349995   16.0349976    1.8E-0006                         *
*   1.5000  17.1458351   17.1458333    1.8E-0006                         *
*   1.5500  18.3810483   18.3810466    1.7E-0006                         *
*   1.6000  19.7456016   19.7456000    1.6E-0006                         *
*   1.6500  21.2447426   21.2447411    1.5E-0006                         *
*   1.7000  22.8839837   22.8839824    1.3E-0006                         *
*   1.7500  24.6690860   24.6690848    1.1E-0006                         *
*   1.8000  26.6060454   26.6060444    9.6E-0007                         *
*   1.8500  28.7010813   28.7010806    7.5E-0007                         *
*   1.9000  30.9606268   30.9606263    5.2E-0007                         *
*   1.9500  33.3913206   33.3913204    2.7E-0007                         *
*   2.0000  36.0000000   36.0000000    5.0E-0014                         *
*                                                                        *
* ---------------------------------------------------------------------- *
* REFERENCE: "Méthode de calcul numérique- Tome 2 - Programmes en Basic  *
*             et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984"*
*             [BIBLI 04].                                                *
*                                                                        *
*                                    TPW Release By J-P Moreau, Paris.   *
*                                            (www.jpmoreau.fr)           *
*************************************************************************}
Program Limits;
Uses WinCrt1;

Const SIZE = 100;

Type
      pMat = ^MAT;
      MAT = Array[0..2,0..SIZE] of Double;
      pVec = ^VEC;
      VEC = Array[0..SIZE] of Double;

Var   X: pVec;
      Y,Z: pMat;
      cc,er,xl,x0,h,ya,yb,yex: Double;
      k,kl,l: Integer;

      {y'=z}
      Function F(x,y,z:Double): Double;
      Begin
        F:=z
      End;

      {l=1: z'=z/x²+3y/x² or l=2: z'=z/x²+3y/x²+5x²}
      Function G(x,y,z:Double): Double;
      Var gg: Double;
      Begin
        gg:=z/x+3.0*y/Sqr(x);
        if l=2 then gg:=gg+5.0*x*x;
        G:=gg
      End;

      {exact solution: y=8/x + 2x3 + x4}
      Function FX(x:Double): Double;
      Begin
        FX:=8.0/x+x*x*x*(x+2.0)
      End;

      {Integrate sytem from x to x+h using Runge-Kutta}
      Procedure RK4(x,y,z,h:Double; var x1,y1,z1:Double);
      Var c1,c2,c3,c4,d1,d2,d3,d4,h2: Double;
      Begin
        c1:=F(x,y,z);
        d1:=G(x,y,z);
        h2:=h/2.0;
        c2:=F(x+h2,y+h2*c1,z+h2*d1);
        d2:=G(x+h2,y+h2*c1,z+h2*d1);
        c3:=F(x+h2,y+h2*c2,z+h2*d2);
        d3:=G(x+h2,y+h2*c2,z+h2*d2);
        c4:=F(x+h,y+h*c3,z+h*d3);
        d4:=G(x+h,y+h*c3,z+h*d3);
        x1:=x+h;
        y1:=y+h*(c1+2.0*c2+2.0*c3+c4)/6.0;
        z1:=z+h*(d1+2.0*d2+2.0*d3+d4)/6.0
      End;

{main program}
BEGIN
  {allocate memory for X,Y,Z}
  New(X); New(Y); New(Z);
  x0:=1.0;       {starting x}
  xl:=2.0;       {ending x}
  kl:=20;        {number of steps in x}
  h:=(xl-x0)/kl; {integration step}
  ya:=11.0;      {starting value for y}
  yb:=36.0;      {final value for y''}
  X^[0]:=x0;
  {integration loop}
  For l:=1 to 2 do
  begin
    Y^[l,0]:=(l-1)*ya; Z^[l,0]:=1.0;
    for k:=0 to kl-1 do
      RK4(X^[k],Y^[l,k],Z^[l,k],h,X^[k+1],Y^[l,k+1],Z^[l,k+1]);
  end;
  {solution is y0=y2+cc*y1, cc is such as y(x1)=yb}
  cc:=(yb-Y^[2,kl])/Y^[1,kl];
  {write header}
  WriteLn;
  WriteLn('     X     Y estimated    Y true    Absolute Error ');
  WriteLn(' --------------------------------------------------');
  {write kl+1 result lines}
  for k:=0 to kl do
  begin
    Y^[0,k]:=Y^[2,k] + cc*Y^[1,k];
    yex:=FX(X^[k]);       {exact value}
    er:=Abs(yex-Y^[0,k]); {current error}
    WriteLn(X^[k]:9:4,Y^[0,k]:12:7,' ',yex:12:7,'   ',er:-5)
  end;
  WriteLn(' --------------------------------------------------');
  ReadKey;
  {free memory and exit program}
  Dispose(X); Dispose(Y); Dispose(Z);
  DoneWinCrt

END.

{end of file limits.pas}