{*******************************************************************
*    Solve Y' = F(X,Y) with Initial Condition Y(X0)=Y0 using the   *
*    Adams-Bashforth Method                                        *
* ---------------------------------------------------------------- *
* REFERENCE: "Méthode de calcul numérique- Tome 2 - Programmes en  *
*             Basic et en Pascal By Claude Nowakowski, Edition du  *
*             P.S.I., 1984" [BIBLI 04].                            * 
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
* (Solve Y' = -Y + X/((1+X)*(1+X))  with Y(0)=1 for X0=0 up to     *
*  X1=1.0, exact solution is Y = 1 / (1+X).                        *
*                                                                  *
*      X           Y        Y exact      Error                     *
* ------------------------------------------------                 *
*   0.000000    1.000000    1.000000   0                           *
*   0.050000    0.952381    0.952381   7.2E-0010                   *
*   0.100000    0.909091    0.909091   1.9E-0009                   *
*   0.150000    0.869525    0.869565   4.0E-0005                   *
*   0.200000    0.833265    0.833333   6.8E-0005                   *    
*   0.250000    0.799910    0.800000   9.0e-0005                   *
*   0.300000    0.769125    0.769231   1.1e-0004                   *
*   0.350000    0.740623    0.740741   1.2e-0004                   *
*   0.400000    0.714160    0.714286   1.3e-0004                   *
*   0.450000    0.689525    0.689655   1.3e-0004                   *
*   0.500000    0.666533    0.666667   1.3e-0004                   *
*   0.550000    0.645026    0.645161   1.4e-0004                   *
*   0.600000    0.624865    0.625000   1.4e-0004                   *
*   0.650000    0.605926    0.606061   1.3e-0004                   *
*   0.700000    0.588103    0.588235   1.3e-0004                   *
*   0.750000    0.571298    0.571429   1.3e-0004                   *
*   0.800000    0.555428    0.555556   1.3e-0004                   *
*   0.850000    0.540416    0.540541   1.2e-0004                   *
*   0.900000    0.526194    0.526316   1.2e-0004                   *
*   0.950000    0.512703    0.512821   1.2e-0004                   *
*   1.000000    0.499886    0.500000   1.1e-0004                   *
* ------------------------------------------------                 *
*                                                                  * 
*                             Pascal Release By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)          *
*******************************************************************}
Program AdamsBashforth;
Uses WinCrt;

const H=0.05;  {integration step}

var
    B,X,Y: Array[0..3] of double;
    I,K: Integer;
    C1,C2,C3,C4,ER,X1,YEX: double;

    {User defined function Y'=F(X,Y) }
    Function F(X,Y:double): double;
    Begin
      F:=-Y + X/((1+X)*(1+X))
    End;

    {Exact solution Y=FX(X) }
    Function FX(X:double): double;
    Begin
      FX:=1/(1+X)
    End;

    {Runge-Kutta method to calculate first points only}
    Procedure RK4;
    Begin
      C1:=F(X[K],Y[K]);
      C2:=F(X[K]+H/2,Y[K]+H/2*C1);
      C3:=F(X[K]+H/2,Y[K]+H/2*C2);
      C4:=F(X[K]+H,Y[K]+H*C3);
      X[K+1]:=X[K]+H;
      Y[K+1]:=Y[K]+H*(C1+2*C2+2*C3+C4)/6
    End;

{main program}
BEGIN
  {Initial conditions}
  X[0]:=0.0;  {starting X} 
  X1:=1.0;    {ending X}
  Y[0]:=1.0;  {initial Y}
  {write header}
  writeln('       X           Y        Y exact      Error    ');
  writeln('  ------------------------------------------------');
  {write initial line}
  writeln(X[0]:12:6,Y[0]:12:6,Y[0]:12:6,'   0');
  {use Runge-Kutta to start}
  For K:=0 to 1 do
  begin
    RK4; YEX:=FX(X[K+1]); ER:=ABS(YEX-Y[K+1]);
    writeln(X[K+1]:12:6,Y[K+1]:12:6,YEX:12:6,'  ',ER:-1)
  end;
  {main integration loop}
  Repeat
    For I:=1 to 3 do B[I]:=F(X[3-I],Y[3-I]);
    X[3]:=X[2]+H;
    Y[3]:=Y[2]+H*(23*B[1]-16*B[2]+5*B[3])/12;        
    YEX:=FX(X[3]); ER:=ABS(Y[3]-YEX);
    writeln(X[3]:12:6,Y[3]:12:6,YEX:12:6,'  ',ER:-1);
    For K:=0 to 2 do
    begin
      X[K]:=X[K+1]; Y[K]:=Y[K+1]
    end    
  Until X[2]>=X1;
  writeln('  ------------------------------------------------');

  ReadKey;
  DoneWinCrt

END.

{end of file adambash.pas}