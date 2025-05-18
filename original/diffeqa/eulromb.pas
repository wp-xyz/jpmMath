{*******************************************************************
*    Solve Y' = F(X,Y) with Initial Condition Y(X0)=Y0 using the   *
*    Euler-Romberg Method                                          *
* ---------------------------------------------------------------- *
* REFERENCE: "Méthode de calcul numérique- Tome 2 - Programmes en  *
*             Basic et en Pascal By Claude Nowakowski, Edition du  *
*             P.S.I., 1984" [BIBLI 04].                            * 
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
* (Solve Y' = X*X*Y with Y(0)=1 for X0=0 up to X1=1.0, exact       *
*  solution is Y = Exp(X^3/3) ).                                   *
*                                                                  *
* Starting X...: 0                                                 *
* Ending X.....: 1                                                 *
* Value Y at X0: 1                                                 *
* Initial step : 0.1                                               *
* Desired error: 1e-6                                              *
*                                                                  *
*    X         Y            Y       Error       Number of          *
*          estimated      exact                subdivisions        *
* ----------------------------------------------------------       *
*   0.1     1.000333    1.000333  0.00000000       4               *
*   0.2     1.002670    1.002670  0.00000001       4               *
*   0.3     1.009041    1.009041  0.00000006       4               *
*   0.4     1.021562    1.021562  0.00000014       4               *
*   0.5     1.042547    1.042547  0.00000027       4               *
*   0.6     1.074654    1.074655  0.00000086       4               *
*   0.7     1.121125    1.121126  0.00000107       4               *
*   0.8     1.186094    1.186095  0.00000126       4               *
*   0.9     1.275067    1.275069  0.00000133       4               *
*   1.0     1.395611    1.395612  0.00000114       4               *
*   1.1     1.558410    1.558411  0.00000047       4               *
* ----------------------------------------------------------       *
*                                                                  * 
*                             Pascal Release By J-P Moreau, Paris. *
*                                     (www.jpmoreau.fr)            *
*******************************************************************}
Program EulerRomberg;
Uses WinCrt;

Const
      Dim = 100;
      LA=10;       {maximum number of subdivisions}

Type  Vector = Array[0..Dim] of Double;
      IVec = Array[0..Dim] of Integer;

Var   X,Y: Vector;
      NL: IVec;
      N,NC:Integer;
      EF,ER,H,X0,X1,Y0,YEX:Double;

{Y' = F(X,Y) }
Function F(X,Y:Double): Double;
Begin
  F:=X*X*Y
End;

{Exact solution FX(X) }
Function FX(X:Double): Double;
Begin
  FX:=Exp(X*X*X/3)
End;

Procedure Euler_Romberg(NC:Integer; H,ER,X0,Y0:Double; Var X,Y:Vector; Var NL:IVEC);
{**************************************************************
* Solve an Ordinary Differential Equation Y'=F(X,Y) using the *
* Euler-Romberg Method.                                       *
* ----------------------------------------------------------- *
* Inputs:                                                     *
*          NC: number of points to calculate                  *
*           H: abscissa integration step                      *
*          ER: desired precision                              *
*          X0: starting abscissa                              *
*          Y0: initial value of Y at X=X0                     *
*           F: external user defined function Y'              *
* Outputs:                                                    *
*           X: Table of NC abscissas                          *
*           Y: Table of NC solution ordinates                 *
*          NL: Table of number of steps per abscissa          *
*                                                             *
* Note:    the ending abscissa X1 is given by X0 + NC*H.      *
**************************************************************}                                                    
Var N:Integer;
    ET,XC,YC:Double;
    T:Array[0..20] of Double;
    J,K,L,LM,M,MM:Integer;
Begin
  {Initial conditions}
  X[0]:=X0; Y[0]:=Y0;
  {main integration loop}
  For N:=0 to NC do
  begin    
    XC:=X[N]; YC:=Y[N];
    T[1]:=Y[N] + H*F(XC,YC);
    L:=1; LM:=2;
    Repeat
      XC:=X[N]; YC:=Y[N];
      For J:=1 to LM do
      begin
        XC:=XC+H/LM; YC:=YC+H/LM*F(XC,YC)
      end;
      T[L+1]:=YC; M:=1; K:=L; MM:=2; ET:=1.0;
      If K>1 Then
        Repeat
          T[K]:=(MM*T[K+1]-T[K])/(MM-1);
          ET:=Abs(T[K]-T[K-1]);
          Inc(M); Dec(K); MM:=MM*2
        Until (ET<ER) or (K=1);
      If K=1 Then
      begin
        Inc(L); LM:=LM*2 
      end
    Until (L=LA) or (ET<ER);
    X[N+1]:=X[N]+H; Y[N+1]:=T[K]; NL[N+1]:=L
  end
End;

{main program}
BEGIN
  {Initial data}
  writeln;
  write(' Starting X...: '); Readln(X0);
  write(' Ending X.....: '); Readln(X1);
  write(' Value Y at X0: '); Readln(Y0);
  write(' Initial step : '); Readln(H);
  write(' Desired error: '); Readln(ER);

  {number of calculations}
  NC:=Round((X1-X0)/H);

  {call integration routine}
  Euler_Romberg(NC,H,ER,X0,Y0,X,Y,NL);

  {write header}
  writeln;
  writeln('   X         Y            Y       Error      Number of   ');
  writeln('         estimated      exact               subdivisions ');
  writeln('---------------------------------------------------------');
  {results loop}
  For N:=0 to NC do
  begin    
    YEX:=FX(X[N+1]);
    EF:=Abs(YEX-Y[N+1]);
    Writeln(X[N+1]:6:2,Y[N+1]:12:6,YEX:12:6,EF:12:8,NL[N+1]:8)
  end;
  writeln('---------------------------------------------------------');
  ReadKey;
  DoneWinCrt
END.

{end of file eulromb.pas}