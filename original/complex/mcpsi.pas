{********************************************************************
*         Calculate the Psi function for a complex argument         *
* ----------------------------------------------------------------- *
* EXPLANATION:                                                      *
*                                                                   *
*      Purpose: This program computes the psi function psi(z)       *
*               for a complex argument using subroutine CPSI        *
*      Input :  x   --- Real part of z                              *
*               y   --- Imaginary part of z                         *
*      Output:  PSR --- Real part of psi(z)                         *
*               PSI --- Imaginary part of psi(z)                    *
*      Examples:                                                    *
*                  x       y      Re[psi(z)]     Im[psi(z)]         *
*                -------------------------------------------        *
*                 3.0     2.0     1.16459152      .67080728         *
*                 3.0    -2.0     1.16459152     -.67080728         *
*                -3.0     2.0     1.39536075     2.62465344         *
*                -3.0    -2.0     1.39536075    -2.62465344         *
*                                                                   *
* NOTE:          Psi(z) = Gamma'(z) / Gamma(z)                      *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*                                                                   *
*  Please enter x and y (z=x+iy): 3 2                               *
*                                                                   *
*     x       y      Re[Psi(z)]      Im[Psi(z)]                     *
*   ----------------------------------------------                  *
*    3.0     2.0     1.16459152      0.67080728                     *
*                                                                   *
* ----------------------------------------------------------------- *
* REFERENCE:"Fortran Routines for Computation of Special Functions, *
*            jin.ece.uiuc.edu/routines/routines.html".              *
*                                                                   *
*                              Pascal Release By J-P Moreau, Paris. *
*                                       (www.jpmoreau.fr)           *
********************************************************************}  
  PROGRAM MCPSI;
  Uses WinCrt;

  Var  aa,E,X,Y,Z,PSR,PSI: Double;
       A, W: Array[1..9] of DOUBLE;

  {used by CPSI}
  Function IPower(x:double; n:integer): double;
  {Calculate x^n}
  var i,m : integer;
  result :double;
  begin
    result := 1.0;
    if n=0 then
    begin
      IPower:=result;
      exit
    end;
    m:=  n;
    if n<0 then m:=-n;
    for i:=1 to m do result :=x*result;
    IPower :=result;
    if n<0 then IPower:=1.0/result
  end;

{*** Begin utility functions for TANH ***}
Procedure WCoeffs(X:DOUBLE);
Var I:Integer;
Begin
  aa:=0.5;
  Z:=X;
  for I:=1 to 9 do
  begin
    W[I]:=0.0;
    if Z>aa then W[I]:=1.0;
    Z:=Z-W[I]*aa;
    aa:=aa/2.0
  end
End;

Procedure Coeffs;
Begin
  E:=Exp(1.0);
  A[1]:=1.648721270700128;
  A[2]:=1.284025416687742;
  A[3]:=1.133148453066826;
  A[4]:=1.064494458917859;
  A[5]:=1.031743407499103;
  A[6]:=1.015747708586686;
  A[7]:=1.007843097206448;
  A[8]:=1.003913889338348;
  A[9]:=1.001955033591003
End;

{Modified cordic exponential subroutine}
Function XExp(X:DOUBLE): DOUBLE;
{ This subroutine takes an input value and returns Y:=EXP(X)
  X may be any positive or negative real value. }
Label fin;
Var Y:DOUBLE;
    I,K:Integer;
Begin
{ Get coefficients }
  Coeffs;
{ Reduce the range of X }
  K:=Round(INT(X));
  X:=X-K;
{ Determine the weighting coeffs, W(I) }
  WCoeffs(X);
{ Calculate products }
  Y:=1.0;
  for I:=1 to 9 do
    if W[I]>0.0 then Y:=Y*A[I];
{ Perform residual multiplication }
  Y:=Y*(1.0+Z*(1.0+Z/2.0*(1.0+Z/3.0*(1.0+Z/4.0))));
{ Account for factor EXP(K) }
  if K<0 then E:=1.0/E;
  if ABS(K)<1 then goto fin;
  for I:=1 to ABS(Round(K)) do Y:=Y*E;
{ Restore X }
  X:=X+K;
fin:XExp:=Y
End;

{---------------------------------------------*
*          Hyperbolic sine Function           *
* ------------------------------------------- *
* This Procedure uses the definition of the   *
* hyperbolic sine and the modified cordic     *
* exponential Function XExp(X) to approximate *
* SINH(X) over the entire range of real X.    *
* ------------------------------------------- }
Function SinH(X:DOUBLE): DOUBLE;
Label 10,fin;
Var Y:DOUBLE;
    I:Integer;
Begin
{ Is X small enough to cause round off error? }
  if ABS(X)<0.35 then goto 10;
{ Calculate SINH(X) using exponential definition
  Get Y:=EXP(X) }
  Y:=XExp(X);
  Y:=(Y-(1.0/Y))/2.0;
  goto fin;
10: {series approximation (for X small) }
  Z:=1.0; Y:=1.0;
  for I:=1 to 8 do
  begin
    Z:=Z*X*X/((2*I)*(2*I+1));
    Y:=Y+Z
  end;
  Y:=X*Y;
fin:SinH:=Y
End;

{ hyperbolic cosine Function }
Function CosH(X:DOUBLE): DOUBLE;
Var Y:DOUBLE;
Begin
  Y:=XExp(X);
  Y:=(Y+(1.0/Y))/2.0;
  CosH:=Y
End;

{*** End utility functions for TANH ***}

{ hyperbolic tangent Function
  TANH(X]:=SINH(X)/COSH(X)   }
Function TanH(X:DOUBLE): DOUBLE;
Var V,Y: DOUBLE;
Begin
  V:=SinH(X);
  Y:=CosH(X);
  TanH:=V/Y
End;


      Procedure CPSI(X,Y:Double; Var PSR,PSI:Double);
{     =============================================
        Purpose: Compute the psi function for a
                 complex argument
        Input :  x   --- Real part of z
                 y   --- Imaginary part of z
        Output:  PSR --- Real part of psi(z)
                 PSI --- Imaginary part of psi(z)
      ============================================= }
      Var A: array[1..8] of Double;
          CT2,RI,RR,TH,TM,TN,X0,X1,Y1,Z0,Z2:Double;
          K,N:Integer;

      Begin

        {Initialize Table A}
        A[1]:=-0.8333333333333E-01;     A[2]:=0.83333333333333333E-02;
        A[3]:=-0.39682539682539683E-02; A[4]:=0.41666666666666667E-02;
        A[5]:=-0.75757575757575758E-02; A[6]:=0.21092796092796093E-01;
        A[7]:=-0.83333333333333333E-01; A[8]:=0.4432598039215686;
        
        IF (Y = 0.0) AND (X = INT(X)) AND (X <= 0.0) THEN
        begin
           PSR:=1.0E+200;
           PSI:=0.0
        end
        ELSE
        begin
           IF X < 0.0 THEN
           begin
             X1:=X;
             Y1:=Y;
             X:=-X;
             Y:=-Y
           end;
           X0:=X;
           IF X < 8.0 THEN
           begin
             N:=8-Round(INT(X));
             X0:=X+1.0*N
           end;
           IF (X0 = 0.0) AND (Y <> 0.0) Then TH:=0.5*PI;
           IF X0 <> 0.0 Then TH:=ArcTan(Y/X0);
           Z2:=X0*X0+Y*Y;
           Z0:=SQRT(Z2);
           PSR:=Ln(Z0)-0.5*X0/Z2;
           PSI:=TH+0.5*Y/Z2;
           For K:=1 to 8 do
           begin
              PSR:=PSR+A[K]*IPower(Z2,-K)*COS(2.0*K*TH);
              PSI:=PSI-A[K]*IPower(Z2,-K)*SIN(2.0*K*TH)
           end;
           IF X < 8.0 THEN
           begin
              RR:=0.0;
              RI:=0.0;
              For K:=1 to N do
              begin
                 RR:=RR+(X0-K)/(Sqr(X0-K)+Y*Y);
                 RI:=RI+Y/(Sqr(X0-K)+Y*Y)
              end;
              PSR:=PSR-RR;
              PSI:=PSI+RI
           end;
           IF X1 < 0.0 THEN
           begin
              TN:=SIN(PI*X)/COS(PI*X);
              TM:=TANH(PI*Y);
              CT2:=TN*TN+TM*TM;
              PSR:=PSR+X/(X*X+Y*Y)+PI*(TN-TN*TM*TM)/CT2;
              PSI:=PSI-Y/(X*X+Y*Y)-PI*TM*(1.0+TN*TN)/CT2;
              X:=X1;
              Y:=Y1
           end
        end
      End; {CPSI}


  {main program}
  BEGIN
    Writeln;
    Write(' Please enter x and y (z=x+iy): '); Readln(X, Y);
    Writeln;
    Writeln('   x       y      Re[Psi(z)]      Im[Psi(z)]');
    Writeln(' ----------------------------------------------');

    CPSI(X,Y,PSR,PSI);

    Writeln(X:5:1,'   ',Y:5:1,'  ',PSR:13:8,'   ',PSI:13:8);
    Writeln;

    ReadKey;
    DoneWinCrt

  END.

{end of file mcpsi.pas}