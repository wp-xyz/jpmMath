{*******************************************************
*    Collection of Chebyshev approximation routines    *
* ---------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific *
*             Computing By W.H. Press, B.P. Flannery,  *
*             S.A. Teukolsky and W.T. Vetterling,      *
*             Cambridge University Press, 1986"        *
*             [BIBLI 08].                              *
*                                                      *
*                TPW Release 1.1 By J-P Moreau, Paris. *
*                         (www.jpmoreau.fr)            *
* ---------------------------------------------------- *
* Release 1.1: added procedures CHINT and CHDER.       *
*******************************************************}
UNIT Chebyshe;

INTERFACE

Uses WinCrt;

Const NMAX = 50;
      HALF = 0.5;
      ONE  = 1.0;
      QUART= 0.25;
      TWO  = 2.0;
      ZERO = 0.0;

Type  Tab = Array[1..NMAX] of Double;

      Procedure CHEBFT(A,B: Double; Var C:Tab; N:Integer);
      Function  CHEBEV(A,B:Double;C:Tab;M:Integer;X:Double): Double;
      Procedure CHINT(A,B:Double; C:Tab; Var CINT:Tab; N:Integer);
      Procedure CHDER(A,B:Double; C:Tab; Var CDER:Tab; N:Integer);
      Function  FUNC(x:Double): Double;


IMPLEMENTATION

{user defined function}
Function FUNC(x:Double): Double;
Begin
  FUNC:=SIN(x)
End;

Procedure CHEBFT(A,B: Double; Var C:Tab; N:Integer);
{*******************************************************
* Chebyshev fit: Given a real function FUNC(X), lower  *
* and upper limits of the interval [A,B] for X, and a  *
* maximum degree N, this routine computes the N Cheby- *
* shev coefficients Ck, such that FUNC(X) is approxima-*
* ted by:  N                                           *
*         [Sum Ck Tk-1(Y)] - C1/2, where X and Y are   *
*         k=1                                          *
* related by:     Y = (X - 1/2(A+B)) / (1/2(B-A))      *
* This routine is to be used with moderately large N   *
* (e.g. 30 or 50), the array of C's subsequently to be *
* truncated at the smaller value m such that Cm+1 and  *
* subsequent elements are negligible.                  *
*******************************************************}
Var SUM:Double; F:Tab;
    BMA,BPA,FAC, Y: Double;
    J,K: Integer;
Begin
  BMA:=HALF*(B-A); BPA:=HALF*(B+A);
  for K:=1 to N do
  begin
    Y:=COS(PI*(K-HALF)/N);
    F[K]:=FUNC(Y*BMA+BPA)
  end;
  FAC:=TWO/N;
  for J:=1 to N do
  begin
    SUM:=ZERO;
    for K:=1 to N do
      SUM:=SUM+F[K]*COS((PI*(J-1))*((K-HALF)/N));
    C[J]:=FAC*SUM
  end
End;

Function CHEBEV(A,B:Double;C:Tab;M:Integer;X:Double): Double;
{*********************************************************
* Chebyshev evaluation: All arguments are input. C is an *
* array of Chebyshev coefficients, of length M, the first*
* M elements of Coutput from subroutine CHEBFT (which    *
* must have been called with the same A and B). The Che- *
* byshev polynomial is evaluated at a point Y determined *
* from X, A and B, and the result FUNC(X) is returned as *
* the function value.                                    *
*********************************************************}
Var D,DD,SV,Y,Y2: Double;
    J:Integer;
Begin
  if (X-A)*(X-B) > ZERO then writeln(' X not in range.');
  D:=ZERO; DD:=ZERO;
  Y:=(TWO*X-A-B)/(B-A);  {change of variable}
  Y2:=TWO*Y;
  for J:=M Downto 2 do
  begin
    SV:=D;
    D:=Y2*D-DD+C[J];
    DD:=SV
  end;
  CHEBEV:=Y*D-DD+HALF*C[1]
End;

Procedure CHINT(A,B:Double; C:Tab; Var CINT:Tab; N:Integer);
{*********************************************************
* Given A,B,C, as output from routine CHEBFT, and given  *
* N, the desired degree of approximation (length of C to *
* be used), this routine returns the array CINT, the Che-*
* byshev coefficients of the integral of the function    *
* whose coefficients are C. The constant of integration  *
* is set so that the integral vanishes at A.             *
*********************************************************}
Var
    CON,FAC,SUM: Double;
    J: Integer;
Begin
  CON:=QUART*(B-A);
  SUM:=ZERO;
  FAC:=ONE;
  for J:=2 to N-1 do
  begin
    CINT[J]:=CON*(C[J-1]-C[J+1])/(J-1);
    SUM:=SUM+FAC*CINT[J];
    FAC:=-FAC
  end;
  CINT[N]:=CON*C[N-1]/(N-1);
  SUM:=SUM+FAC*CINT[N];
  CINT[1]:=TWO*SUM      {set the constant of integration}
End;

Procedure CHDER(A,B:Double; C:Tab; Var CDER:Tab; N:Integer);
{*********************************************************
* Given A,B,C, as output from routine CHEBFT, and given  *
* N, the desired degree of approximation (length of C to *
* be used), this routine returns the array CDER, the Che-*
* byshev coefficients of the derivative of the function  *
* whose coefficients are C.                              *
*********************************************************} 
Var
    CON: Double;
    J: Integer;
Begin
  CDER[N]:=ZERO;
  CDER[N-1]:=TWO*(N-1)*C[N];
  if N >= 3 then
    for J:=N-2 Downto 1 do
      CDER[J]:=CDER[J+2]+TWO*J*C[J+1];
  CON:=TWO/(B-A);
  for J:=1 to N do          {normalize to interval B - A}
    CDER[J]:=CDER[J]*CON
End;

END. {of unit

End of file Chebyshe.pas}