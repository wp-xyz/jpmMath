{**************************************************
*  Evaluate a Legendre Polynomial P(x) of Order n *
*  for argument x by using Horner's Rule          *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
*                                                 *
*  give order: 7                                  *
*                                                 *
*  give argument x: 0.5                           *
*                                                 *
*  P(x) =  2.23144531250000E-0001                 *
*                                                 *
* ----------------------------------------------- *
*               TPW Release By J-P Moreau, Paris. *
*                      (www.jpmoreau.fr)          *
**************************************************}
PROGRAM Eval_Legendre;
Uses WinCrt;

Const  NMAX = 12;

Type
       pVEC = ^VEC;
       VEC  = Array[1..NMAX] of Double;

VAR
       A : Array[0..NMAX] of Double;
       B : Array[0..NMAX,0..NMAX] of DOUBLE;
       A1, R: pVec;
       k, n : INTEGER;
       x: double;

{*****************************************************
* Legendre series coefficients evaluation subroutine *
* by means of recursion relation. The order of the   *
* polynomial is n. The coefficients are returned in  *
* A(i) by increasing powers.                         *
*****************************************************}
PROCEDURE Legendre_Coeff;
Var i,j : integer;
Begin
  {Establish p0 and p1 coefficients}
  B[0,0]:=1.0 ; B[1,0]:=0.0 ; B[1,1]:=1.0;
  {Return if order is less then 2}
  if n > 1 then
  begin 
    for i:=2 to n do
    begin
      B[i,0] := -(i-1)*B[i-2,0]/i;
      for j:=1 to i do
      begin
        {Basic recursion relation}
        B[i,j]:=(i+i-1)*B[i-1,j-1]-(i-1)*B[i-2,j];
        B[i,j]:=B[i,j]/i
      end
    end;
    for i:=0 to n do A[i]:=B[n,i]
  end
End;

Procedure HORNER (N:Integer; A:pVEC; X0:Double; K:Integer; B:Boolean; Var R:pVEC);
{---------------------------------------------------------------------
!     EVALUATE A POLYNOMIAL AND ITS DERIVATIVES BY HORNER'S METHOD
!
!     INPUTS:
!     N       ORDER OF POLYNOMIAL
!     A       VECTOR OF SIZE N+1 STORING THE COEFFICIENTS OF
!             POLYNOMIAL IN DECREASING ORDERS OF POWERS
!             I.E. P(X) = A(1)*X**N+A(2)*X**(N-1)+...+A(N)*X+A(N+1)
!     X0      GIVEN ARGUMENT
!     K       MAXIMUM ORDER OF DERIVATIVES ASKED FOR
!     B       FLAG
!             = .TRUE.  EVALUATION ONLY
!             = .FALSE. DETERMINATION OF POLYNOMIAL COEFFICIENTS
!                NEAR X0
!     R       VECTOR OF SIZE K+1 CONTAINS:
!             -  VALUES P(X0), P'(X0), P''(X0),..PK(X0),
!                IF B IS TRUE AND IF K < N
!             -  THE N+1 COEFFICIENTS OF POLYNOMIAL
!                Q(X-X0) = R(1)+R(2)*(X-X0)+...+R(N+1)*(X-X0)**(N),
!                IF B IS FALSE AND IF K = N.
!
!     REFERENCE:
!     ALGORITHM 337, COLLECTED ALGORITHMS FROM CACM, W.PANKIEWICKZ
!--------------------------------------------------------------------}
Var RR: Double;
    I,J,L,NMJ: Integer;
Begin
  RR := A^[1];
  For I := 1 to K+1 do R^[I] := RR;
  For J := 2 to N+1 do
  begin
    R^[1] := R^[1]*X0+A^[J];
    NMJ := N-J+1;
    IF NMJ > K THEN
      L := K
    ELSE
      L := NMJ;
    For I := 2 to L+1 do R^[I] := R^[I]*X0+R^[I-1]
  end;
  IF (B) THEN
  begin
    L := 1;
    For I := 2 to K+1 do
    begin
      L := L*(I-1);
      R^[I] := R^[I]*L
    end
  end
End;

{Evaluate P(x) of order n for argument x}
Function Eval(x:double): double;
Var i:Integer;
Begin
  For i:=n downto 0 do A1^[n-i+1] := A[i];
  Horner(n,A1,x,0,True,R);
  Eval:=R^[1]
End;

{main program}
BEGIN

  New(A1); New(R);

  writeln;
  write(' give order: '); readln(n);
  writeln;
  write(' give argument x: '); readln(x);
  writeln;

  Legendre_Coeff;

  writeln(' P(x) = ', Eval(x));
  writeln;

  ReadKey;
  Dispose(A1); Dispose(R);
  DoneWinCrt

END.

{end of file eval_leg.pas}