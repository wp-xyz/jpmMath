{*******************************************************
*      Program to demonstrate synthetic division       *
*              of polynomials subroutine               *
* ---------------------------------------------------- *
*   Reference: BASIC Scientific Subroutines, Vol. II   *
*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                      *
*                Pascal Version By J-P Moreau, Paris.  *
*                         (www.jpmoreau.fr)            *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
* ( Example: Divide (x+1)^5 by (x+1) )                 *
*                                                      *
* What is the degree of the first polynomial: 5        *
*                                                      *
* Input the first polynomial coefficients as prompted: *
*                                                      *
*   C( 0) = 1                                          *
*   C( 1) = 5                                          *
*   C( 2) = 10                                         *
*   C( 3) = 10                                         *
*   C( 4) = 5                                          *
*   C( 5) = 1                                          *
*                                                      *
* Input the second polynomial coefficients as prompted:*
*                                                      *
*   B( 0) = 1                                          *
*   B( 1) = 1                                          *
*                                                      *
* The coefficients of the resulting polynomial are:    *
*                                                      *
*   A( 0) =   1                                        *
*   A( 1) =   4                                        *
*   A( 2) =   6                                        *
*   A( 3) =   4                                        *
*   A( 4) =   1                                        *
*                                                      *
*******************************************************}
PROGRAM Rsyndiv;
Uses WinCrt;

CONST
        SIZE = 25;

VAR
        i,n1,n2 : INTEGER;
        A, B, C : Array[0..SIZE] of double;


{***********************************************
* Synthetic Division of polynomials Subroutine *
* -------------------------------------------- *
* It is assumed that polynomial coefficients   *
* are real.  Calculates A(x)=C(x)/B(x).        *
* C(x) is of order n1,  B(x) is of order n2,   *
* A(x) is at most of order n1-n2 (n1>n2).      *
***********************************************}
PROCEDURE Division;
Label 100;
Var i,j:integer;
Begin
 for i:=n1 downto n2 do
 begin
   A[i-n2]:=C[i]/B[n2];
   if i=n2 then goto 100;
   for j:=0 to n2 do
     C[i-j]:=C[i-j]-A[i-n2]*B[n2-j];
100: end
End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln('         DIVISION OF POLYNOMIALS');
  writeln;
  write(' What is the degree of the first polynomial: '); read(n1);
  writeln;
  writeln(' Input the first polynomial coefficients as prompted:');
  writeln;
  for i:=0 to n1 do
  begin
    write('   C(',i,') = '); read(C[i])
  end;
  writeln;
  write(' What is the degree of the second polynomial: '); read(n2);
  writeln;
  writeln(' Input the second polynomial coefficients as prompted:');
  writeln;
  for i:=0 to n2 do
  begin
    write('   B(',i,') = '); read(B[i])
  end;
  writeln;

  Division;

  writeln;
  writeln(' The coefficients of the resulting polynomial are:');
  writeln;
  for i:=0 to n1-n2 do
    writeln('   A(',i,') = ',A[i]:3:0);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file rsyndiv.pas}