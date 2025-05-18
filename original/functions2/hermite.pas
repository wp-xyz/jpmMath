{***************************************************
*   Program to demonstrate Hermite coefficients    *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Hermite polynomial coefficients for order 2      *
*                                                  *
*   A( 0) =        -2                              *
*   A( 1) =         0                              *
*   A( 2) =         4                              *
*                                                  *
* Hermite polynomial coefficients for order 3      *
*                                                  *
*   A( 0) =         0                              *
*   A( 1) =       -12                              *
*   A( 2) =         0                              *
*   A( 3) =         8                              *
*                                                  *
* Hermite polynomial coefficients for order 4      *
*                                                  *
*   A( 0) =        12                              *
*   A( 1) =         0                              *
*   A( 2) =       -48                              *
*   A( 3) =         0                              *
*   A( 4) =        16                              *
*                                                  *
* Hermite polynomial coefficients for order 5      *
*                                                  *
*   A( 0) =         0                              *
*   A( 1) =       120                              *
*   A( 2) =         0                              *
*   A( 3) =      -160                              *
*   A( 4) =         0                              *
*   A( 5) =        32                              *
*                                                  *
* Hermite polynomial coefficients for order 6      *
*                                                  *
*   A( 0) =      -120                              *
*   A( 1) =         0                              *
*   A( 2) =       720                              *
*   A( 3) =         0                              *
*   A( 4) =      -480                              *
*   A( 5) =         0                              *
*   A( 6) =        64                              *
*                                                  *
* Hermite polynomial coefficients for order 7      *
*                                                  *
*   A( 0) =         0                              *
*   A( 1) =     -1680                              *
*   A( 2) =         0                              *
*   A( 3) =      3360                              *
*   A( 4) =         0                              *
*   A( 5) =     -1344                              *
*   A( 6) =         0                              *
*   A( 7) =       128                              *
*                                                  *
* Hermite polynomial coefficients for order 8      *
*                                                  *
*   A( 0) =      1680                              *
*   A( 1) =         0                              *
*   A( 2) =    -13440                              *
*   A( 3) =         0                              *
*   A( 4) =     13440                              *
*   A( 5) =         0                              *
*   A( 6) =     -3584                              *
*   A( 7) =         0                              *
*   A( 8) =       256                              *
*                                                  *
* Hermite polynomial coefficients for order 9      *
*                                                  *
*   A( 0) =         0                              *
*   A( 1) =     30240                              *
*   A( 2) =         0                              *
*   A( 3) =    -80640                              *
*   A( 4) =         0                              *
*   A( 5) =     48384                              *
*   A( 6) =         0                              *
*   A( 7) =     -9216                              *
*   A( 8) =         0                              *
*   A( 9) =       512                              *
*                                                  *
* Hermite polynomial coefficients for order 10     *
*                                                  *
*   A( 0) =    -30240                              *
*   A( 1) =         0                              *
*   A( 2) =    302400                              *
*   A( 3) =         0                              *
*   A( 4) =   -403200                              *
*   A( 5) =         0                              *
*   A( 6) =    161280                              *
*   A( 7) =         0                              *
*   A( 8) =    -23040                              *
*   A( 9) =         0                              *
*   A(10) =      1024                              *
*                                                  *
****************************************************
 Explanations
 ------------
 
 Hermite polynomials are defined over the range -inf. < x < inf. The weight- 
 ing function is w(x) = e^(-x^2): 

    inf.
    Sum   e^(-x^2) H (x) H (x) dx = 0     for n <> m         (3.9.8)
    -inf.           n     m
                                  = f(n)  for n = m

 The corresponding recursion relation is
 
    H   (x) = 2x H (x) - 2n H   (x)                          (3.9.9)
     n+1          n          n-1
 
 where H (x) = 1  and  H (x) = 2x
        0               1
 
 As with the other polynomials, a simple subroutine for evaluating the coefficients 
 be written (see program HERMITE). 
 
 Note that Hermite polynomials are either even or odd, depending on N, and the co- 
 efficients are integers. 
-----------------------------------------------------------------------------------}
PROGRAM Hermite;
Uses WinCrt;

Type
        real_ar = DOUBLE;

Var
        A : Array[0..10] of real_ar;
        B : Array[0..10,0..10] of real_ar;
        n, k : INTEGER;


{*************************************************
* Hermite polynomial coefficients evaluation by  *
* means of recursion relation. The order of the  *
* polynomial is n. The coefficients are returned *
* in A(i).                                       *
*************************************************}
PROCEDURE Hermite_Coeff;
Var i,j : integer;
Begin
  {Establish l0 and l1 coefficients}
  B[0,0]:=1.0 ; B[1,0]:=0.0 ; B[1,1]:=2.0;
  {Return if order is less than two}
  if n>1 then
  begin
    for i:=2 to n do
    begin
      B[i,0]:=-2.0*(i-1)*B[i-2,0];
      for j:=1 to i do
        {Basic recursion relation}
        B[i,j]:=2.0*B[i-1,j-1]-2.0*(i-1)*B[i-2,j];
    end;
    for i:=0 to n do A[i]:=B[n,i]
  end
End;

{main program}
BEGIN
  clrscr;
  for n:=2 to 10 do
  begin
    writeln;
    writeln(' Hermite polynomial coefficients for order ',n);
    writeln;

    Hermite_Coeff;

    for k:=0 to n do writeln('   A(',k:2,') = ',A[k]:9:0);
    if n<10 then ReadKey
  end;
  ReadKey; DoneWinCrt
END.

{End of file hermite.pas}