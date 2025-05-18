{***************************************************
*   Program to demonstrate Legendre coefficients   *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Legendre polynomial coefficients for order 2     *
*                                                  *  
* A( 0) =   -0.500000                              *
* A( 1) =    0.000000                              *
* A( 2) =    1.500000                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 3     *
*                                                  *  
* A( 0) =    0.000000                              *
* A( 1) =   -1.500000                              *
* A( 2) =    0.000000                              *
* A( 3) =    2.500000                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 4     *
*                                                  *  
* A( 0) =    0.375000                              *
* A( 1) =    0.000000                              *
* A( 2) =   -3.750000                              *
* A( 3) =    0.000000                              *
* A( 4) =    4.375000                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 5     *
*                                                  *  
* A( 0) =    0.000000                              *
* A( 1) =    1.875000                              *
* A( 2) =    0.000000                              *
* A( 3) =   -8.750000                              *
* A( 4) =    0.000000                              *
* A( 5) =    7.875000                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 6     *
*                                                  *  
* A( 0) =   -0.312500                              *
* A( 1) =    0.000000                              *
* A( 2) =    6.562500                              *
* A( 3) =    0.000000                              *
* A( 4) =  -19.687500                              *
* A( 5) =    0.000000                              *
* A( 6) =   14.437500                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 7     *
*                                                  *  
* A( 0) =    0.000000                              *
* A( 1) =   -2.187500                              *
* A( 2) =    0.000000                              *
* A( 3) =   19.687500                              *
* A( 4) =    0.000000                              *
* A( 5) =  -43.312500                              *
* A( 6) =    0.000000                              *
* A( 7) =   26.812500                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 8     *
*                                                  *  
* A( 0) =    0.273438                              *
* A( 1) =    0.000000                              *
* A( 2) =   -9.843750                              *
* A( 3) =    0.000000                              *
* A( 4) =   54.140625                              *
* A( 5) =    0.000000                              *
* A( 6) =  -93.843750                              *
* A( 7) =    0.000000                              *
* A( 8) =   50.273438                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 9     *
*                                                  *  
* A( 0) =    0.000000                              *
* A( 1) =    2.460938                              *
* A( 2) =    0.000000                              *
* A( 3) =  -36.093750                              *
* A( 4) =    0.000000                              *
* A( 5) =  140.765625                              *
* A( 6) =    0.000000                              *
* A( 7) = -201.093750                              *
* A( 8) =    0.000000                              *
* A( 9) =   94.960938                              *
*                                                  *  
* <Enter> to continue...                           *
*                                                  *  
* Legendre polynomial coefficients for order 10    * 
*                                                  *  
* A( 0) =   -0.246094                              *
* A( 1) =    0.000000                              *
* A( 2) =   13.535156                              *
* A( 3) =    0.000000                              *
* A( 4) = -117.304688                              *
* A( 5) =    0.000000                              *
* A( 6) =  351.914063                              *
* A( 7) =    0.000000                              *
* A( 8) = -427.324219                              *
* A( 9) =    0.000000                              *
* A(10) =  180.425781                              *
*                                                  * 
***************************************************}
PROGRAM Legendre;
Uses WinCrt;

VAR
        A : Array[0..10] of DOUBLE;
        B : Array[0..10,0..10] of DOUBLE;
        k, n : INTEGER;

{*****************************************************
* Legendre series coefficients evaluation subroutine *
* by means of recursion relation. The order of the   *
* polynomial is n. The coefficients are returned in  *
* A(i).                                              *
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


{main program}
BEGIN
  clrscr;
  for n:=2 to 10 do
  begin
    writeln;
    writeln(' Legendre polynomial coefficients for order ',n);
    writeln;

    Legendre_Coeff;

    for k:=0 to n do
      writeln('   A(',k:2,') = ',A[k]:11:6);
    if n <> 10 then readkey;
  end;
  ReadKey; DoneWinCrt
END.

{End of file legendre.pas}

