{***************************************************
*         Test program for Horner's rule           *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Input the five coefficients:                     *
*                                                  *
*    A( 0) =  1                                    *
*    A( 1) =  4                                    *
*    A( 2) =  6                                    *
*    A( 3) =  4                                    *
*    A( 4) =  1                                    *
*                                                  *
* What is the expansion point: 1                   *
*                                                  *
* The shifted coefficients are:                    *
*                                                  *
*    B( 0) =  16                                   *
*    B( 1) =  32                                   *
*    B( 2) =  24                                   *
*    B( 3) =  8                                    *
*    B( 4) =  1                                    *
*                                                  *
***************************************************}
PROGRAM Horner;
Uses Wincrt;

Var
     C : Array[0..4,0..5] of DOUBLE;
     A, B :  Array[0..10] of DOUBLE;
     x0 : DOUBLE;
     i  : integer;


{**************************************************
*      Horner's shifting rule subroutine          *
* ----------------------------------------------- *
* This routine takes a given quartic polynomial   *
* and converts it to a Taylor expansion.          *
* The input series coefficients are A(i), the     *
* expansion point is x0. The shifted coefficients *
* are returned in B(i).                           *
**************************************************}
PROCEDURE Horner_Shift;
Label 100, 200;
Var i, j : integer;
Begin
     for j:=0 to 4 do
       C[j,0] := A[4-j];
     for i:=0 to 4 do
     begin
       C[0,i+1] := C[0,i];
       j := 1;
100:   if j > 4-i then goto 200;
       C[j,i+1] := x0*C[j-1,i+1]+C[j,i];
       j := j+1;
       goto 100;
200: end;
     for i:=0 to 4 do
       B[4-i] := C[i,4-i+1]
End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the five coefficients:');
  writeln;
  for i := 0 to 4 do
  begin
    write('  A(',i,') = '); read(A[i])
  end;
  writeln;
  write(' What is the expansion point: '); read(x0);
  writeln;

  Horner_Shift;

  writeln(' The shifted coefficients are:');
  writeln;
  for i := 0 to 4 do
    writeln('  B(',i,') = ',B[i]:3:0);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file horner.pas}