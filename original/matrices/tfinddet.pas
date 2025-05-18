{************************************************************
 !* This program calculates the determinant of a real square *
 !* matrix using Function FindDet (see ref. below).          *
 !* -------------------------------------------------------- *
 !* SAMPLE RUN:                                              *
 !*                                                          *
 !* Calculate the deteminant of square matrix:               *
 !*                                                          *
 !*           ! 0  -1  -1  -1 !                              *
 !*           ! 1   0  -1  -1 !                              *
 !*           ! 1   1   0  -1 !                              *
 !*           ! 1   1   1   0 !                              *
 !*                                                          *
 !*  Determinant =  1.000000E+00                             *
 !*                                                          *
 !* -------------------------------------------------------- *
 !*              TPW Version By Jean-Pierre Moreau, Paris.   *
 !***********************************************************} 
Program Tfinddet;

Uses WinCrt;

Type MAT = Array[1..4,1..4] of double;

Var

  N: Integer;
  A: MAT;
  D: double;


 {----------------------------------------------------------------------------------------------------
 !Function to find the determinant of a square matrix
 !Author : Louisda16th a.k.a Ashwith J. Rego
 !Description: The subroutine is based on two key points:
 !1] A determinant is unaltered when row operations are performed: Hence, using this principle,
 !row operations (column operations would work as well) are used
 !to convert the matrix into upper traingular form
 !2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
 !----------------------------------------------------------------------------------------------------}
 Function FindDet(matrix: MAT; n: integer): Double;
 Label 10;
 Var m, tmp: double;
     i, j, k, l: integer;
     DetExists: Boolean;
 Begin
     DetExists := TRUE;
     l := 1;
     {Convert to upper triangular form}
     for k := 1 to n-1 do
     begin
         if (matrix[k,k] = 0) then
         begin
             DetExists := FALSE;
             for i := k+1 to n do
                 if (matrix[i,k] <> 0) then
                 begin
                     for j := 1 to n do
                     begin
                         tmp := matrix[i,j];
                         matrix[i,j] := matrix[k,j];
                         matrix[k,j] := tmp
                     end;
                     DetExists := TRUE;
                     l:=-l;
                     goto 10
                 end;
10:          if (DetExists = FALSE) then
                 FindDet := 0.0
         end;
         for j := k+1 to n do
         begin
             m := matrix[j,k]/matrix[k,k];
             for i := k+1 to n do 
                 matrix[j,i] := matrix[j,i] - m*matrix[k,i]
         end
     end; {of k loop}                        
     
     {Calculate determinant by finding product of diagonal elements}
     tmp := l;
     for i := 1 to n do
         tmp := tmp * matrix[i,i];

     FindDet := tmp
     
 End; {FindDet}

 
BEGIN

 N:=4;

 A[1][1] :=  0; A[1][2] :=  1; A[1][3] :=  1; A[1][4] := 1;
 A[2][1] := -1; A[2][2] :=  0; A[2][3] :=  1; A[2][4] := 1;
 A[3][1] := -1; A[3][2] := -1; A[3][3] :=  0; A[3][4] := 1;
 A[4][1] := -1; A[4][2] := -1; A[4][3] := -1; A[4][4] := 0;

 D := FindDet(A, N);

 Writeln;
 Writeln(' Determinant = ', D);

 ReadKey;
 DoneWinCrt

END.

{end of file TFindDet.pas}