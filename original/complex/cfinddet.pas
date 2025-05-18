{************************************************************
 !*  This program calculates the determinant of a complex    *
 !* square matrix using Function CFindDet (see ref. below).  *
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
 !*  Determinant =  (1.000000, 0.000000)                     *
 !*                                                          *
 !* -------------------------------------------------------- *
 !*                                                          *
 !*                   Pascal Version By J-P Moreau, Paris.   *
 !***********************************************************} 
Program TCfinddet;

Uses WinCrt;

Const NMAX = 25;

  Type
       Complex = Record
         R, I: Real
       End;

       Matc = Array[1..NMAX,1..NMAX] of Complex;

  Var

   N: integer;
   A: Matc;
   D: Complex;


  Function CABS(Z: COMPLEX): Real;
  Begin
    CABS := sqrt(Z.R*Z.R + Z.I*Z.I)
  End;

  Procedure CADD(Z1: COMPLEX; Z2: COMPLEX; Var Z: COMPLEX);
  Begin
    Z.R := Z1.R + Z2.R;
    Z.I := Z1.I + Z2.I;
  End;

  Procedure CDIF(Z1: COMPLEX; Z2: COMPLEX; Var Z: COMPLEX);
  Begin
    Z.R := Z1.R - Z2.R;
    Z.I := Z1.I - Z2.I;
  End;

  Procedure CMUL(Z1: COMPLEX; Z2: COMPLEX; Var Z: COMPLEX);
  Begin
    Z.R := Z1.R*Z2.R - Z1.I*Z2.I;
    Z.I := Z1.R*Z2.I + Z1.I*Z2.R;
  End;

  Procedure CDIV(Z1: COMPLEX; Z2: COMPLEX; Var Z: COMPLEX);
  Var
    d: double; C: COMPLEX;
  Begin
    d := Z2.R*Z2.R+Z2.I*Z2.I;
    if d<1E-12 then
      writeln(' Complex Divide by zero!')
    else
    begin
      C.R:=Z2.R; C.I:=-Z2.I;
      CMUL(Z1,C,Z);
      Z.R:=Z.R/d; Z.I:=Z.I/d
    end
  End;
 
 {----------------------------------------------------------------------------------------------------
 !Procedure to find the determinant of a complex square matrix
 !Author : Louisda16th a.k.a Ashwith J. Rego
 !Description: The subroutine is based on two key points:
 !1] A determinant is unaltered when row operations are performed: Hence, using this principle,
 !row operations (column operations would work as well) are used
 !to convert the matrix into upper traingular form
 !2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
 !----------------------------------------------------------------------------------------------------}
 Procedure CFindDet(matrix: Matc; n: integer; Var cdet: COMPLEX);
 Label 10, 20;
 Var
     m, temp, CONE, CZERO: Complex;
     i, j, k, l: integer;
     DetExists: Boolean;
 Begin
     DetExists := TRUE;
     l := 1;
     CZERO.R:=0.0; CZERO.I:=0.0;
     CONE.R:=1.0; CONE.I:=0.0;
     {Convert to upper triangular form}
     for k := 1 to n-1 do
     begin
         if (matrix[k,k].R = CZERO.R) and (matrix[k,k].I = CZERO.I) then
         begin
             DetExists := FALSE;
             for i := k+1 to n do
             begin
                 if (matrix[i,k].R <> CZERO.R) or (matrix[i,k].I <> CZERO.I) then
                 begin
                     for j := 1 to n do
                     begin
                         temp.R := matrix[i][j].R;
                         temp.I := matrix[i][j].I;
                         matrix[i][j].R := matrix[k][j].R;
                         matrix[i][j].I := matrix[k][j].I;
                         matrix[k][j].R := temp.R;
                         matrix[k][j].I := temp.I
                     end;
                     DetExists := TRUE;
                     l:=-l;
                     goto 20
                 end;
             end;
20:          if (DetExists = FALSE) then
             begin
                 cdet.R := CZERO.R;
                 cdet.I := CZERO.I;
                 goto 10
             end
         end;
         for j := k+1 to n do
         begin
             {m = matrix(j,k)/matrix(k,k) }
	     CDIV(matrix[j,k],matrix[k,k],m);
             for i := k+1 to n do
             begin
                 {matrix(j,i) = matrix(j,i) - m*matrix(k,i) }
                 CMUL(m,matrix[k,i],temp);
		 matrix[j][i].R := matrix[j][i].R - temp.R;
                 matrix[j][i].I := matrix[j][i].I - temp.I
             end
         end
     end;
     
     {Calculate determinant by finding product of diagonal elements}
     cdet.R := 1.0; cdet.I := 0.0;
     for i := 1 to n do
     begin
         {CDet = CDet * matrix(i,i) }
	 CMUL(cdet,matrix[i][i],temp);
	 cdet.R := temp.R;
         cdet.I := temp.I
     end;
     
10: End;


BEGIN

 N:=4;

 A[1][1].R := 0.0; A[1][1].I := 0.0;
 A[2][1].R := 1.0; A[2][1].I := 0.0;
 A[3][1].R := 1.0; A[3][1].I := 0.0;
 A[4][1].R := 1.0; A[4][1].I := 0.0;

 A[1][2].R := -1.0; A[1][2].I := 0.0;
 A[2][2].R :=  0.0; A[2][2].I := 0.0;
 A[3][2].R :=  1.0; A[3][2].I := 0.0;
 A[4][2].R :=  1.0; A[4][2].I := 0.0;

 A[1][3].R := -1.0; A[1][3].I := 0.0;
 A[2][3].R := -1.0; A[2][3].I := 0.0;
 A[3][3].R :=  0.0; A[3][3].I := 0.0;
 A[4][3].R :=  1.0; A[4][3].I := 0.0;

 A[1][4].R := -1.0; A[1][4].I := 0.0;
 A[2][4].R := -1.0; A[2][4].I := 0.0;
 A[3][4].R := -1.0; A[3][4].I := 0.0;
 A[4][4].R :=  0.0; A[4][4].I := 0.0;

 
 CFindDet(A, N, D);

 Writeln;
 Writeln(' Determinant = (', D.R, ',', D.I, ')');
 Writeln;

 ReadKey;
 DoneWinCrt

End.

{end of file cfinddet.pas}