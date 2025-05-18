{************************************************************
*         Determinant of a complex square matrix            *
*          By Gauss Method with full pivoting               *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
* Calculate the determinant of complex matrix:              *
*  ( 47,-15) ( 62,5) (  0,-72) (61, 20)                     *
*  (  6, 14) (-17,3) (-102,91) ( 7,-12)                     *
*  ( 13, 55) ( 32,8) (  41, 7) (25,  1)                     *
*  (111,25)  ( 40,0) ( 12,-82) (58,-30)                     *
*                                                           *
*  Det =  1.74165640000000E+0007 - 1.05983200000305E+0007 I *
*                                                           *
* --------------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal       *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988".         *
************************************************************} 
Program CDetMat;

Uses WinCrt;

Const NMAX = 20;

  Type
       Complex = Record
         R, I: Real
       End;

       Matc = Array[1..NMAX,1..NMAX] of Complex;
       Veci = Array[1..NMAX] of Integer;

  Var
       n: Integer;
       eps: Real;

       A: Matc; det:Complex;


  Function CABS(Z:Complex): Real;
  Begin
    CABS := sqrt(Z.R*Z.R+Z.I*Z.I)
  End;

  Procedure CADD(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R + Z2.R;
    Z.I := Z1.I + Z2.I
  End;

  Procedure CDIF(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R - Z2.R;
    Z.I := Z1.I - Z2.I
  End;

  Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
  Begin
    Z.R := Z1.R*Z2.R - Z1.I*Z2.I;
    Z.I := Z1.R*Z2.I + Z1.I*Z2.R
  End;

  Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
  Var d:Real; C:Complex;
  Begin
    d := Z2.R*Z2.R + Z2.I*Z2.I;
    if d<1E-10 then
      writeln(' Complex Divide by zero !')
    else
    begin
      C.R:=Z2.R; C.I:=-Z2.I;
      CMUL(Z1,C,Z);
      Z.R:=Z.R/d; Z.I:=Z.I/d
    end
  End;

  {****************************************************************
  * TSCGT procedure implements the triangularization algorithm of *
  * Gauss with full pivoting at each step for a complex matrix, A *
  * and saves the made transformations in KP and LP.              *
  * ------------------------------------------------------------- *
  * INPUTS:                                                       *
  *          N:   size of complex matrix A                        *
  *          A:   complex matrix of size N x N                    *
  * OUTPUTS;                                                      *
  *          it:  =0 if A is singular, else =1.                   *
  *           C:  contains the upper triangular matrix and the    *
  *               multipliers used during the process.            *
  *          KP:  contains the column exchanges.                  *
  *          LP:  contains the line exchanges.                    *
  ****************************************************************}
  Procedure TSCGT(eps:Real; N:integer; A:Matc; Var it:integer; Var C:Matc;
                    Var KP,LP:Veci);
  Var i,j,k,k0,l0:integer;
      C0,C1,P0,T0:Complex;
  Begin
    C:=A;
    it:=1; K:=1;
    While (it=1) and (k<N) do
    begin
      P0:=C[k,k]; l0:=k; k0:=k;
      For i:=k to N do
        For j:=1 to N do
          if CABS(C[i,j]) > CABS(P0) then
          begin
            P0:=C[i,j];
            l0:=i; k0:=j
          end;
      LP[k]:=l0; KP[k]:=k0;
      if CABS(P0) < eps then
        it:=0
      else
      begin
        if l0<>k then
          For j:=k to N do
          begin
            T0:=C[k,j];
            C[k,j]:=C[l0,j];
            C[l0,j]:=T0
          end;
        if k0<>k then
          For i:=1 to N do
          begin
            T0:=C[i,k];
            C[i,k]:=C[i,k0];
            C[i,k0]:=T0
          end;
        For i:=k+1 to N do
        begin
          C0:=C[i,k]; CDIV(C0,P0,C[i,k]);
          For j:=k+1 to N do
          begin
            C0:=C[i,j];
            CMUL(C[i,k],C[k,j],C1);
            CDIF(C0,C1,C[i,j])
          end
        end;
        Inc(k)
      end
    end;
    if (it=1) and (CABS(C[N,N]) < eps) then it:=0
  End; {TSCGT}

  {***************************************************************
  * The DCGT procedure calculates the complex determinant of a   *
  * complex square matrix by the Gauss method with full pivoting *
  * ------------------------------------------------------------ *
  * INPUTS:                                                      *
  *         eps:  required precision                             *
  *          N :  size of matrix A                               *
  *          A :  complex matrix of size N x N                   *
  * OUTPUT:                                                      *
  *         det:  complex determinant.                           *                  
  ***************************************************************}
  Procedure DCGT(eps:Real; N:integer; A:Matc; Var det:Complex);
  Var it,k,l: integer;
      C0,Z0,Z1: Complex;
      C:Matc;
      KP,LP:Veci;
  Begin
    Z0.R:=0.0; Z0.I:=0.0;
    Z1.R:=1.0; Z1.I:=0.0;
    TSCGT(eps,N,A,it,C,KP,LP);
    if it=0 then
      det:=Z0
    else
    begin
      det:=Z1;
      For k:=1 to N do
      begin
        C0:=det; CMUL(C0,C[k,k],det)
      end;
      l:=0;
      For K:=1 to N-1 do
      begin
        if LP[k]<>k then Inc(l);
        if KP[k]<>k then Inc(l)
      end;
      if Odd(l) then
      begin
        C0:=det; det.R:=-C0.R; det.I:=-C0.I
      end
    end
  end; {DCGT}

{main program}
BEGIN

(* Example #1
  N:=3;   {size of complex matrix A}

  A[1,1].R:=1.0; A[1,1].I:=0.0;
  A[1,2].R:=0.0; A[1,2].I:=1.0;
  A[1,3].R:=0.0; A[1,3].I:=1.0;

  A[2,1].R:=0.0; A[2,1].I:=1.0;
  A[2,2].R:=1.0; A[2,2].I:=0.0;
  A[2,3].R:=1.0; A[2,3].I:=0.0;

  A[3,1].R:= 0.0; A[3,1].I:=1.0;
  A[3,2].R:= 1.0; A[3,2].I:=0.0;
  A[3,3].R:=-1.0; A[3,3].I:=0.0;

  ( Det = -4.0 + 0 I )   *)

{ Example #2 }
  N:=4;

  A[1,1].R:=47.0; A[1,1].I:=-15.0;
  A[1,2].R:=62.0; A[1,2].I:=  5.0;
  A[1,3].R:= 0.0; A[1,3].I:=-72.0;
  A[1,4].R:=61.0; A[1,4].I:= 20.0;

  A[2,1].R:=   6.0; A[2,1].I:= 14.0;
  A[2,2].R:=- 17.0; A[2,2].I:=  3.0;
  A[2,3].R:=-102.0; A[2,3].I:= 91.0;
  A[2,4].R:=   7.0; A[2,4].I:=-12.0;

  A[3,1].R:= 13.0; A[3,1].I:=-55.0;
  A[3,2].R:= 32.0; A[3,2].I:=  8.0;
  A[3,3].R:= 41.0; A[3,3].I:=  7.0;
  A[3,4].R:= 25.0; A[3,4].I:=  1.0;

  A[4,1].R:=111.0; A[4,1].I:= 25.0;
  A[4,2].R:= 40.0; A[4,2].I:=  0.0;
  A[4,3].R:= 12.0; A[4,3].I:=-82.0;
  A[4,4].R:= 58.0; A[4,4].I:=-30.0;

  eps:=1E-10;

  DCGT(eps,N,A,det);

  writeln;
  write(' Det = ',det.R);
  if det.I>=0.0 then write(' + ')
                else write(' - ');
  writeln(abs(det.I),' I');
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file cdetmat.pas}