{*******************************************************************
* Solve a Complex Linear System By Gauss Method with full pivoting *
* and correction process.                                          *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*  Solve Complex Linear System AX = B with:                        *
*                                                                  *
*      ( 47,-15) ( 62,5) (  0,-72) (61, 20)         ( 629,988)     *
*  A = (  6, 14) (-17,3) (-102,91) ( 7,-12) and B = (-180,825)     *
*      ( 13, 55) ( 32,8) (  41, 7) (25,  1)         ( 877,441)     *
*      (111,25)  ( 40,0) ( 12,-82) (58,-30)         ( 734,-88)     *
*                                                                  *
* System solutions:                                                *
*  X1 = -1.00000000000000E+0000 +  3.00000000000000E+0000 I        *
*  X2 =  2.00000000000000E+0000 +  9.99999999999999E+0000 I        *
*  X3 =  7.00000000000000E+0000 -  5.00000000000000E+0000 I        *
*  X4 =  1.70000000000000E+0001 +  6.00000000000000E+0000 I        *
*                                                                  *
* ---------------------------------------------------------------- *
* Reference:  "Algebre, Algorithmes et programmes en Pascal        *
*              By Jean-Louis Jardrin, DUNOD Paris, 1988".          *
*******************************************************************}
Program CSysLin;

Uses WinCrt;

  Const NMAX = 20;

  Type
       Complex = Record
         R, I: Double
       End;

       Matc = Array[1..NMAX,1..NMAX] of Complex;
       Vecc = Array[1..NMAX] of Complex;
       Veci = Array[1..NMAX] of Integer;

  Var
       I,it,M,N: Integer;
       dta, eps: Double;

       A: Matc; B,X: Vecc;


  Function CABS(Z:Complex): Double;
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
  Var d:Double; C:Complex;
  Begin
    d := Z2.R*Z2.R+Z2.I*Z2.I;
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
  Procedure TSCGT(eps:double; N:integer; A:Matc; Var it:integer; Var C:Matc;
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

  {*************************************************************
  * Function BSCGT calculates the solution of upper triangular * 
  * system [S(n-1)] by the back substitution method and deduces*
  * from it the solution of system [S]. The call must be made  *
  * only after a call to TSCGT and the matrix of [S] must be   *
  * regular.                                                   *
  * ---------------------------------------------------------- *
  * INPUTS:                                                    *
  *         N : size of matrix C                               *
  *         C:  contains the upper triangular matrix and the   *
  *             multipliers used during the triangularization  *
  *             process (in output of TSCGT).     .            *
  *         W : contains the right-side coefficients           *
  *         KP: contains the column exchanges.                 *
  *         LP: contains the line exchanges.                   *
  * OUTPUT:                                                    *
  *         Z : system solution complex vector.                *
  *************************************************************}
  Procedure BSCGT(N:integer; C:Matc; W:Vecc; KP,LP:Veci; Var Z:Vecc);
  Var i,j,k,k0,l0: integer;
      C0,C1,S,T0,Z0: Complex;
  Begin
    Z0.R:=0.0; Z0.I:=0.0;
    For k:=1 to N-1 do
    begin
      l0:=LP[k];
      if l0<>k then
      begin
        T0:=W[k]; W[k]:=W[l0]; W[l0]:=T0
      end;
      For i:=k+1 to N do
      begin
        C0:=W[i];
        CMUL(C[i,k],W[k],C1); CDIF(C0,C1,W[i])
      end
    end;
    CDIV(W[N],C[N,N],Z[N]);
    For i:=N-1 Downto 1 do
    begin
      S:=Z0;
      For j:=i+1 to N do
      begin
        C0:=S;
        CMUL(C[i,j],Z[j],C1); CADD(C0,C1,S)
      end;
      CDIF(W[i],S,C0);
      CDIV(C0,C[i,i],Z[i])
    end;
    For k:=N-1 Downto 1 do
    begin
      K0:=KP[k];
      if k0 <> k then
      begin
        T0:=Z[k]; Z[k]:=Z[k0]; Z[k0]:=T0
      end
    end
  End; {BSCGT}

  {************************************************************
  * Solve a Complex Linear System AX = B By Gauss Method with *
  * full pivoting and a correction process.                   *
  * --------------------------------------------------------- *
  * INPUTS:                                                   *
  *         eps, dta : absolute and relative precisions       *
  *         M : maximum number of iterations                  *
  *         N : size of linear system                         *
  *         A : complex matrix                                *
  *         B : right-hand side (complex vector)              *
  * OUTPUTS:                                                  *
  *         it: flag, =-1 if no convergence, =0 if matrix A   *
  *             is singular, =1 convergence ok.               *
  *         X : contains system solution (X1,X2,...,Xn)       *
  *                                                           *
  ************************************************************}
  Procedure RSLCGTC(eps,dta:Double; M,N:integer; A:Matc; B:Vecc;
                    var it:integer; Var X:Vecc);
  Var I,J,L: integer;
      phi1,phi2:Double;
      C0,C1,S,Z0: Complex;
      C:Matc; W,Z:Vecc; KP,LP:Veci;
  Begin
    Z0.R:=0.0; Z0.I:=0.0;
    TSCGT(eps,N,A,it,C,KP,LP);

    if it=1 then
    begin

      BSCGT(N,C,B,KP,LP,X);

      it:=-1; L:=1;
      While(it=-1) and (L<=M) do
      begin
        phi1:=0.0;
        For I:=1 to N do
          if CABS(X[I]) > phi1 then phi1:=CABS(X[I]);
        For I:=1 to N do
        begin
          S:=Z0;
          For J:=1 to N do
          begin
            C0:=S;
            CMUL(A[I,J],X[J],C1);
            CADD(C0,C1,S)
          end;
          CDIF(B[I],S,W[I])
        end;

        BSCGT(N,C,W,KP,LP,Z);

        For I:=1 to N do
        begin
          C0:=X[I];
          CADD(C0,Z[I],X[I])
        end;
        phi2:=0.0;
        For I:=1 to N do
          if CABS(Z[I]) > phi2 then phi2:=CABS(Z[I]);
        if (phi2/phi1) < dta then it:=1
                             else Inc(L)
      end {while}
    end {if it=1}
  End; {RSLCGTC}


{main program}
BEGIN

  N:=4;

  A[1,1].R:= 47.0; A[1,1].I:=-15.0;
  A[1,2].R:= 62.0; A[1,2].I:=  5.0;
  A[1,3].R:=  0.0; A[1,3].I:=-72.0;
  A[1,4].R:= 61.0; A[1,4].I:= 20.0;

  A[2,1].R:=   6.0; A[2,1].I:= 14.0;
  A[2,2].R:= -17.0; A[2,2].I:=  3.0;
  A[2,3].R:=-102.0; A[2,3].I:= 91.0;
  A[2,4].R:=   7.0; A[2,4].I:=-12.0;

  A[3,1].R:=13.0; A[3,1].I:=-55.0;
  A[3,2].R:=32.0; A[3,2].I:=  8.0;
  A[3,3].R:=41.0; A[3,3].I:=  7.0;
  A[3,4].R:=25.0; A[3,4].I:=  1.0;

  A[4,1].R:=111.0; A[4,1].I:= 25.0;
  A[4,2].R:= 40.0; A[4,2].I:=  0.0;
  A[4,3].R:= 12.0; A[4,3].I:=-82.0;
  A[4,4].R:= 58.0; A[4,4].I:=-30.0;

  B[1].R:= 629.0; B[1].I:=988.0;
  B[2].R:=-180.0; B[2].I:=825.0;
  B[3].R:= 877.0; B[3].I:=441.0;
  B[4].R:= 734.0; B[4].I:=-88.0;

  eps:=1E-10; dta:=1E-10; M:=1;

  RSLCGTC(eps,dta,M,N,A,B,it,X);

  writeln;
  writeln(' System solutions:');
  For I:=1 to N do
  begin
    write('  X',I,' = ', X[I].R);
    if X[I].I >= 0.0 then write(' + ')
                     else write(' - ');
    writeln(ABS(X[I].I),' I')
  end;
  writeln;

  Readkey;
  DoneWinCrt

END.

{end of file trslcgtc.pas}