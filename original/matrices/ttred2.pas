{***********************************************************************
*   Find Eigenvalues and Eigenvectors of a real symmetric matrix by    *
*   using Householder reduction and  QL method.                        *
* -------------------------------------------------------------------- *
* SAMPLE RUNS:                                                         *
*                                                                      *
* Example #1:    (symmetric matrix)                                    *
*                                                                      *
*          4 -2 -1  0                                                  *
*     A = -2  4  0 -1                                                  *
*         -1  0  4 -2                                                  *
*          0 -1 -2  4                                                  *
*                                                                      *
*  Output file ttred2.lst contains:                                    *
*                                                                      *
* EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC MATRIX                   *
* BY HOUSEHOLDER REDUCTION AND QL METHOD.                              *
*                                                                      *
* N = 4                                                                *
*                                                                      *
* Input symmetric matrix:                                              *
*  4.000000 -2.000000 -1.000000  0.000000                              *
* -2.000000  4.000000  0.000000 -1.000000                              *
* -1.000000  0.000000  4.000000 -2.000000                              *
*  0.000000 -1.000000 -2.000000  4.000000                              *
*                                                                      *
* Error code: 0                                                        *
*                                                                      *
* Eigenvalues:                                                         *
*  1     1.00000000                                                    *
*  2     3.00000000                                                    *
*  3     5.00000000                                                    *
*  4     7.00000000                                                    *
*                                                                      *
* Eigenvectors (in lines):                                             *
*  1.000000  1.000000  1.000000  1.000000                              *
*  1.000000  1.000000 -1.000000 -1.000000                              *
*  1.000000 -1.000000  1.000000 -1.000000                              *
*  1.000000 -1.000000 -1.000000  1.000000                              *
*                                                                      *
*                                                                      *
* Example #2:    (symmetric matrix)                                    *
*                                                                      *
*          1  2  3 -7    12                                            *
*          2  4  7  3    -1                                            *
*     A =  3  7 10  8     4                                            *
*         -7  3  8 -0.75 -9                                            *
*         12 -1  4 -9    10                                            *
*                                                                      *
*                                                                      *
* EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC MATRIX                   *
* BY HOUSEHOLDER REDUCTION AND QL METHOD.                              *
*                                                                      *
* N = 5                                                                *
*                                                                      *
* Input symmetric matrix:                                              *
*  1.000000  2.000000  3.000000 -7.000000 12.000000                    *
*  2.000000  4.000000  7.000000  3.000000 -1.000000                    *
*  3.000000  7.000000 10.000000  8.000000  4.000000                    *
* -7.000000  3.000000  8.000000 -0.750000 -9.000000                    *
* 12.000000 -1.000000  4.000000 -9.000000 10.000000                    *
*                                                                      *
* Error code: 0                                                        *
*                                                                      *
* Eigenvalues:                                                         *
*  1   -10.48654517                                                    *
*  2    -7.77457973                                                    *
*  3     0.46334953                                                    *
*  4    18.29182060                                                    *
*  5    23.75595477                                                    *
*                                                                      *
* Eigenvectors (in lines):                                             *
*  0.467172 -0.007947 -0.507827  1.000000  0.264432                    *
*  1.000000 -0.326100  0.209620 -0.147689 -0.815422                    *
*  0.208812  1.000000 -0.478027 -0.275000 -0.216915                    *
*  0.093053  0.594911  1.000000  0.455666  0.050740                    *
*  0.705820 -0.012677  0.131486 -0.527499  1.000000                    *
*                                                                      *
* -------------------------------------------------------------------- *
* Reference:   From Numath Library By Tuan Dang Trong in Fortran 77    *
*              [BIBLI 18].                                             *
*                                                                      *
*                               TPW Release 1.1 By J-P Moreau, Paris   *
*                                        (www.jpmoreau.fr)             *
* -------------------------------------------------------------------- *
* Release 1.1:  added results of example #2.                           *
***********************************************************************}
PROGRAM TEST_TRED2;

Uses WinCrt;

Const

      NMAX = 50;

Type
      pMat = ^MAT;
      MAT = Array[1..NMAX,1..NMAX] of Double;
      pVec = ^VEC;
      VEC = Array[1..NMAX] of Double;

Var
      D, E: pVec;
      Z: pMat;
      zmax: Double;
      I,IERR,J,N: Integer;

      fp: TEXT;

      FUNCTION Sign(a,b : Double) : Double;
      BEGIN
        IF (b <0.0) THEN Sign := - Abs(a)
                     ELSE Sign :=   Abs(a)
      END;

    Procedure TQL2(N:Integer; Var D,E:pVec; Var Z:pMat; Var IER:Integer);
{---------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               TQL2(NM,N,D,E,Z,IER);
!     INPUTS:
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!--------------------------------------------------------------------------}
    Label 10,12,18,20,26,30,34,36,38;
    Var
        B,C,F,G,H,P,R,S,EPS,EPS1: Double;
        I,J,K,L,M,JM: Integer;
    Begin

      EPS:=0.0; JM:=30;
      IER := 0;
      IF N = 1 then  GOTO 38;

{     MACHINE EPSILON }

      IF EPS <> 0.0 then GOTO 12;
      EPS := 1.0;
10:   EPS := EPS/2.0;
      EPS1 := 1.0+EPS;
      IF EPS1 > 1.0 then GOTO 10;

12:   For I := 2 to N do E^[I-1] := E^[I];
      E^[N] := 0.0;
      F := 0.0;
      B := 0.0;

      For L := 1 to N do
      begin
        J := 0;
        H := EPS*(ABS(D^[L])+ABS(E^[L]));
        IF B < H then B := H;

{     SEEK SMALLEST ELEMENT OF SUBDIAGONAL }

        For M := L to N do
          IF ABS(E^[M]) <= B then GOTO 18;

18:     IF M = L then GOTO 26;

{     START ITERATION }

20:     IF J = JM then  GOTO 36;
        J := J+1;

{     SHIFT }

        G := D^[L];
        P := (D^[L+1]-G)/(2.0*E^[L]);
        R := SQRT(P*P+1.0);
        D^[L] := E^[L]/(P+SIGN(R,P));
        H := G-D^[L];
        For I := L+1 to N do D^[I] := D^[I] - H;
        F := F + H;

{     QL TRANSFORMATION }

        P := D^[M];
        C := 1.0;
        S := 0.0;
        For I := M-1 Downto L do
        begin
          G := C*E^[I];
          H := C*P;
          IF ABS(P) >= ABS(E^[I]) THEN
          begin
            C := E^[I]/P;
            R := SQRT(C*C+1.0);
            E^[I+1] := S*P*R;
            S := C/R;
            C := 1.0/R
          end
          ELSE
          begin
            C := P/E^[I];
            R := SQRT(C*C+1.0);
            E^[I+1] := S*E^[I]*R;
            S := 1.0/R;
            C := C*S
          end;
          P := C*D^[I]-S*G;
          D^[I+1] := H+S*(C*G+S*D^[I]);

{     ELEMENTS OF EIGENVECTORS }

          For K := 1 to N do
          begin
            H := Z^[K,I+1];
            Z^[K,I+1] := S*Z^[K,I]+C*H;
            Z^[K,I] := Z^[K,I]*C-S*H
          end
        end;
        E^[L] := S*P;
        D^[L] := C*P;
        IF ABS(E^[L]) > B then GOTO 20;

{     CONVERGENCE }

26:     D^[L] := D^[L] + F
      end;

{     SORT EIGENVALUES AND EIGENVECTORS
      IN ASVENDING ORDER  }

      For L := 2 to N do
      begin
        I := L-1;
        K := I;
        P := D^[I];
        For J := L to N do
        begin
          IF D^[J] >= P then GOTO 30;
          K := J;
          P := D^[J];
30:     end;
        IF K = I then  GOTO 34;
        D^[K] := D^[I];
        D^[I] := P;
        For J := 1 to N do
        begin
          P := Z^[J,I];
          Z^[J,I] := Z^[J,K];
          Z^[J,K] := P
        end;
34:   end;
      GOTO 38;

{     NO CONVERGENCE }

36:   IER := L;
38: End; {TQL2}

Procedure TRED2(Var Z:pMAT; N:Integer; Var D, E:pVEC);
{-----------------------------------------------------------------
   HOUSEHOLDER REDUCTION OF MATRIX Z TO A TRIDIAGONAL FORM.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.
 -----------------------------------------------------------------
  INPUTS:
    Z  (R*8)  TABLE (N,N) STORING THE REAL ELEMENTS OF INPUT
              SYMMETRIC MATRIX
    N   (I4)  SIZE OF Z
  OUTPUTS:
    D  (R*8)  MAIN DIAGONAL(N) OF THE TRIDIAGONAL MATRIX
    E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
    Z  (R*8)  TABLE (N,N) STORING THE ELEMENTS OF THE
              TRANSFORMATION MATRIX.
 -----------------------------------------------------------------}
Var
    j,k,L: Integer;
    F,G,H,hh,Xscale: Double;
Begin

  For i := N DownTo 2 do
  begin
    L := i - 1;
    H := 0.0;
    Xscale := 0.0;
    If L > 1 Then
    begin
       For k := 1 To L do
         Xscale := Xscale + Abs(Z^[i, k]);
       if Abs(Xscale) < 2E-16 Then
          E^[i] := Z^[i, L]
       else
       begin
          For k := 1 To L do
          begin
            Z^[i, k] := Z^[i, k] / Xscale;
            H := H + Z^[i, k] * Z^[i, k]
          end;
          F := Z^[i, L];
          If F > 0.0 Then
            G := -Sqrt(H)
          Else
            G := Sqrt(H);
          E^[i] := Xscale * G;
          H := H - F * G;
          Z^[i, L] := F - G;
          F := 0.0;
          For j := 1 To L do
          begin
            Z^[j, i] := Z^[i, j] / H;
            G := 0.0;
            For k := 1 To j do
              G := G + Z^[j, k] * Z^[i, k];
            For k := j + 1 To L do
              G := G + Z^[k, j] * Z^[i, k];
            E^[j] := G / H;
            F := F + E^[j] * Z^[i, j]
          end;
          hh := F / (H + H);
          For j := 1 To L do
          begin
            F := Z^[i, j];
            E^[j] := E^[j] - hh * F;
            G := E^[j];
            For k := 1 To j do
              Z^[j, k] := Z^[j, k] - F * E^[k] - G * Z^[i, k]
          end
       end
    end
    else
       E^[i] := Z^[i, L];
    D^[i] := H
  end; {i loop}
  D^[1] := 0.0;
  E^[1] := 0.0;

  For i := 1 To n do
  begin
    L := i - 1;
    If D^[i] <> 0.0 Then
    begin
      For j := 1 To L do
      begin
        G := 0.0;
        For k := 1 To L do
          G := G + Z^[i, k] * Z^[k, j];
        For k := 1 To L do
          Z^[k, j] := Z^[k, j] - G * Z^[k, i]
      end
    end;
    D^[i] := Z^[i, i];
    Z^[i, i] := 1.0;
    For j := 1 To L do
    begin
      Z^[j, i] := 0.0;
      Z^[i, j] := 0.0
    end
  end

End;


{main program}
BEGIN

    N := 4;

    New(D); New(E); New(Z);

    {define input symmetric matrix}
    Z^[1,1]:= 4.0; Z^[1,2]:=-2.0; Z^[1,3]:=-1.0; Z^[1,4]:= 0.0;
    Z^[2,1]:=-2.0; Z^[2,2]:= 4.0; Z^[2,3]:= 0.0; Z^[2,4]:=-1.0;
    Z^[3,1]:=-1.0; Z^[3,2]:= 0.0; Z^[3,3]:= 4.0; Z^[3,4]:=-2.0;
    Z^[4,1]:= 0.0; Z^[4,2]:=-1.0; Z^[4,3]:=-2.0; Z^[4,4]:= 4.0;

    Assign(fp,'ttred2.lst'); Rewrite(fp);

    writeln(fp);
    writeln(fp,' EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC MATRIX');
    writeln(fp,' BY HOUSEHOLDER REDUCTION AND QL METHOD.');
    writeln(fp);;
    writeln(fp,' N = ', N);
    writeln(fp);
    writeln(fp,' Input symmetric matrix:');
    for I:=1 to N do
    begin
      for J:=1 to N do write(fp,Z^[I,J]:10:6);
      writeln(fp);
    end;
    writeln(fp);

    TRED2(Z,N,D,E);
    TQL2(N,D,E,Z,IERR);

    writeln(fp,' Error code: ', IERR);
    writeln(fp);
    {print eigenvalues}
    writeln(fp,' Eigenvalues:');
    for I:=1 to N do
      writeln(fp,I:3, D^[I]:15:8);
    writeln(fp);
    {print normalized eigenvectors (with respect to unity) }
    writeln(fp,' Eigenvectors (in lines):');
    for J:=1 to N do
    begin
      zmax:=Z^[1,j];
      for I:=2 to N do
	if abs(Z^[I,J]) > abs(zmax)  then zmax:=Z^[I,J];
      for I:=1 to N do
        write(fp,Z^[I,J]/zmax:10:6);
      writeln(fp)
    end;
  
    close(fp);
    writeln;
    writeln(' Results in file ttred2.lst');
    writeln;
      
    ReadKey; DoneWinCrt
 
END.

{end of file ttred2.pas}