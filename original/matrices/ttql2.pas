{*******************************************************************************
*     Find Eigenvalues and Eigenvectors of a symmetric tridiagonal matrix      *
*     using QL method                                                          *
* ---------------------------------------------------------------------------- *
* SAMPLE RUNS                                                                  *
*                                                                              *
* Example #1:    (tridiagonal matrix)                                          *
*          1  2  0  0     0                                                    *
*          2  4  7  0     0                                                    *
*     A =  0  7 10  8     0                                                    *
*          0  0  8 -0.75 -9                                                    *
*          0  0  0 -9    10                                                    *
*                                                                              *
*  Output file ttql2.lst contains:                                             *
*                                                                              *
*  EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC                                 *
*  TRIDIAGONAL MATRIX BY QL METHOD.                                            *
*                                                                              *
*  N =           5                                                             *    
*                                                                              *
*  Input tridiagonal matrix:                                                   *
*   1.000000  0.000000  0.000000  0.000000  0.000000                           *
*   2.000000  4.000000  2.000000  0.000000  0.000000                           *
*   0.000000  7.000000 10.000000  7.000000  0.000000                           *
*   0.000000  0.000000  8.000000 -0.750000  8.000000                           *
*   0.000000  0.000000  0.000000 -9.000000 10.000000                           *
*                                                                              *
*  Error code: 0                                                               *
*                                                                              *
*  Eigenvalues:                                                                *
*   1    -9.15659229                                                           *
*   2    -0.78071442                                                           *
*   3     2.53046878                                                           *
*   4    12.53989066                                                           *
*   5    19.11694726                                                           *
*                                                                              *
*  Eigenvectors (in lines]:                                                    *
*  -0.056408  0.286458 -0.522285  1.000000  0.469812                           *
*   1.000000 -0.890357  0.322363  0.344649  0.287721                           *
*   1.000000  0.765234 -0.446362 -0.252815 -0.304616                           *
*   0.097161  0.560616  0.656182 -0.282210  1.000000                           *
*   0.051876  0.469920  1.000000  0.728439 -0.719095                           *
*                                                                              *
* ---------------------------------------------------------------------------- *
* Reference:   From Numath Library By Tuan Dang Trong in Fortran 77            *
*              [BIBLI 18].                                                     *
*                                                                              *
*                                       TPW Release 1.0 By J-P Moreau, Paris   *
*                                                (www.jpmoreau.fr)             *
*******************************************************************************}
PROGRAM TEST_TQL2;

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

    Procedure TQL2(NM,N:Integer; Var D,E:pVec; Var Z:pMat; Var IER:Integer);
{---------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
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


{main program}
BEGIN

    N := 5;

    New(D); New(E); New(Z);

    {define main diagonal}
    D^[1]:=1.0; D^[2]:=4.0; D^[3]:=10.0; D^[4]:=-0.75; D^[5]:=10.0;
    {define subdiagonal}
    E^[1]:=0.0; E^[2]:=2.0; E^[3]:=7.0; E^[4]:=8.0; E^[5]:=-9.0;

    {Initialize matrix Z to unity matrix}
    For I:=1 to N do
      For J:=1 to N do
       if J=I then Z^[I,J] := 1.0
       else  Z^[I,J] := 0.0;

    Assign(fp,'ttql2.lst'); Rewrite(fp);

    writeln(fp);
    writeln(fp,' EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC');
    writeln(fp,' TRIDIAGONAL MATRIX BY QL METHOD.');
    writeln(fp);;
    writeln(fp,' N = ', N);
    writeln(fp);
    writeln(fp,' Input tridiagonal matrix:');
    for I:=1 to N do
    begin
      for J:=1 to N do
      begin
	if J < I-1 then write(fp,0.0:10:6);
	if J = I-1 then write(fp,E^[I]:10:6);
        if J = I+1 then write(fp,E^[I+1]:10:6);
	if J = I then write(fp,D^[I]:10:6);
	if (J>I+1) and (J<=N) then write(fp,0.0:10:6)
      end;
      writeln(fp);
    end;
    writeln(fp);

    TQL2(N,N,D,E,Z,IERR);

    writeln(fp,' Error code: ', IERR);
    writeln(fp);
    {print eigenvalues}
    writeln(fp,' Eigenvalues:');
    for I:=1 to N do
      writeln(fp,I:3, D^[I]:15:8);
    writeln(fp);
    {print normalized eigenvectors (with respect to unity) }
    writeln(fp,' Eigenvectors (in lines]:');
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
    writeln(' Results in file ttql2.lst');
    writeln;
      
    ReadKey; DoneWinCrt
 
END.

{end of file ttql2.pas}