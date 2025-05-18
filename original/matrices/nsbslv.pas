{*************************************************************
*   Solve banded linear system AX = B By LU decomposition    *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Solve banded linear system AX = B, where:                 *
*           1 10  0  0  0             1                      *
*           9  2 20  0  0             1                      *
*       A = 0 19  3 30  0         B = 1                      *
*           0  0 29  4 40             1                      *
*           0  0  0 39  5             1  )                   *
*                                                            *
* Note: A is given as follows:                               *
*           0  0  0  0  0                                    *
*           0 10 20 30 40                                    *
*           1  2  3  4  5                                    *
*           9 19 29 39  0                                    *
*           0  0  0  0  0                                    *
*                                                            *
*  SOLUTION:                                                 *
*                                                            *
*   0.40664649                                               *
*   0.05933535                                               *
*  -0.13892446                                               *
*   0.00964672                                               *
*   0.12475556                                               *
*                                                            *
*  Error code: 0                                             *
*                                                            *
* ---------------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77 *
*       [BIBLI 18].                                          *
*                                                            *
*                     TPW Release 1.0 By J-P Moreau, Paris.  *
*                              (www.jpmoreau.fr)             *
*************************************************************}
PROGRAM TEST_NSBSLV;

Uses WinCrt;

Const  NMAX = 25;

Type
       pMAT = ^MAT;
       MAT = Array[1..NMAX,1..NMAX] of Double;
       pVEC = ^VEC;
       VEC = Array[1..NMAX] of Double;
       pIVEC = ^IVEC;
       IVEC = Array[1..NMAX] of Integer;

Var
       A: pMAT;      {pointer to banded matrix}
       B: pVEC;      {pointer to second member vector}
       X: pVEC;      {pointer to solution vector}
       IPVT:pIVEC;   {pointer to auxiliary integer vector}
       IND:Integer;  {error code}

       I,J,N,MU,ML: Integer;


       Function Max(a,b:Integer): Integer;
       Begin
         if a>=b then Max:=a
         else Max:=b
       End;

       Function Min(a,b:Integer): Integer;
       Begin
         if a<=b then Min:=a
         else Min:=b
       End;

       FUNCTION IAMAX(A:VEC; N:Integer): Integer;
       Var T:Double; I:Integer;
       Begin
         T:=0.0;
         For I:=1 to N do
           IF ABS(A[I]) > T THEN
           begin
             T:=ABS(A[I]);
             IAMAX:=I
           end
       End;

       Procedure SCALE(N:Integer;T:Double; Var A:VEC);
       Var I: Integer;
       Begin
         For I:=1 to N do A[I]:=T*A[I]
       End;

       Procedure DAXPY(N:Integer; A:Double; X:VEC; Var Y:VEC);
       Var I: Integer;
       Begin
         For I:=1 to N do Y[I]:=Y[I] + A*X[I]
       End;


    Procedure NSBFAC(Var B:pMAT; N,ML,MU: Integer; Var IPVT: pIVEC; Var IND:Integer);
{-------------------------------------------------------------------
!     LU factorization of a band matrix (non symmetric) with partial
!     pivoting.
!     INPUTS:
!     B  : banded matrix. The correspondance between full matrix
!     A(i,j) and band matrix B(k,l) is given by following sequence:
!     m:=ml+mu+1;
!     for j:=1 to n do
!     begin
!       i1:=max(1,j-mu);
!       i2:=min(n,j+ml);
!       for i:=i1 to i2 do
!       begin
!         k:=i-j+m;
!         B[k,j]:=A[i,j]
!       end
!     end;
!     N   : size of B
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     OUTPUTS:
!     B   : banded matrix storing the LU elements, the ML first lines
!           of which are used for pivoting.
!     IPVT: integer vector of size N storing the pivoting indices.
!     IND : flag = 0,  B is non singular
!                = k,  B may be singular
!------------------------------------------------------------------}
    Label  10, 20;
    Var T: Double;
        I0,I,J,J0,J1,JU,JZ,K,KP1,L,LM,M,MM,NM1: Integer;
        TMP,TMP1: VEC;
    Begin
      M:=ML+MU+1;
      IND:=0;

      J0:=MU+2;
      J1:=MIN(N,M)-1;
      IF J1 >= J0 THEN
        For JZ:=J0 to J1 do
        begin
          I0:=M+1-JZ;
          For I:=I0 to ML do B^[I,JZ]:=0.0
        end;
      JZ:=J1;
      JU:=0;

      NM1:=N-1;
      IF NM1 >= 1 THEN
        For K:=1 to NM1 do
        begin
          KP1:=K+1;
          JZ:=JZ+1;
          IF (JZ <= N) AND (ML >= 1) THEN
            For I:=1 to ML do B^[I,JZ]:=0.0;

          LM:=MIN(ML,N-K);
          For I:=1 to LM+1 do TMP[I]:=B^[I+M-1,K];
          L:=IAMAX(TMP,LM+1)+M-1;
          IPVT^[K]:=L+K-M;
          IF B^[L,K] = 0.0 THEN GOTO 10;
          IF L <> M THEN
          begin
            T:=B^[L,K];
            B^[L,K]:=B^[M,K];
            B^[M,K]:=T
          end;
          T:=-1.0/B^[M,K];

          For I:=1 to LM do TMP[I]:=B^[I+M,K];
          SCALE(LM,T,TMP);
          For I:=1 to LM do B^[I+M,K]:=TMP[I];

          JU:=MIN(MAX(JU,MU+IPVT^[K]),N);
          MM:=M;
          IF JU >= KP1 THEN
            For J:=KP1 to JU do
            begin
              L:=L-1;
              MM:=MM-1;
              T:=B^[L,J];
              IF L <> MM THEN
              begin
                B^[L,J]:=B^[MM,J];
                B^[MM,J]:=T
              end;
              For I:=1 to LM do
              begin
                TMP[I]:=B^[I+M,K];
                TMP1[I]:=B^[I+MM,J]
              end;
              DAXPY(LM,T,TMP,TMP1);
              For I:=1 to LM do
                B^[I+MM,J]:=TMP1[I]
            end;
          GOTO 20;
10:       ind:=K;
20:     end;
      IPVT^[N]:=N;
      IF B^[M,N] = 0.0 THEN IND:=N
    End;


    Procedure NSBSLV(A:pMAT; N,ML,MU: Integer; IPVT:pIVEC; B:pVEC; Var X:pVEC);
{--------------------------------------------------------------------
!     Solve banded linear system AX = B
!     INPUTS:
!     A   : banded matrix as output of LU factorization by NSBFAC
!           (see storing mode in NSBFAC subroutine).
!     LDA : 1st dimension of A in calling program (lda >= 2ml+mu+1)
!     N   : order of A                             ---------------
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     IPVT: integer vector of size N storing the pivoting indices
!           as output of NSBFAC.
!     B   : second member vector
!     OUTPUT:
!     X   : solution vector
!--------------------------------------------------------------------}
    Var T: Double;
        I,K,L,LA,LB,LM,M,NM1: Integer;
        TMP,TMP1:VEC;
    Begin
      For I:=1 to N do X^[I]:=B^[I];
      M:=ML+MU+1;
      NM1:=N-1;
{     solve L*y = B }
      IF (ML <> 0) AND (NM1 >= 1) THEN
        For K:=1 to NM1 do
        begin
          LM:=MIN(ML,N-K);
          L:=IPVT^[K];
          T:=X^[L];
          IF L <> K THEN
          begin
            X^[L]:=X^[K];
            X^[K]:=T
          end;
          For I:=1 to LM do
          begin
            TMP[I]:=A^[I+M,K];
            TMP1[I]:=X^[I+K]
          end;
          DAXPY(LM,T,TMP,TMP1);
          For I:=1 to LM do X^[I+K]:=TMP1[I]
        end;
{     solve U*y = X }
      For K:=N Downto 1 do
      begin
        X^[K]:=X^[K]/A^[M,K];
        LM:=MIN(K,M)-1;
        LA:=M-LM;
        LB:=K-LM;
        T:=-X^[K];
        For I:=1 to LM do
        begin
          TMP[I]:=A^[I+LA-1,K];
          TMP1[I]:=X^[I+LB-1]
        end;
        DAXPY(LM,T,TMP,TMP1);
        For I:=1 to LM do X^[I+LB-1]:=TMP1[I]
      end
    End;


{main program}
BEGIN

  New(A); New(B); New(X); New(IPVT);

  N:=5;            {size of system}

  For I:=1 to N do
    For J:=1 to N do
      A^[I,J] := 0.0;

  MU:=1;           {one upper subdiagonal}
  ML:=1;           {one lower subdiagonal}

{ Define banded A matrix }
  A^[2,2]:=10.0; A^[2,3]:=20.0; A^[2,4]:=30.0; A^[2,5]:=40.0;
  A^[3,1]:= 1.0; A^[3,2]:= 2.0; A^[3,3]:= 3.0; A^[3,4]:= 4.0; A^[3,5]:=5.0;
  A^[4,1]:= 9.0; A^[4,2]:=19.0; A^[4,3]:=29.0; A^[4,4]:=39.0;

{ Example #1 (B(1)..B(5) := ONE) }
  For I:=1 to N do B^[I] := 1.0;

{ Example #2 (results are X1..X5 := ONE)
  B^[1]:=11.0; B^[2]:=31.0; B^[3]:=52.0; B^[4]:=73.0; B^[5]:=44.0; } 

{ call LU factorization }  
  NSBFAC(A,N,ML,MU,IPVT,IND);

{ call appropriate solver }
  NSBSLV(A,N,ML,MU,IPVT,B,X);

{ print solution and error code (must be zero) }
  writeln;
  writeln(' SOLUTION:');
  writeln;
  For I:=1 to N do writeln(X^[I]:12:8);
  writeln;
  writeln(' Error code: ', IND);
  writeln;

  ReadKey;
  Dispose(A); Dispose(B); Dispose(X); Dispose(IPVT);
  DoneWinCrt

END.

{ end of file nsbslv.pas}