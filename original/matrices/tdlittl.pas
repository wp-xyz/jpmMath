{*************************************************************
*       SOLVE A LINEAR SYSTEM BY DIRECT FACTORIZATION        *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Solve linear system A X = B, where:                       *
*                                                            *
*       1  0  0  0  0  1              1                      *
*       1  1  0  0  0 -1              0                      *
*  A = -1  1  1  0  0  1          B = 1                      *
*       1 -1  1  1  0 -1              0                      *
*      -1  1 -1  1  1  1              1                      *
*       1 -1  1 -1  1 -1              0 )                    *
*                                                            *
* LINEAR SYSTEM AX = B:                                      *
*                                                            *
*  1.0000  0.0000  0.0000  0.0000  0.0000  1.0000    1.0000  *
*  1.0000  1.0000  0.0000  0.0000  0.0000 -1.0000    0.0000  *
* -1.0000  1.0000  1.0000  0.0000  0.0000  1.0000    1.0000  *
*  1.0000 -1.0000  1.0000  1.0000  0.0000 -1.0000    0.0000  *
* -1.0000  1.0000 -1.0000  1.0000  1.0000  1.0000    1.0000  *
*  1.0000 -1.0000  1.0000 -1.0000  1.0000 -1.0000    0.0000  *
*                                                            *
* SOLUTION:                                                  *
*  3.43750000000000E-0001                                    *
*  3.12500000000000E-0001                                    *
*  3.75000000000000E-0001                                    *
*  2.50000000000000E-0001                                    *
*  5.00000000000000E-0001                                    *
*  6.56250000000000E-0001                                    *
*                                                            *
* ---------------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77 *
*       [BIBLI 18].                                          *
*                                                            *
*                           TPW Release By J-P Moreau, Paris *
*                                  (www.jpmoreau.fr)         *
*************************************************************}
PROGRAM TEST_DLITTL;

Uses WinCrt;

Const
       NMAX=50;

Type
       pMAT = ^MAT;
        MAT = Array[1..NMAX,1..NMAX] of Double;
       pVEC = ^VEC;
        VEC = Array[1..NMAX] of Double;

Var
       A, W: pMAT;
       B, X, Z: pVEC;

       I,J,N: Integer;


    Procedure DLITTL(A:pMAT; B:pVEC;N:Integer;Var X:pVEC; W:pMAT;Z:pVEC);
{-----------------------------------------------------------------------
!     SOLVE A LINEAR SYSTEM BY DIRECT FACTORIZATION
!     A*X = L*U*X = B
!     L(LOWER TRIANGULAR MATRIX WITH L(I,I):=1)
!     U(UPPER TRIANGULAR MATRIX),THE SUB-DIAGONAL TERMS OF L AND THE
!     TERMS U(I,J) ARE STORED IN THE SAME MATRIX, W. THE RESOLUTION
!     IS OBTAINED BY SOLVING TWO TRIANGULAR SYSTEMS:
!     L*Z = B  AND U*X = Z.
!     THE PARTIAL PIVOTING IS MADE BY CHOOSING THE BIGGEST ELEMENT OF
!     EACH COLUMN OF THE TRANSFORMATION MATRIX.
!     INPUTS:
!     A = TABLE OF SIZE (N,N)                                      R*8
!     B = SECOND MEMBER VECTOR (N)                                 R*8
!     N = ORDER OF LINEAR SYSTEM                                   I*4
!     OUTPUTS:
!     X = SOLUTION VECTOR (N)                                      R*8
!     WORKING ZONE:
!     W = TABLE OF SIZE (LDA,N)                                    R*8
!     Z = AUXILIARY VECTOR (N)                                     R*8
!     NOTE:
!     MESSAGE '** DLITTL ** NO UNIQUE SOLUTION' IS GIVEN WHEN A 
!     NULL PIVOT IS FOUND.
!-----------------------------------------------------------------------}
    Label 2,10,32;
    Var
        SUM,AMAX,T: Double;
        I,J,K,P: Integer;
    Begin
      P:=1;
      AMAX:=ABS(A^[1,1]);
{     SEEK MAX. ELEMENT OF FIRST COLUMN }
      For J:=2 to N do
      begin
        IF ABS(A^[J,1]) < AMAX THEN GOTO 2;
        AMAX:=A^[J,1];
        P:=J;
2:    end;
      IF AMAX = 0.0 THEN GOTO 32;
      IF P <> 1 THEN
{     EXCHANGE LINES OF MATRIX A AND ORDER OF UNKNOWNS
      LINKED TO B  }
      begin
        For J:=1 to N do
        begin
          T:=A^[1,J];
          A^[1,J]:=A^[P,J];
          A^[P,J]:=T
        end;
        T:=B^[1];
        B^[1]:=B^[P];
        B^[P]:=T
      end;
{     FIRST LINE OF U AND FIRST COLUMN OF L }
      W^[1,1]:=A^[1,1];
      For J:=2 to N do
      begin
        W^[1,J]:=A^[1,J];
        W^[J,1]:=A^[J,1]/W^[1,1]
      end;
{     SECOND LINE OF U AND SECOND COLUMN OF L
      (N-1)TH LINE OF U AND (N-1)TH COLUMN OF L }
      For I:=2 to N-1 do
      begin
        For J:=I to N do
        begin
{     SEEK PIVOT }
          AMAX:=0.0;
          SUM:=0.0;
          For K:=1 to I-1 do
            SUM:=SUM+W^[J,K]*W^[K,I];
          T:=ABS(A^[J,I]-SUM);
          IF T < AMAX THEN GOTO 10;
          AMAX:=T;
          P:=J;
10:     end;
        IF AMAX = 0.0 THEN GOTO 32;
        IF P <> I THEN
{     EXCHANGE LINES OF MATRICES A , W AND ORDER OF UNKNOWNS
      LINKED TO B  }
        begin
          For J:=1 to N do
          begin
            T:=A^[I,J];
            A^[I,J]:=A^[P,J];
            A^[P,J]:=T;
            T:=W^[I,J];
            W^[I,J]:=W^[P,J];
            W^[P,J]:=T
          end;
          T:=B^[I];
          B^[I]:=B^[P];
          B^[P]:=T
        end;
{     CALCULATE THE U(I,I) }
        SUM:=0.0;
        For K:=1 to I-1 do
          SUM:=SUM+W^[I,K]*W^[K,I];
        W^[I,I]:=A^[I,I]-SUM;
        IF W^[I,I] = 0.0 THEN GOTO 32;
        For J:=I+1 to N do
        begin
          SUM:=0.0;
          For K:=1 to I-1 do
            SUM:=SUM+W^[I,K]*W^[K,J];
{     CALCULATE THE U(I,J) }
          W^[I,J]:=A^[I,J]-SUM;
          SUM:=0.0;
          For K:=1 to I-1 do
            SUM:=SUM+W^[J,K]*W^[K,I];
{     CALCULATE THE L(I,J) }
          W^[J,I]:=(A^[J,I]-SUM)/W^[I,I]
        end
      end;
      SUM:=0.0;
      For K:=1 to N-1 do
        SUM:=SUM+W^[N,K]*W^[K,N];
{     CALCULATE  U(N,N) }
      W^[N,N]:=A^[N,N]-SUM;
      IF W^[N,N] = 0.0 THEN GOTO 32;
{     SOLVE SYSTEM  Z = L(-1)*B }
      Z^[1]:=B^[1];
      For I:=2 to N do
      begin
        SUM:=0.0;
        For K:=1 to I-1 do
          SUM:=SUM+W^[I,K]*Z^[K];
        Z^[I]:=B^[I]-SUM
      end;
{     SOLVE SYSTEM  X = U(-1)*Z }
      X^[N]:=Z^[N]/W^[N,N];
      For I:=N-1 Downto 1 do
      begin
        SUM:=0.0;
        For K:=I+1 to N do
          SUM:=SUM+W^[I,K]*X^[K];
        X^[I]:=(Z^[I]-SUM)/W^[I,I]
      end;
{     END OF RESOLUTION OF  A*X = L*(U*X) = B }
      exit;
32:   WRITELN(' ** DLITTL ** NO UNIQUE SOLUTION.')
    End;


{main program}
BEGIN

  N:=6;     {size of system}

  {allocate memory}
  New(A); New(W); New(B);  New(X); New(Z);

  For I:=1 to N do
  begin
    For J:=1 to N do A^[I,J]:=0.0;
    B^[I]:=0.0
  end;

  {define elements of A <> 0}
  A^[1,1]:= 1.0; A^[1,6]:= 1.0;
  A^[2,1]:= 1.0; A^[2,2]:= 1.0; A^[2,6]:=-1.0;
  A^[3,1]:=-1.0; A^[3,2]:= 1.0; A^[3,3]:= 1.0; A^[3,6]:= 1.0;
  A^[4,1]:= 1.0; A^[4,2]:=-1.0; A^[4,3]:= 1.0; A^[4,4]:= 1.0; A^[4,6]:=-1.0;
  A^[5,1]:=-1.0; A^[5,2]:= 1.0; A^[5,3]:=-1.0; A^[5,4]:= 1.0; A^[5,5]:= 1.0; A^[5,6]:= 1.0;
  A^[6,1]:= 1.0; A^[6,2]:=-1.0; A^[6,3]:= 1.0; A^[6,4]:=-1.0; A^[6,5]:= 1.0; A^[6,6]:=-1.0;

  {define elements of B <> 0}
  B^[1]:=1.0; B^[3]:=1.0; B^[5]:=1.0;

  {print data}
  writeln;
  writeln(' LINEAR SYSTEM AX=B:');
  writeln;
  For I:=1 to N do
  begin
    For J:=1 to N do write(A^[I,J]:8:4);
    writeln('   ', B^[I]:8:4)
  end;

  {call main procedure}  
  DLITTL(A,B,N,X,W,Z);

  {print results}
  writeln;
  writeln(' SOLUTION:');
  For I:=1 to N do writeln(' ',X^[I]);

  {free memory}
  Dispose(A); Dispose(W); Dispose(B);  Dispose(X); Dispose(Z);

  ReadKey; DoneWinCrt
  
END.


{ End of file tdlittl.pas}