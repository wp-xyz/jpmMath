{********************************************************************
* THIS PROGRAM CALCULATES THE EIGENVALUES/EIGENVECTORS OF A COMPLEX *
* MATRIX C = A+I*Z BY THE JACOBI METHOD, THIS METHOD IS ALSO RECOM- *
* MANDED TO DEAL WITH REAL MATRICES THE EIGENVALUES OF WHICH ARE    *
* COMPLEX.                                                          *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
* Example #1:                                                       *
* (Find all eigenvalues and eigenvectors of 5x5 matrix:             *
*           1.0   2.0   3.0  -7.0   12.0                            *
*           2.0   4.0   7.0   3.0  -10.0                            *
*           3.0   7.0  10.0   8.0    4.0                            *
*          -7.0   3.0   8.0  -0.75  -9.0                            *
*          12.0 -10.0   4.0  -9.0   10.0  )                         *
*                                                                   *
*  Matrix C (real parts):                                           *
*  1.000000  2.000000  3.000000 -7.000000 12.000000                 *
*  2.000000  4.000000  7.000000  3.000000-10.000000                 *
*  3.000000  7.000000 10.000000  8.000000  4.000000                 *
* -7.000000  3.000000  8.000000 -0.750000 -9.000000                 *
* 12.000000-10.000000  4.000000 -9.000000 10.000000                 *
*  Matrix C (imaginary parts):                                      *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*                                                                   *
*  Eigenvalues (real parts):                                        *
*  3.349223 -10.152415 -12.544478  17.027498  26.570173             *
*  Eigenvalues (imaginary parts):                                   *
*  0.000000   0.000000   0.000000   0.000000   0.000000             *
*  Error code (must be zero): 0                                     *
*                                                                   *
*  Eigenvectors in columns (real parts):                            *
*  0.671571  0.935133 -0.563486  0.360472  0.568860                 *
*  1.000000 -0.231511  0.802269  0.418719 -0.516435                 *
* -0.432515 -0.427944 -0.567937  1.000000 -0.152070                 *
* -0.622247  1.000000  0.609801  0.229197 -0.576420                 *
* -0.290043 -0.140178  1.000000  0.295367  1.000000                 *
*  Eigenvectors in columns (imaginary parts):                       *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*  0.000000  0.000000  0.000000  0.000000  0.000000                 *
*                                                                   *
* Example #2:                                                       *
*  Matrix C (real parts):                                           *
*  1.000000  2.000000  3.000000 -7.000000 12.000000                 *
*  2.000000  4.000000  7.000000  3.000000 -1.000000                 *
*  3.000000  7.000000 10.000000  8.000000  4.000000                 *
* -7.000000  3.000000  8.000000 -0.750000 -9.000000                 *
* 12.000000 -1.000000  4.000000 -9.000000 10.000000                 *
*  Matrix C (imaginary parts):                                      *
*  1.000000  1.000000  1.000000  1.000000  1.000000                 *
*  1.000000  1.000000  1.000000  1.000000  1.000000                 *
*  1.000000  1.000000  1.000000  1.000000  1.000000                 *
*  1.000000  1.000000  1.000000  1.000000  1.000000                 *
*  1.000000  1.000000  1.000000  1.000000  1.000000                 *
*  Eigenvalues (real parts):                                        *
*  0.468367 -7.775245-10.354718 18.634772 23.276824                 *
*  Eigenvalues (imaginary parts):                                   *
*  0.040316  0.003317  0.942489  3.320560  0.693318                 *
*  Error code (must be zero):           0                           *
*                                                                   *
*  Eigenvectors in columns (real parts):                            *
*  0.209070  1.000000  0.460630  0.199619  0.693336                 *
*  1.000000 -0.326405 -0.010671  0.590353 -0.101481                 *
* -0.481799  0.206172 -0.514461  1.000000  0.001810                 *
* -0.273820 -0.142932  1.000000  0.380689 -0.610896                 *
* -0.217450 -0.813120  0.265266  0.198153  1.000000                 *
*  Eigenvectors in columns (imaginary parts):                       *
*  0.001998  0.000000 -0.003799 -0.132150  0.023902                 *
*  0.000000  0.001431 -0.059816  0.069630  0.183347                 *
* -0.033383  0.016553 -0.054470  0.000000  0.276292                 *
*  0.009639 -0.022869  0.000000  0.218001  0.166455                 *
* -0.005084 -0.011084 -0.036704 -0.220804  0.000000                 *
* ----------------------------------------------------------------- *
* REFERENCE:                                                        *
*            From Numath Library By Tuan Dang Trong in Fortran 77   *
*            [BIBLI 18].                                            *
*                                                                   *
*                         Pascal Release 1.1 By J-P Moreau, Paris.  *
*                                    (www.jpmoreau.fr)              *
*                                                                   *
* Release 1.1: added example #2.                                    *
********************************************************************}   
PROGRAM TEST_COMEIG;

Uses WinCrt;

  Const NMAX=25;
	  
  Type
       pMAT = ^MAT;
        MAT = Array[1..NMAX,1..NMAX] of Double;
       pVEC = ^VEC;
        VEC = Array[1..2*NMAX] of Double;

  Var  AR, AI, T, U: pMAT;
       EN: pVEC;
       i,j, IER, N: Integer;


    FUNCTION ABSCD(AR,AI:Double): Double;
{     ABSOLUTE VALUE OF A COMPLEX NUMBER C=AR+I*AI
      ABSCD=SQRT(AR^2+AI^2)   }
    Var
      XR,XI,W: Double;
    Begin
      XR:=ABS(AR);
      XI:=ABS(AI);
      IF XR <= XI THEN
      begin
        W:=XR;
        XR:=XI;
        XI:=W
      end;
      IF XI = 0.0 THEN
        ABSCD:=XR
      ELSE
        ABSCD:=XR*SQRT(1.0+(XI/XR)*(XI/XR))
    End;

    Procedure DIVCD(AR,AI,BR,BI:Double; Var ZR,ZI:Double);
{     COMPLEX DIVISION Z=ZR+I*ZI=(AR+I*AI)/(BR+I*BI)
      DO NOT USE IF BR=BI=0  }
    Var
      YR,YI,W: Double;
    Begin
      YR:=BR;
      YI:=BI;
      IF ABS(YR) > ABS(YI) THEN
      begin
        W:=YI/YR;
        YR:=W*YI+YR;
        ZR:=(AR+W*AI)/YR;
        ZI:=(AI-W*AR)/YR
      end
      ELSE
      begin
        W:=YR/YI;
        YI:=W*YR+YI;
        ZR:=(W*AR+AI)/YI;
        ZI:=(W*AI-AR)/YI
      end
    End;

    Procedure COMEIG (NDIM,N,NN:Integer; A,Z,T,U:pMAT; Var IER:Integer; EN:pVEC);
{----------------------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE EIGENVALUES/EIGENVECTORS OF A COMPLEX MATRIX 
!     C = A+I*Z BY THE JACOBI METHOD, THIS METHOD IS ALSO RECOMMANDED TO DEAL WITH
!     REAL MATRICES THE EIGENVALUES OF WHICH ARE COMPLEX.

!     DATA:
!     NDIM    1ST DIMENSION OF TABLES A, Z, T, U IN MAIN PROGRAM (HERE NDIM=N)
!     N       REAL SIZE OF COMPLEX MATRIX C = A+I*Z
!     NN      MAXIMUM NUMBER OF ITERATIONS
!     A       TABLE STORING THE REAL PART OF GIVEN MATRIX C
!     Z       TABLE STORING THE IMAGINARY PART OF GIVEN MATRIX C

!     OUTPUTS:
!     A(J,J),Z(J,J),J:=1,N   IN MAIN DIAGONALS OF TABLES A AND Z, YOU
!             HAVE NOW RESPECTIVELY THE REAL AND IMAGINARY PARTS OF EIGENVALUES.
!     T,U     THESE TABLES CONTAIN NOW RESPECTIVELY THE REAL AND IMAGINARY PARTS
!             OF THE EIGENVECTORS MATRIX X = T+I*U (STORED IN COLUMNS).
!     IER     ERROR CODE
!             = 0  CONVERGENCE OK
!             = 1  NO CONVERGENCE AFTER NN ITERATIONS

!     WORKING ZONE:
!     EN      TABLE OF SIZE 2*N

!     NOTES:
!     1/      IN CASE OF CONVERGENCE (IER = 0), THE MATRIX EQUATION  C*X = LAMDA*X
!             IS VERIFIED TO THE MACHINE PRECISION.
!
!     REFERENCE:
!     P.J.EBERLEIN, NUMER.MATH 14, PP 232-245 (1970)
!----------------------------------------------------------------------------------}
    Label 10,20, 65,70, 90,95,97, Return;
    Var
        MARK: Boolean;
        B,DN,E,EMAX,EPS,EPS1,TAU,W1,W2: Double;
        BR,BI,DR,DI,ER,EI,G,HI,HJ,HR,T1,T2: Double;
        C,CA,CX,D,DE,DX,ROOT,ROOT1,ROOT2,S,SA,SB,SH,SIG,SX,SW: Double;
        CB,CH,CN,COS2A,COTX,COT2X,ETA,SIN2A,TANH,TI,TR,TSE: Double;
        C1I,C2I,C1R,C2R,S1I,S2I,S1R,S2R,VI,VR,ZI,ZM: Double;
        AKI,AMI,ZKI,ZMI,AIK,AIM,ZIK,ZIM,TIK,TIM,UIK,UIM: Double;
        I,IM,IT,J,K,M:Integer;
    Begin

{     CHECK MACHINE EPSILON (AROUND 1.2E-16 FOR PC) }
      EPS := 1.0;
10:   EPS := 0.5*EPS;
      EPS1 := EPS+1.0;
      IF EPS1 > 1.0 Then Goto 10;
      MARK := FALSE;

{     INITIALIZE EIGENVECTORS }

      For I := 1 to N do
      begin
        T^[I,I] := 1.0;
        U^[I,I] := 0.0;
        For J := I+1 to N do
        begin
          T^[I,J] := 0.0;
          T^[J,I] := 0.0;
          U^[I,J] := 0.0;
          U^[J,I] := 0.0
        end
      end;
      IT := 0;
20:   IT := IT+1;

{     SAFETY TEST IN CASE OF NO CONVERGENCE }

      IF IT > NN Then Goto 90;
      IF MARK Then Goto 95;

{     DEFINE CONVERGENCE CRITERIUM }

      TAU := 0.0;
      For K := 1 to N do
      begin
        W1 := 0.0;
        For I := 1 to N do
          IF I <> K Then W1 := W1+ABS(A^[I,K])+ABS(Z^[I,K]);
        TAU := TAU+W1;
        EN^[K] := W1+ABS(A^[K,K])+ABS(Z^[K,K])
      end;

{     PERMUTE  LINES AND COLUMNS }

      For K := 1 to N-1 do
      Begin
        EMAX := EN^[K];
        I := K;
        For J := K+1 to N do
          IF EN^[J] > EMAX THEN
          begin
            EMAX := EN^[J];
            I := J
          end;
        IF I <> K THEN
        begin
          EN^[I] := EN^[K];
          For J := 1 to N do
          begin
            W2 := A^[K,J];
            A^[K,J] := A^[I,J];
            A^[I,J] := W2;
            W2 := Z^[K,J];
            Z^[K,J] := Z^[I,J];
            Z^[I,J] := W2
          end;
          For J := 1 to N do
          begin
            W2 := A^[J,K];
            A^[J,K] := A^[J,I];
            A^[J,I] := W2;
            W2 := Z^[J,K];
            Z^[J,K] := Z^[J,I];
            Z^[J,I] := W2;
            W2 := T^[J,K];
            T^[J,K] := T^[J,I];
            T^[J,I] := W2;
            W2 := U^[J,K];
            U^[J,K] := U^[J,I];
            U^[J,I] := W2
          end
        end
      End;

{     CONVERGENCE IF TAU < 100*EPS }

      IF TAU < 100.0*EPS Then Goto 95;

{     BEGIN ITERATIONS }

      MARK := TRUE;
      For K := 1 to N-1 do
      begin
        For M := K+1 to N do
        begin
          G := 0.0;
          HR := 0.0;
          HJ := 0.0;
          HI := 0.0;
          For I := 1 to N do
            IF (I <> K) AND (I <> M) THEN
            begin
              HR := HR+A^[K,I]*A^[M,I]+Z^[K,I]*Z^[M,I]-A^[I,K]*A^[I,M]-Z^[I,K]*Z^[I,M];
              HI := HI+Z^[K,I]*A^[M,I]-A^[K,I]*Z^[M,I]-A^[I,K]*Z^[I,M]+Z^[I,K]*A^[I,M];
              T1 := A^[I,K]*A^[I,K]+Z^[I,K]*Z^[I,K]+A^[M,I]*A^[M,I]+Z^[M,I]*Z^[M,I];
              T2 := A^[I,M]*A^[I,M]+Z^[I,M]*Z^[I,M]+A^[K,I]*A^[K,I]+Z^[K,I]*Z^[K,I];
              G := G+T1+T2;
              HJ := HJ-T1+T2
            end;
          BR := A^[K,M]+A^[M,K];
          BI := Z^[K,M]+Z^[M,K];
          ER := A^[K,M]-A^[M,K];
          EI := Z^[K,M]-Z^[M,K];
          DR := A^[K,K]-A^[M,M];
          DI := Z^[K,K]-Z^[M,M];
          T1 := BR*BR+EI*EI+DR*DR;
          T2 := BI*BI+ER*ER+DI*DI;

          IF T1 >= T2 THEN
          begin
            SW := 1.0;
            C := BR;
            S := EI;
            D := DR;
            DE := DI;
            ROOT2 := SQRT(T1)
          end
          ELSE
          begin
            SW :=-1.0;
            C := BI;
            S :=-ER;
            D := DI;
            DE := DR;
            ROOT2 := SQRT(T2)
          end;
          ROOT1 := SQRT(S*S+C*C);
          SIG := 1.0;
          IF D < 0.0 Then SIG := -1.0;
          CA := 1.0;
          IF C < 0.0 Then CA := -1.0;
          SA := 0.0;

          IF ROOT1 < EPS THEN
          begin
            SX := 0.0;
            SA := 0.0;
            CX := 1.0;
            CA := 1.0;
            IF SW > 0.0 THEN
            begin
              E := ER;
              B := BI
            end
            ELSE
            begin
              E := EI;
              B :=-BR
            end;
            DN := D*D+DE*DE;
            Goto 65
          end;

          IF ABS(S) > EPS THEN
          begin
            CA := C/ROOT1;
            SA := S/ROOT1
          end;
          COT2X := D/ROOT1;
          COTX := COT2X+SIG*SQRT(1.0+COT2X*COT2X);
          SX := SIG/SQRT(1.0+COTX*COTX);
          CX := SX*COTX;
          ETA := (ER*BR+BI*EI)/ROOT1;
          TSE := (BR*BI-ER*EI)/ROOT1;
          T1 := SIG*(TSE*D-ROOT1*DE)/ROOT2;
          T2 := (D*DE+ROOT1*TSE)/ROOT2;
          DN := ROOT2*ROOT2+T2*T2;
          T2 := HJ*CX*SX;
          COS2A := CA*CA-SA*SA;
          SIN2A := 2.0*CA*SA;
          W1 := HR*COS2A+HI*SIN2A;
          W2 := HI*COS2A-HR*SIN2A;
          HR := CX*CX*HR-SX*SX*W1-CA*T2;
          HI := CX*CX*HI+SX*SX*W2-SA*T2;
          B := SW*T1*CA+ETA*SA;
          E := CA*ETA-SW*T1*SA;

{     ROOT1 < EPS }

65:       S := HR-SIG*ROOT2*E;
          C := HI-SIG*ROOT2*B;
          ROOT := SQRT(C*C+S*S);

          IF ROOT < EPS THEN
          begin
            CB := 1.0;
            CH := 1.0;
            SB := 0.0;
            SH := 0.0;
            Goto 70
          end;
          CB := -C/ROOT;
          SB :=  S/ROOT;
          T2 := CB*B-E*SB;
          CN := T2*T2;
          TANH := ROOT/(G+2.0*(CN+DN));
          CH := 1.0/SQRT(1.0-TANH*TANH);
          SH := CH*TANH;

{     ROOT < EPS }

70:       W1 := SX*SH*(SA*CB-SB*CA);
          C1R := CX*CH-W1;
          C2R := CX*CH+W1;
          C1I :=-SX*SH*(CA*CB+SA*SB);
          C2I := C1I;
          W2 := SX*CH*CA;
          W1 := CX*SH*SB;
          S1R := W2-W1;
          S2R :=-W2-W1;
          W2 := SX*CH*SA;
          W1 := CX*SH*CB;
          S1I := W2+W1;
          S2I := W2-W1;
          W1 := SQRT(S1R*S1R+S1I*S1I);
          W2 := SQRT(S2R*S2R+S2I*S2I);

          IF (W1 > EPS) OR (W2 > EPS) THEN
          begin
{     BEGIN TRANSFORMATIONS }
            MARK := FALSE;
            For I := 1 to N do
            begin
              AKI := A^[K,I];
              AMI := A^[M,I];
              ZKI := Z^[K,I];
              ZMI := Z^[M,I];
              A^[K,I] := C1R*AKI-C1I*ZKI+S1R*AMI-S1I*ZMI;
              Z^[K,I] := C1R*ZKI+C1I*AKI+S1R*ZMI+S1I*AMI;
              A^[M,I] := S2R*AKI-S2I*ZKI+C2R*AMI-C2I*ZMI;
              Z^[M,I] := S2R*ZKI+S2I*AKI+C2R*ZMI+C2I*AMI
            end;
            For I := 1 to N do
            begin
              AIK := A^[I,K];
              AIM := A^[I,M];
              ZIK := Z^[I,K];
              ZIM := Z^[I,M];
              TIK := T^[I,K];
              TIM := T^[I,M];
              UIK := U^[I,K];
              UIM := U^[I,M];
              A^[I,K] := C2R*AIK-C2I*ZIK-S2R*AIM+S2I*ZIM;
              Z^[I,K] := C2R*ZIK+C2I*AIK-S2R*ZIM-S2I*AIM;
              A^[I,M] := C1R*AIM-C1I*ZIM-S1R*AIK+S1I*ZIK;
              Z^[I,M] := C1R*ZIM+C1I*AIM-S1R*ZIK-S1I*AIK;
              T^[I,K] := C2R*TIK-C2I*UIK-S2R*TIM+S2I*UIM;
              U^[I,K] := C2R*UIK+C2I*TIK-S2R*UIM-S2I*TIM;
              T^[I,M] := C1R*TIM-C1I*UIM-S1R*TIK+S1I*UIK;
              U^[I,M] := C1R*UIM+C1I*TIM-S1R*UIK-S1I*TIK
            end
          end {if w1>EPS...}

{     END TRANSFORMATIONS }

        End  {of M loop}
      End;   {of K loop}

{     GO TO NEXT ITERATION }

      Goto 20;

{     NO CONVERGENCE ! }

90:   IER := 1;
      Goto return;

{     CONVERGENCE OK }

95:   IER := 0;

{     SORT SOLUTIONS IN INCREASING ORDER }
      For J:=2 to N do
      begin
        VR:=A^[J,J];
        VI:=Z^[J,J];
        For K:=1 to N do
        begin
          EN^[K]:=T^[K,J];
          EN^[K+N]:=U^[K,J]
        end;
        For I:=J-1 Downto 1 do
        begin
          IF ABSCD(A^[I,I],Z^[I,I]) <=  ABSCD(VR,VI) Then Goto 97;
          A^[I+1,I+1]:=A^[I,I];
          Z^[I+1,I+1]:=Z^[I,I];
          For K:=1 to N do
          begin
            T^[K,I+1]:=T^[K,I];
            U^[K,I+1]:=U^[K,I]
          end
        end;
        I:=0;
97:     A^[I+1,I+1]:=VR;
        Z^[I+1,I+1]:=VI;
        For K:=1 to N do
        begin
          T^[K,I+1]:=EN^[K];
          U^[K,I+1]:=EN^[K+N]
        end
      end;

{     NORMALIZE VECTORS (BIGGEST COMPONENT TO UNITY) }
      For J:=1 to N do
      begin
        ZM:=0.0;
        For I:=1 to N do
        begin
          ZI:=ABS(T^[I,J]) + ABS(U^[I,J]);
          IF ZI >= ZM THEN
          begin
            IM:=I;
            ZM:=ZI
          end
        end;
        ZM:=T^[IM,J];
        ZI:=U^[IM,J];
        For I:=1 to N do
        begin
          DIVCD(T^[I,J],U^[I,J],ZM,ZI,TR,TI);
          T^[I,J]:=TR;
          U^[I,J]:=TI
        end
      end;
Return:End;

  {print square matrix of type MAT with caption}
  Procedure PrintMat(caption:string; N:Integer; A:pMAT);
  Var i,j:integer;
  Begin
    writeln(caption);
    For i:=1 to N do
    begin
      For j:=1 to N do write(' ',A^[i,j]:10:6);
      writeln
    end
  End;

{main program}
BEGIN

    New(AR); New(AI); New(T); New(U); New(EN);

    N:=5;   {size of given matrix}

    { Define matrix A := AR + I * AI }
    AR^[1,1]:= 1.0; AR^[1,2]:=  2.0; AR^[1,3]:= 3.0; AR^[1,4]:=-7.0;  AR^[1,5]:= 12.0;
    AR^[2,1]:= 2.0; AR^[2,2]:=  4.0; AR^[2,3]:= 7.0; AR^[2,4]:= 3.0;  AR^[2,5]:=-10.0;
    AR^[3,1]:= 3.0; AR^[3,2]:=  7.0; AR^[3,3]:=10.0; AR^[3,4]:= 8.0;  AR^[3,5]:=  4.0;
    AR^[4,1]:=-7.0; AR^[4,2]:=  3.0; AR^[4,3]:= 8.0; AR^[4,4]:=-0.75; AR^[4,5]:= -9.0;
    AR^[5,1]:=12.0; AR^[5,2]:=-10.0; AR^[5,3]:= 4.0; AR^[5,4]:=-9.0;  AR^[5,5]:= 10.0;

    For i:=1 to N do
      For j:=1 to N do
        AI^[i,j]:=0.0;   {here matrix is real, symmetric}

    writeln;
    PrintMat(' Matrix C (real parts):',N,AR);
    PrintMat(' Matrix C (imaginary parts):',N,AI);
    ReadKey;

    COMEIG (N,N,25,AR,AI,T,U,IER,EN);

    writeln;
    writeln(' Eigenvalues (real parts):');
    For i:=1 to N do write(' ',AR^[i,i]:10:6);
    writeln;
    writeln(' Eigenvalues (imaginary parts):');
    For i:=1 to N do write(' ',AI^[i,i]:10:6);
    writeln;
    writeln(' Error code (must be zero): ', IER);
    writeln;
    PrintMat(' Eigenvectors in columns (real parts):',N,T);
    PrintMat(' Eigenvectors in columns (imaginary parts):',N,U); 
    writeln;

    ReadKey;
    Dispose(AR); Dispose(AI); Dispose(T); Dispose(U); Dispose(EN);
    DoneWinCrt

END.

{end of file tcomeig.pas}